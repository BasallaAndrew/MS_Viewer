"""Main window: toolbar, three stacked panes, status bar, and selection state."""

from __future__ import annotations

import concurrent.futures
import traceback
from pathlib import Path
from typing import Callable

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets

from msview import index as ix
from msview.data import SUPPORTED_SUFFIXES, BackendUnavailable, RawSource, open_source
from msview.method_pane import MethodPane
from msview.panes import MS1_PEN, ChromatogramPane, MS2Pane, SpectrumPane

MAX_AVERAGE_SCANS = 500
AVERAGE_PPM = 5.0
MS2_SEARCH_RT_WINDOW = 1.0  # min either side, when clicking an MS1 peak
MS2_SEARCH_PPM = 20.0
XIC_DEFAULT_PPM = 10.0
XIC_DEFAULT_DA = 0.01


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("msview")
        self.resize(1200, 900)
        self.setAcceptDrops(True)

        self.source: RawSource | None = None
        self.idx: np.ndarray | None = None
        self.current_ms1: int | None = None  # set in single-scan mode
        self.rt_range: tuple[float, float] | None = None  # set in range mode
        self.ms2_list = np.empty(0, dtype=np.int32)
        self.ms2_pos = 0
        self._ms2_where = ""
        self._generation = 0  # invalidates stale background results
        self._xic_generation = 0
        self._method_state = "empty"  # empty / loading / loaded, for the current file
        self._xic_cache: dict[tuple, tuple] = {}  # (mz, tol_da) -> (rt, intensity, scans)
        # Python threads rather than QThreadPool: pythonnet (Thermo backend)
        # only works on threads Python created.
        self._executor = concurrent.futures.ThreadPoolExecutor(max_workers=2)
        # Finished futures are collected on the GUI thread by polling; emitting
        # Qt signals from executor threads can deadlock with PySide.
        self._futures: dict[concurrent.futures.Future, tuple] = {}
        self._poll = QtCore.QTimer(self, interval=25)
        self._poll.timeout.connect(self._collect_finished)

        self._build_ui()
        self._build_actions()

    # --- construction -----------------------------------------------------

    def _build_ui(self):
        self.chrom = ChromatogramPane()
        self.ms1 = SpectrumPane(MS1_PEN)
        self.ms2 = MS2Pane()
        self.ms1.set_spectrum(None, "MS1")
        self.ms2.set_spectrum(None, "MS2")
        self.ms2.set_nav_enabled(False, False)

        self.method = MethodPane()

        # Below the chromatogram: either the MS1/MS2 spectra or the Method view.
        spectra = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        spectra.setChildrenCollapsible(False)
        spectra.addWidget(self.ms1)
        spectra.addWidget(self.ms2)
        spectra.setSizes([350, 350])
        self.lower = QtWidgets.QStackedWidget()
        self.lower.addWidget(spectra)
        self.lower.addWidget(self.method)

        splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        splitter.setChildrenCollapsible(False)
        splitter.addWidget(self.chrom)
        splitter.addWidget(self.lower)
        splitter.setSizes([300, 700])
        self.setCentralWidget(splitter)
        for pane in (self.chrom, self.ms1, self.ms2, self.method.gradient):
            pane.hover.connect(self.statusBar().showMessage)
            pane.activated.connect(lambda p=pane: self.set_active_pane(p))
        self.active_pane = self.chrom
        self.chrom.set_active(True)
        self.setFocusPolicy(QtCore.Qt.FocusPolicy.StrongFocus)

        self.chrom.clicked.connect(self._on_chrom_clicked)
        self.chrom.rangeSelected.connect(self.select_range)
        self.ms1.mzClicked.connect(self._on_ms1_clicked)
        self.ms1.xicRequested.connect(self.show_xic)
        self.ms2.xicRequested.connect(self.show_xic)
        self.ms2.prevRequested.connect(lambda: self.step_ms2(-1))
        self.ms2.nextRequested.connect(lambda: self.step_ms2(1))

        self._status_style = self.statusBar().styleSheet()
        self.statusBar().showMessage("Open a Thermo .raw or .mzML file (⌘O or drag and drop)")
        self.chrom.title.setText("No file open: click Open (⌘O) or drag a .raw / .mzML file onto this window")

    def _build_actions(self):
        open_action = QtGui.QAction("Open…", self)
        open_action.setShortcut(QtGui.QKeySequence.StandardKey.Open)
        open_action.triggered.connect(self.choose_file)
        quit_action = QtGui.QAction("Quit", self)
        quit_action.setShortcut(QtGui.QKeySequence.StandardKey.Quit)
        quit_action.triggered.connect(self.close)

        file_menu = self.menuBar().addMenu("&File")
        file_menu.addAction(open_action)
        file_menu.addSeparator()
        file_menu.addAction(quit_action)

        toolbar = self.addToolBar("Main")
        toolbar.setMovable(False)
        toolbar.addAction(open_action)
        toolbar.addSeparator()
        self.kind_box = QtWidgets.QComboBox()
        self.kind_box.addItems(["TIC", "Base peak", "XIC"])
        self.kind_box.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        self.kind_box.currentIndexChanged.connect(self._on_kind_changed)
        toolbar.addWidget(self.kind_box)
        self.xic_action = toolbar.addWidget(self._build_xic_fields())
        self.xic_action.setVisible(False)
        toolbar.addSeparator()
        toolbar.addWidget(QtWidgets.QLabel(" View: "))
        self.view_box = QtWidgets.QComboBox()
        self.view_box.addItems(["MS1 / MS2 spectra", "Method & run info"])
        self.view_box.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        self.view_box.setToolTip("What to show below the chromatogram")
        self.view_box.currentIndexChanged.connect(self._on_view_changed)
        toolbar.addWidget(self.view_box)
        toolbar.addSeparator()
        self.file_label = QtWidgets.QLabel("No file")
        self.file_label.setContentsMargins(8, 0, 8, 0)
        toolbar.addWidget(self.file_label)


    def _build_xic_fields(self) -> QtWidgets.QWidget:
        box = QtWidgets.QWidget()
        box.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed)
        row = QtWidgets.QHBoxLayout(box)
        row.setContentsMargins(6, 0, 0, 0)
        row.setSpacing(4)
        self.xic_mz = QtWidgets.QLineEdit()
        self.xic_mz.setPlaceholderText("m/z")
        self.xic_mz.setToolTip("m/z to extract; press Return to update")
        self.xic_mz.setFixedWidth(100)
        validator = QtGui.QDoubleValidator(0.0, 1e6, 6)
        validator.setLocale(QtCore.QLocale.c())
        self.xic_mz.setValidator(validator)
        self.xic_tol = QtWidgets.QDoubleSpinBox()
        self.xic_tol.setDecimals(3)
        self.xic_tol.setRange(0.001, 1000.0)
        self.xic_tol.setValue(XIC_DEFAULT_PPM)
        self.xic_tol.setKeyboardTracking(False)
        self.xic_unit = QtWidgets.QComboBox()
        self.xic_unit.addItems(["ppm", "Da"])
        self.xic_unit.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        extract = QtWidgets.QPushButton("Extract")
        extract.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        for widget in (QtWidgets.QLabel("m/z"), self.xic_mz, QtWidgets.QLabel("±"),
                       self.xic_tol, self.xic_unit, extract):
            row.addWidget(widget)

        self.xic_mz.returnPressed.connect(self.extract_xic)
        self.xic_tol.valueChanged.connect(lambda _v: self.extract_xic(quiet=True))
        self.xic_unit.currentIndexChanged.connect(self._on_xic_unit_changed)
        extract.clicked.connect(self.extract_xic)
        return box

    # --- keyboard target ---------------------------------------------------

    def keyPressEvent(self, event):
        # Arrow keys reach here only when no text box / table has focus, so
        # typing an m/z or browsing the parameter table keeps its own arrows.
        handlers = {
            QtCore.Qt.Key.Key_Left: lambda: self.step_horizontal(-1),
            QtCore.Qt.Key.Key_Right: lambda: self.step_horizontal(1),
            QtCore.Qt.Key.Key_Up: lambda: self.step_ms2(-1),
            QtCore.Qt.Key.Key_Down: lambda: self.step_ms2(1),
        }
        handler = handlers.get(event.key())
        if handler is None:
            super().keyPressEvent(event)
            return
        event.accept()
        handler()

    def set_active_pane(self, pane):
        """The last pane clicked decides what ← / → do (highlighted header)."""
        self.setFocus()  # take the arrow keys back from any text box
        if pane is self.active_pane:
            return
        self.active_pane.set_active(False)
        self.active_pane = pane
        pane.set_active(True)

    def step_horizontal(self, step: int):
        """← / →: step MS2s when the MS2 pane is active, otherwise MS1 scans."""
        if self.active_pane is self.ms2:
            self.step_ms2(step)
        else:
            self.step_ms1(step)

    # --- method view ---------------------------------------------------------

    def _method_view(self) -> bool:
        return self.view_box.currentIndex() == 1

    def _on_view_changed(self, index: int):
        self.lower.setCurrentIndex(index)
        if not self._method_view() and self.active_pane is self.method.gradient:
            self.set_active_pane(self.chrom)
        self._load_method_info()

    def _load_method_info(self):
        """Read method/run info in the background the first time it's needed."""
        if not self._method_view() or self.source is None or self._method_state != "empty":
            return
        self._method_state = "loading"
        self.method.set_info(None, "Reading method information…")
        source = self.source

        def done(info):
            if source is self.source:
                self._method_state = "loaded"
                self.method.set_info(info)
                self._sync_method_cursor()

        def failed(exc):
            if source is self.source:
                self._method_state = "loaded"
                self.method.set_info(None, f"Could not read method information: {exc}")

        self._run(source.method_info, done, failed)

    def _sync_method_cursor(self):
        if self.current_ms1 is not None:
            self.method.gradient.show_cursor(float(self.idx["rt"][ix.row_of(self.idx, self.current_ms1)]))
        elif self.rt_range is not None:
            self.method.gradient.show_range(*self.rt_range)

    # --- background work --------------------------------------------------

    def _run(self, fn: Callable[[], object], on_done, on_failed=None):
        self._futures[self._executor.submit(fn)] = (on_done, on_failed or self._report_error)
        self._poll.start()

    def _collect_finished(self):
        for future in [f for f in self._futures if f.done()]:
            on_done, on_failed = self._futures.pop(future)
            exc = future.exception()
            if exc is None:
                on_done(future.result())
                continue
            try:
                exc.formatted = "".join(traceback.format_exception(exc))
            except Exception:
                pass
            on_failed(exc)
        if not self._futures:
            self._poll.stop()

    def _wait_for_tasks(self):
        concurrent.futures.wait(list(self._futures))

    def _report_error(self, exc: Exception):
        self.statusBar().showMessage(f"Error: {exc}")
        print(getattr(exc, "formatted", exc))

    def flash_status(self, message: str):
        self.statusBar().showMessage(message, 4000)
        self.statusBar().setStyleSheet("QStatusBar { background: #f8d7da; }")
        QtCore.QTimer.singleShot(600, lambda: self.statusBar().setStyleSheet(self._status_style))

    # --- opening files ----------------------------------------------------

    def choose_file(self):
        start = QtCore.QSettings().value("last_dir", str(Path.home()))
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open LC-MS file", start,
            "MS data (*.raw *.RAW *.mzML *.mzml);;All files (*)",
        )
        if path:
            self.open_file(path)

    def open_file(self, path: str | Path):
        path = Path(path)
        QtCore.QSettings().setValue("last_dir", str(path.parent))
        self._generation += 1
        generation = self._generation
        # Let in-flight work on the old file finish before closing it.
        self._wait_for_tasks()
        if self.source is not None:
            self.source.close()
        self.source, self.idx = None, None
        self._xic_cache.clear()
        self._method_state = "empty"
        self.method.set_info(None, f"Loading {path.name}…")
        self.current_ms1, self.rt_range = None, None
        self.ms2_list = np.empty(0, dtype=np.int32)
        self.chrom.clear()
        self.ms1.set_spectrum(None, "MS1")
        self.ms1.set_precursor(None, None)
        self.ms2.set_spectrum(None, "MS2")
        self.ms2.set_nav_enabled(False, False)
        self.file_label.setText(f"{path.name} (loading…)")
        self.setWindowTitle(f"msview: {path.name}")
        self.statusBar().showMessage(f"Indexing {path.name}…")
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)

        if path.suffix.lower() == ".raw":
            # Start .NET on the GUI thread; pythonnet is fragile when the
            # runtime is first loaded from a worker thread.
            try:
                from msview.data.thermo import ensure_runtime

                ensure_runtime()
            except BackendUnavailable as exc:
                QtWidgets.QApplication.restoreOverrideCursor()
                self._show_open_error(path, exc)
                return

        def load():
            source = open_source(path)
            try:
                return source, source.scan_index()
            except Exception:
                source.close()
                raise

        def loaded(result):
            QtWidgets.QApplication.restoreOverrideCursor()
            source, idx = result
            if generation != self._generation:
                source.close()
                return
            self._on_loaded(path, source, idx)

        def failed(exc):
            QtWidgets.QApplication.restoreOverrideCursor()
            if generation == self._generation:
                self._show_open_error(path, exc)

        self._run(load, loaded, failed)

    def _show_open_error(self, path: Path, exc: Exception):
        self.file_label.setText("No file")
        self.setWindowTitle("msview")
        self.statusBar().showMessage(f"Could not open {path.name}")
        box = QtWidgets.QMessageBox(self)
        box.setIcon(QtWidgets.QMessageBox.Icon.Critical)
        box.setWindowTitle("Could not open file")
        box.setText(str(exc))
        if not isinstance(exc, (BackendUnavailable, FileNotFoundError, ValueError)):
            box.setDetailedText(getattr(exc, "formatted", ""))
        box.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        box.exec()

    def _on_loaded(self, path: Path, source: RawSource, idx: np.ndarray):
        self.source, self.idx = source, idx
        n_ms1 = int((idx["ms_level"] == 1).sum())
        self.file_label.setText(
            f"{path.name}   |   {len(idx):,} scans ({n_ms1:,} MS1, {len(idx) - n_ms1:,} MSn)   |   "
            f"RT {idx['rt'][0]:.1f}–{idx['rt'][-1]:.1f} min"
        )
        self.statusBar().showMessage(f"Loaded {path.name}", 3000)
        self._update_chromatogram(keep_view=False)
        self._load_method_info()
        if n_ms1:
            # Start at the most intense point of the run.
            rt, tic, scans = source.chromatogram("tic")
            self.select_ms1(int(scans[int(np.argmax(tic))]))

    def _xic_mode(self) -> bool:
        return self.kind_box.currentIndex() == 2

    def _on_kind_changed(self, *_):
        xic = self._xic_mode()
        self.xic_action.setVisible(xic)
        if xic and not self.xic_mz.text():
            self._prefill_xic_mz()
        self._update_chromatogram()
        if xic and not self.xic_mz.hasAcceptableInput():
            self.xic_mz.setFocus()

    def _prefill_xic_mz(self):
        """Start from the current MS2 precursor, if there is one."""
        if self.idx is not None and len(self.ms2_list):
            row = self.idx[ix.row_of(self.idx, int(self.ms2_list[self.ms2_pos]))]
            if np.isfinite(row["precursor_mz"]):
                self.xic_mz.setText(f"{row['precursor_mz']:.4f}")

    def _update_chromatogram(self, *_, keep_view: bool = True):
        if self.source is None:
            return
        if self._xic_mode():
            if self.xic_mz.hasAcceptableInput():
                self.extract_xic(keep_view=keep_view)
            else:
                self.statusBar().showMessage("Type an m/z in the XIC box and press Return")
            return
        kind = "tic" if self.kind_box.currentIndex() == 0 else "bpc"
        rt, inten, scans = self.source.chromatogram(kind)
        self.chrom.set_data(rt, inten, scans, "TIC" if kind == "tic" else "Base peak", keep_view=keep_view)

    def _on_xic_unit_changed(self, *_):
        default = XIC_DEFAULT_PPM if self.xic_unit.currentText() == "ppm" else XIC_DEFAULT_DA
        self.xic_tol.setValue(default)  # triggers a re-extract

    def show_xic(self, mz: float):
        """Switch the chromatogram to an XIC of mz (used by the spectrum context menu)."""
        self.xic_mz.setText(f"{mz:.4f}")
        if self._xic_mode():
            self.extract_xic()
        else:
            self.kind_box.setCurrentIndex(2)  # _on_kind_changed extracts

    def extract_xic(self, *_, quiet: bool = False, keep_view: bool = True):
        if self.source is None or not self._xic_mode():
            return
        if not self.xic_mz.hasAcceptableInput():
            if not quiet:
                self.flash_status("Enter a valid m/z for the XIC")
            return
        self.xic_mz.clearFocus()  # hand the arrow keys back to scan navigation
        mz = float(self.xic_mz.text())
        tol = self.xic_tol.value()
        ppm = self.xic_unit.currentText() == "ppm"
        tol_da = mz * tol * 1e-6 if ppm else tol
        label = f"{tol:g} ppm" if ppm else f"{tol:g} Da"
        title = f"XIC   |   m/z {mz:.4f} ± {label}   |   MS1"
        export_name = f"xic_{mz:.4f}.csv"

        def show(result):
            rt, inten, scans = result
            self.chrom.set_data(rt, inten, scans, "XIC", keep_view=keep_view,
                                title=title, export_name=export_name)

        key = (round(mz, 6), round(tol_da, 9))
        self._xic_generation += 1
        generation = self._xic_generation
        if key in self._xic_cache:
            show(self._xic_cache[key])
            return
        self.chrom.title.setText(f"Extracting XIC for m/z {mz:.4f} ± {label}…")
        source = self.source

        def done(result):
            if source is self.source:
                self._xic_cache[key] = result
            if generation == self._xic_generation and source is self.source and self._xic_mode():
                show(result)

        self._run(lambda: source.xic(mz, tol_da), done)

    # --- drag and drop ----------------------------------------------------

    @staticmethod
    def _dropped_path(event) -> str | None:
        urls = event.mimeData().urls()
        if urls and urls[0].isLocalFile():
            path = urls[0].toLocalFile()
            if Path(path).suffix.lower() in SUPPORTED_SUFFIXES:
                return path
        return None

    def dragEnterEvent(self, event):
        if self._dropped_path(event):
            event.acceptProposedAction()

    def dropEvent(self, event):
        path = self._dropped_path(event)
        if path:
            event.acceptProposedAction()
            self.open_file(path)

    # --- selection --------------------------------------------------------

    def _on_chrom_clicked(self, rt: float):
        if self.idx is None:
            return
        scan = ix.nearest_ms1(self.idx, rt)
        if scan is not None:
            self.select_ms1(scan)

    def select_ms1(self, scan: int):
        """Single-scan mode: show one MS1 and the MS2s of its duty cycle."""
        self._generation += 1  # cancel any pending average
        row = self.idx[ix.row_of(self.idx, scan)]
        self.current_ms1, self.rt_range = scan, None
        self.chrom.show_cursor(float(row["rt"]))
        self.method.gradient.show_cursor(float(row["rt"]))
        spec = self.source.spectrum(scan)
        self.ms1.set_spectrum(spec, f"MS1   |   scan {scan}   |   RT {row['rt']:.3f} min",
                              export_name=f"ms1_scan{scan}.csv")
        self._set_ms2_list(ix.ms2_for_cycle(self.idx, scan), "in this cycle")

    def select_range(self, rt0: float, rt1: float):
        """Range mode: average the MS1s in [rt0, rt1] and list every MS2 in it."""
        if self.idx is None:
            return
        ms1_scans = ix.ms1_in_range(self.idx, rt0, rt1)
        if not len(ms1_scans):
            scan = ix.nearest_ms1(self.idx, (rt0 + rt1) / 2)
            if scan is not None:
                self.select_ms1(scan)
            return
        self.current_ms1, self.rt_range = None, (rt0, rt1)
        self.chrom.show_range(rt0, rt1)
        self.method.gradient.show_range(rt0, rt1)
        self._set_ms2_list(ix.ms2_in_range(self.idx, rt0, rt1), "in this range")

        used = ix.evenly_subsample(ms1_scans, MAX_AVERAGE_SCANS)
        warning = ""
        if len(used) < len(ms1_scans):
            warning = f" (capped: {len(used)} of {len(ms1_scans)} scans, evenly spaced)"
            self.flash_status(f"Range has {len(ms1_scans)} MS1 scans; averaging {len(used)} evenly spaced scans")
        self._generation += 1
        generation = self._generation
        source = self.source
        self.ms1.title.setText(f"Averaging {len(used)} MS1 scans…")

        def done(spec):
            if generation != self._generation:
                return
            title = (f"MS1   |   average of {len(used)} scans   |   RT {rt0:.3f}–{rt1:.3f} min{warning}")
            self.ms1.set_spectrum(spec, title, export_name=f"ms1_avg_{rt0:.2f}-{rt1:.2f}min.csv")

        self._run(lambda: source.average(used, ppm=AVERAGE_PPM), done)

    def step_ms1(self, step: int):
        if self.idx is None:
            return
        if self.current_ms1 is not None:
            origin = self.current_ms1
        elif self.rt_range is not None:
            origin = ix.nearest_ms1(self.idx, self.rt_range[0] if step < 0 else self.rt_range[1])
        else:
            return
        scan = ix.step_ms1(self.idx, origin, step)
        if scan is not None and (scan != self.current_ms1):
            self.select_ms1(scan)

    # --- MS2 --------------------------------------------------------------

    def _set_ms2_list(self, scans: np.ndarray, where: str):
        self.ms2_list = scans
        self.ms2_pos = 0
        self._ms2_where = where
        self._show_ms2()

    def step_ms2(self, step: int):
        if not len(self.ms2_list):
            return
        pos = int(np.clip(self.ms2_pos + step, 0, len(self.ms2_list) - 1))
        if pos != self.ms2_pos:
            self.ms2_pos = pos
            self._show_ms2()

    def _show_ms2(self):
        n = len(self.ms2_list)
        self.ms2.set_nav_enabled(self.ms2_pos > 0, self.ms2_pos < n - 1)
        if not n:
            self.ms2.set_spectrum(None, f"No MS2 scans {self._ms2_where}")
            self.ms1.set_precursor(None, None)
            self.chrom.show_ms2_marker(None)
            return
        scan = int(self.ms2_list[self.ms2_pos])
        row = self.idx[ix.row_of(self.idx, scan)]
        self.ms2.set_spectrum(self.source.spectrum(scan), self._ms2_title(row, self.ms2_pos, n),
                              export_name=f"ms{row['ms_level']}_scan{scan}.csv")
        self.ms1.set_precursor(float(row["precursor_mz"]), ix.isolation_window(row))
        self.chrom.show_ms2_marker(float(row["rt"]))

    @staticmethod
    def _ms2_title(row, pos: int, n: int) -> str:
        parts = [f"MS{row['ms_level']} {pos + 1} of {n}", f"scan {row['scan']}", f"RT {row['rt']:.3f}"]
        if np.isfinite(row["precursor_mz"]):
            charge = f"z={row['charge']}" if row["charge"] else "z=?"
            parts.append(f"precursor {row['precursor_mz']:.4f} ({charge})")
        activation = str(row["activation"])
        if np.isfinite(row["collision_energy"]):
            activation = f"{activation} {row['collision_energy']:g}".strip()
        if activation:
            parts.append(activation)
        return "   |   ".join(parts)

    def _on_ms1_clicked(self, x: float):
        """Find the MS2 that fragmented the clicked MS1 peak.

        First look in the current MS2 list; otherwise search the run within
        ±MS2_SEARCH_RT_WINDOW of the MS1 being shown, and make the MS2 list
        every MS2 of that precursor so the arrows step through them.
        """
        if self.idx is None:
            return
        mz = self.ms1.peak_near(x, px=10)
        if mz is None:
            self.flash_status("Click on a peak to find its MS2")
            return

        if len(self.ms2_list):
            rows = np.searchsorted(self.idx["scan"], self.ms2_list)
            exact, in_window = ix.precursor_matches(self.idx, rows, mz, MS2_SEARCH_PPM)
            hits = exact if exact.any() else in_window
            if hits.any():
                distance = np.where(hits, np.abs(self.idx["precursor_mz"][rows] - mz), np.inf)
                self.ms2_pos = int(np.argmin(distance))
                self._show_ms2()
                return

        if self.current_ms1 is not None:
            rt = float(self.idx["rt"][ix.row_of(self.idx, self.current_ms1)])
            window = MS2_SEARCH_RT_WINDOW
        else:
            rt0, rt1 = self.rt_range
            rt, window = (rt0 + rt1) / 2, (rt1 - rt0) / 2 + MS2_SEARCH_RT_WINDOW
        scans, pos = ix.ms2_for_precursor(self.idx, mz, rt, window, MS2_SEARCH_PPM)
        if pos < 0:
            self.flash_status(f"No MS2 of m/z {mz:.4f} within ±{MS2_SEARCH_RT_WINDOW:g} min")
            return
        self.ms2_list, self.ms2_pos = scans, pos
        self._ms2_where = f"for m/z {mz:.4f}"
        self._show_ms2()
        found_rt = float(self.idx["rt"][ix.row_of(self.idx, int(scans[pos]))])
        self.statusBar().showMessage(
            f"{len(scans)} MS2 scan(s) of m/z {mz:.4f} near RT {rt:.2f}; showing the closest "
            f"(RT {found_rt:.3f}). Use ↑/↓ to step through them.", 8000)

    def closeEvent(self, event):
        self._wait_for_tasks()
        if self.source is not None:
            self.source.close()
        super().closeEvent(event)
