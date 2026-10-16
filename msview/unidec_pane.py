"""UniDec view: the selected MS1 (top) and its deconvolved mass spectrum (bottom).

The input is whatever the MS1 pane shows: a single scan or a dragged average.
"""

from __future__ import annotations

import sys
from typing import Callable

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

from msview.data.base import Spectrum
from msview.deconv import DeconResult, DeconSettings, charge_ladder, deconvolve
from msview.panes import MS1_PEN, SpectrumPane

MASS_PEN = pg.mkPen("#1b7a3e", width=1.2)
FIT_PEN = pg.mkPen("#d62728", width=1, style=QtCore.Qt.PenStyle.DashLine)
PEAK_COLOURS = ["#d62728", "#9467bd", "#ff7f0e", "#17becf", "#8c564b"]
MAX_MARKED_PEAKS = len(PEAK_COLOURS)
LADDER_MIN_FRACTION = 0.02  # only mark charge states with at least this much signal


def _spin(lo, hi, value, decimals=0, step=1.0, suffix="", width=None):
    box = QtWidgets.QDoubleSpinBox() if decimals else QtWidgets.QSpinBox()
    if decimals:
        box.setDecimals(decimals)
    box.setRange(lo, hi)
    box.setSingleStep(step)
    box.setValue(value)
    box.setKeyboardTracking(False)  # update on Return / focus-out, not every keystroke
    box.setGroupSeparatorShown(True)
    if suffix:
        box.setSuffix(suffix)
    if width:
        box.setFixedWidth(width)
    return box


class UniDecPane(QtWidgets.QWidget):
    def __init__(self, runner: Callable, parent=None):
        """runner(fn, on_done, on_failed) runs fn in the background (MainWindow._run)."""
        super().__init__(parent)
        self._runner = runner
        self._input: Spectrum | None = None
        self._input_label = ""
        self._result: DeconResult | None = None
        self._shown = False
        self._busy = False
        self._pending = False
        self._stale = False  # input/settings changed since the last result

        defaults = DeconSettings()
        self.z_min = _spin(1, 500, defaults.z_min, width=60)
        self.z_max = _spin(1, 500, defaults.z_max, width=60)
        self.mass_min = _spin(1, 10_000_000, defaults.mass_min, decimals=1, step=1000, width=100)
        self.mass_max = _spin(1, 10_000_000, defaults.mass_max, decimals=1, step=1000, width=100)
        self.mass_bin = _spin(0.01, 1000, defaults.mass_bin, decimals=2, step=0.5, width=70)
        self.auto_fwhm = QtWidgets.QCheckBox("auto")
        self.auto_fwhm.setChecked(True)
        self.fwhm = _spin(0.001, 100, 0.85, decimals=3, step=0.05, width=75)
        self.fwhm.setEnabled(False)
        self.peak_window = _spin(0.1, 100_000, defaults.peak_window, decimals=1, step=5, width=75)
        self.threshold = _spin(0, 100, defaults.peak_threshold * 100, decimals=1, step=1, suffix=" %", width=75)
        self.zoom_only = QtWidgets.QCheckBox("Only zoomed m/z range")
        self.zoom_only.setToolTip("Deconvolve just the m/z range visible in the top pane")
        self.auto_run = QtWidgets.QCheckBox("Auto-run")
        self.auto_run.setChecked(True)
        self.auto_run.setToolTip("Deconvolve whenever the MS1 selection or a setting changes")
        self.run_button = QtWidgets.QPushButton("Deconvolve")
        self.run_button.clicked.connect(self.request_run)
        for w in (self.auto_fwhm, self.zoom_only, self.auto_run, self.run_button):
            w.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)

        tips = {
            self.z_min: "Lowest charge state to consider", self.z_max: "Highest charge state to consider",
            self.mass_min: "Lowest mass to consider (Da)", self.mass_max: "Highest mass to consider (Da)",
            self.mass_bin: "Mass resolution of the deconvolved spectrum (Da)",
            self.fwhm: "Peak width in m/z used by UniDec (auto estimates it from the data)",
            self.peak_window: "Mass peaks closer than this are merged when picking (Da)",
            self.threshold: "Minimum peak height to pick, as % of the tallest",
        }
        for w, tip in tips.items():
            w.setToolTip(tip)

        controls = QtWidgets.QHBoxLayout()
        controls.setContentsMargins(8, 4, 8, 4)
        controls.setSpacing(6)
        for item in (QtWidgets.QLabel("Charge"), self.z_min, QtWidgets.QLabel("–"), self.z_max,
                     QtWidgets.QLabel("  Mass"), self.mass_min, QtWidgets.QLabel("–"), self.mass_max,
                     QtWidgets.QLabel("Da   Bin"), self.mass_bin, QtWidgets.QLabel("Da   Peak FWHM"),
                     self.auto_fwhm, self.fwhm, QtWidgets.QLabel("m/z   Peak window"), self.peak_window,
                     QtWidgets.QLabel("Da   Threshold"), self.threshold):
            controls.addWidget(item)
        controls.addStretch(1)
        controls2 = QtWidgets.QHBoxLayout()
        controls2.setContentsMargins(8, 0, 8, 4)
        controls2.setSpacing(10)
        self.status = QtWidgets.QLabel("")
        self.status.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        for item in (self.zoom_only, self.auto_run, self.run_button):
            controls2.addWidget(item)
        controls2.addWidget(self.status, 1)

        self.spectrum = SpectrumPane(MS1_PEN)
        self.masses = SpectrumPane(MASS_PEN, x_unit="Da", label_decimals=1)
        self.masses.plot.getPlotItem().setLabel("left", "Relative intensity (%)")

        self.fit_curve = pg.PlotCurveItem(pen=FIT_PEN)
        self.fit_curve.setZValue(5)
        self.spectrum.plot.addItem(self.fit_curve)
        self.ladders = []
        for colour in PEAK_COLOURS:
            ladder = pg.ScatterPlotItem(symbol="t", size=9, pen=None, brush=pg.mkBrush(colour))
            ladder.setZValue(6)
            self.spectrum.plot.addItem(ladder, ignoreBounds=True)
            peak_marker = pg.ScatterPlotItem(symbol="t", size=11, pen=None, brush=pg.mkBrush(colour))
            peak_marker.setZValue(6)
            self.masses.plot.addItem(peak_marker, ignoreBounds=True)
            self.ladders.append((ladder, peak_marker))
        self.z_labels = []  # charge labels for the most intense mass

        splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        splitter.setChildrenCollapsible(False)
        splitter.addWidget(self.spectrum)
        splitter.addWidget(self.masses)
        splitter.setSizes([400, 400])

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addLayout(controls)
        layout.addLayout(controls2)
        layout.addWidget(splitter, 1)

        self.auto_fwhm.toggled.connect(lambda on: self.fwhm.setEnabled(not on))
        for box in (self.z_min, self.z_max, self.mass_min, self.mass_max, self.mass_bin, self.fwhm,
                    self.peak_window, self.threshold):
            box.valueChanged.connect(self._settings_changed)
        for check in (self.auto_fwhm, self.zoom_only):
            check.toggled.connect(self._settings_changed)
        self.auto_run.toggled.connect(lambda on: on and self._settings_changed())
        self.clear("Select an MS1 scan or drag a range on the chromatogram")

    # --- input ------------------------------------------------------------------

    def clear(self, message: str = "") -> None:
        self._input, self._result, self._stale = None, None, False
        self.spectrum.set_spectrum(None, "Input spectrum")
        self._show_result(None, message or "Deconvolved mass spectrum")
        self.status.setText(message)

    def set_input(self, spec: Spectrum | None, label: str) -> None:
        """Called whenever the MS1 selection changes."""
        self._input, self._input_label = spec, label
        self.spectrum.set_spectrum(spec, f"Input: {label}", export_name="unidec_input.csv")
        self._mark_stale()

    def set_shown(self, shown: bool) -> None:
        """Only deconvolve while the UniDec view is on screen."""
        self._shown = shown
        if shown and self._stale and self.auto_run.isChecked():
            self.request_run()

    def _settings_changed(self, *_):
        if self._input is not None:
            self._mark_stale()

    def _mark_stale(self) -> None:
        self._stale = True
        self._show_result(None, "Deconvolved mass spectrum: "
                          + ("running…" if self.auto_run.isChecked() and self._shown else "click Deconvolve"))
        if self.auto_run.isChecked() and self._shown:
            self.request_run()

    def settings(self) -> DeconSettings:
        mz_min = mz_max = None
        if self.zoom_only.isChecked() and self._input is not None and len(self._input):
            mz_min, mz_max = self.spectrum.vb.viewRange()[0]
        return DeconSettings(
            z_min=self.z_min.value(), z_max=self.z_max.value(),
            mass_min=self.mass_min.value(), mass_max=self.mass_max.value(), mass_bin=self.mass_bin.value(),
            peak_fwhm=None if self.auto_fwhm.isChecked() else self.fwhm.value(),
            mz_min=mz_min, mz_max=mz_max,
            peak_window=self.peak_window.value(), peak_threshold=self.threshold.value() / 100,
        )

    # --- running -------------------------------------------------------------------

    def request_run(self) -> None:
        if self._input is None or not len(self._input):
            self.status.setText("Nothing to deconvolve yet: select an MS1 scan or range")
            return
        if self._busy:  # one run at a time; rerun with the latest input when it finishes
            self._pending = True
            return
        self._busy, self._pending, self._stale = True, False, False
        spec, settings = self._input, self.settings()
        first_use = "unidec.engine" not in sys.modules
        self.status.setText("Loading UniDec (the first run can take ~30 s)…" if first_use else "Deconvolving…")
        self.run_button.setEnabled(False)
        self._runner(lambda: deconvolve(spec.mz, spec.intensity, settings), self._on_done, self._on_failed)

    def _finish_run(self) -> bool:
        """Common bookkeeping; returns True if the result is still current."""
        self._busy = False
        self.run_button.setEnabled(True)
        if self._pending:
            self.request_run()
            return False
        return True

    def _on_done(self, result: DeconResult) -> None:
        if not self._finish_run():
            return
        self._result = result
        self._show_result(result)
        self.status.setText(
            f"UniDec: {result.seconds:.2f} s   ·   R² {result.r_squared:.3f}   ·   "
            f"peak FWHM {result.peak_fwhm:.3g} m/z   ·   {len(result.peaks)} peak(s)")

    def _on_failed(self, exc: Exception) -> None:
        if not self._finish_run():
            return
        self._show_result(None, "Deconvolution failed")
        self.status.setText(f'<span style="color:#b00020">{exc}</span>')
        print(getattr(exc, "formatted", exc))

    # --- display ------------------------------------------------------------------

    def _show_result(self, result: DeconResult | None, message: str = "") -> None:
        for ladder, peak_marker in self.ladders:
            ladder.setData([], [])
            peak_marker.setData([], [])
        for label in self.z_labels:
            self.spectrum.plot.removeItem(label)
        self.z_labels = []
        if result is None:
            self.fit_curve.setData([], [])
            self.masses.set_spectrum(None, message or "Deconvolved mass spectrum", export_name="unidec_masses.csv")
            return

        self.masses.label_min_fraction = self.threshold.value() / 100
        self.masses.set_spectrum(Spectrum(result.mass, result.intensity, centroid=False),
                                 self._masses_title(result), export_name="unidec_masses.csv")

        # Overlay UniDec's fit on the input, scaled to the input's intensity.
        spec = self._input
        lo, hi = float(result.data_mz[0]), float(result.data_mz[-1])
        window = (spec.mz >= lo) & (spec.mz <= hi)
        top_in = float(spec.intensity[window].max()) if window.any() else 1.0
        top_data = float(np.nanmax(result.data_intensity)) or 1.0
        self.fit_curve.setData(result.data_mz, np.nan_to_num(result.fit) * top_in / top_data)

        # Charge-state ladders for the main masses, placed on the data.
        settings = self.settings()
        for k, peak in enumerate(result.peaks[:MAX_MARKED_PEAKS]):
            ladder, peak_marker = self.ladders[k]
            peak_marker.setData([peak.mass], [peak.height])
            z, mz = charge_ladder(peak.mass, settings.z_min, settings.z_max, lo, hi)
            if not len(mz):
                continue
            heights = self._local_max(spec, mz, max(result.peak_fwhm, 0.01))
            real = heights >= LADDER_MIN_FRACTION * top_in
            z, mz, heights = z[real], mz[real], heights[real]
            ladder.setData(mz, heights)
            if k == 0:
                for zi, x, y in zip(z, mz, heights):
                    label = pg.TextItem(f"{zi}+", color=PEAK_COLOURS[0], anchor=(0.5, 1.6))
                    label.setPos(x, y)
                    self.spectrum.plot.addItem(label, ignoreBounds=True)
                    self.z_labels.append(label)

    @staticmethod
    def _local_max(spec: Spectrum, centres: np.ndarray, width: float) -> np.ndarray:
        lo = np.searchsorted(spec.mz, centres - width / 2)
        hi = np.searchsorted(spec.mz, centres + width / 2, side="right")
        return np.array([spec.intensity[a:b].max() if b > a else 0.0 for a, b in zip(lo, hi)])

    @staticmethod
    def _masses_title(result: DeconResult) -> str:
        if not result.peaks:
            return "Deconvolved mass spectrum: no peaks above the threshold"
        parts = [f'<span style="color:{PEAK_COLOURS[k]}">&#9650;</span>&nbsp;{p.mass:,.1f} Da ({p.height:.0f}%)'
                 for k, p in enumerate(result.peaks[:MAX_MARKED_PEAKS])]
        more = f"&nbsp;&nbsp;+{len(result.peaks) - MAX_MARKED_PEAKS} more" if len(result.peaks) > MAX_MARKED_PEAKS else ""
        return "Masses:&nbsp;&nbsp;" + "&nbsp;&nbsp;&nbsp;".join(parts) + more
