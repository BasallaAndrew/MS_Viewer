"""Plot panes: chromatogram, MS1 spectrum and MS2 spectrum.

Mouse conventions shared by every pane (see NavViewBox):
  scroll = zoom x, shift+scroll = zoom y, double-click = reset view,
  right-click = context menu (copy image / export CSV).
Chromatogram: left-drag selects an RT range, right/middle-drag pans.
Spectra: any drag pans m/z; the y-axis follows the tallest visible peak.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets

from msview.data.base import Spectrum

pg.setConfigOptions(background="w", foreground="#222", antialias=False)

AXIS_WIDTH = 64
Y_HEADROOM = 1.15  # leave room above the tallest peak for its label
LABEL_COUNT = 10
LABEL_MIN_PX = 58  # minimum horizontal spacing between peak labels

CHROM_PEN = pg.mkPen("#333333", width=1)
MS1_PEN = pg.mkPen("#1f4e9c", width=1)
MS2_PEN = pg.mkPen("#8c1c13", width=1)
CURSOR_PEN = pg.mkPen("#d62728", width=2)
PRECURSOR_PEN = pg.mkPen("#e67e00", width=1.5, style=QtCore.Qt.PenStyle.DashLine)
MS2_MARKER_PEN = pg.mkPen("#e67e00", width=1.5, style=QtCore.Qt.PenStyle.DashLine)
RANGE_BRUSH = pg.mkBrush(31, 119, 180, 50)
ACTIVE_HEADER_STYLE = "#headerBar { background: #d6e4f5; border-top: 2px solid #1f4e9c; }"
INACTIVE_HEADER_STYLE = "#headerBar { border-top: 2px solid transparent; }"
ISOLATION_BRUSH = pg.mkBrush(230, 126, 0, 45)


class NavViewBox(pg.ViewBox):
    sigClicked = QtCore.Signal(float, float)
    sigRangeDrag = QtCore.Signal(float, float, bool)  # x0, x1, finished
    sigContextMenu = QtCore.Signal(QtCore.QPoint, float)  # screen pos, data x
    sigUserZoomed = QtCore.Signal()

    def __init__(self, select_drag: bool = False):
        super().__init__(enableMenu=False)
        self.select_drag = select_drag
        self.setMouseEnabled(x=True, y=True)

    def wheelEvent(self, ev, axis=None):
        delta = ev.delta()
        if not delta:
            ev.ignore()
            return
        s = 1.02 ** (-delta / 8)
        center = self.mapToView(ev.pos())
        if ev.modifiers() & QtCore.Qt.KeyboardModifier.ShiftModifier or axis == 1:
            # Zoom intensity about the baseline.
            self.scaleBy(y=s, center=QtCore.QPointF(center.x(), 0))
        else:
            self.scaleBy(x=s, center=center)
            self.sigUserZoomed.emit()
        ev.accept()

    def mouseDragEvent(self, ev, axis=None):
        if axis is not None:  # dragging on an axis: default behaviour
            super().mouseDragEvent(ev, axis)
            return
        ev.accept()
        if self.select_drag and ev.button() == QtCore.Qt.MouseButton.LeftButton:
            x0 = self.mapToView(ev.buttonDownPos()).x()
            x1 = self.mapToView(ev.pos()).x()
            self.sigRangeDrag.emit(min(x0, x1), max(x0, x1), ev.isFinish())
            return
        dx = self.mapToView(ev.lastPos()).x() - self.mapToView(ev.pos()).x()
        self.translateBy(x=dx)
        self.sigUserZoomed.emit()

    def mouseClickEvent(self, ev):
        button = ev.button()
        if button == QtCore.Qt.MouseButton.LeftButton:
            ev.accept()
            if not ev.double():  # double-clicks are handled by PlotPane.eventFilter
                p = self.mapToView(ev.pos())
                self.sigClicked.emit(p.x(), p.y())
        elif button == QtCore.Qt.MouseButton.RightButton:
            ev.accept()
            self.sigContextMenu.emit(ev.screenPos().toPoint(), self.mapToView(ev.pos()).x())


class PlotPane(QtWidgets.QWidget):
    """A header row above a pyqtgraph plot, with shared zoom/export behaviour."""

    hover = QtCore.Signal(str)
    activated = QtCore.Signal()  # any mouse press in the plot
    export_header = "x,y"
    export_name = "data.csv"

    def __init__(self, x_label: str, select_drag: bool = False, parent=None):
        super().__init__(parent)
        self.vb = NavViewBox(select_drag=select_drag)
        self.plot = pg.PlotWidget(viewBox=self.vb)
        self.plot.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        pi = self.plot.getPlotItem()
        pi.hideButtons()
        pi.setLabel("bottom", x_label)
        pi.getAxis("left").setWidth(AXIS_WIDTH)
        pi.getAxis("left").enableAutoSIPrefix(False)
        self.vb.disableAutoRange()

        # The header bar is highlighted when this pane is the keyboard target.
        self.header_bar = QtWidgets.QWidget(objectName="headerBar")
        self.header_bar.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground, True)
        self.header = QtWidgets.QHBoxLayout(self.header_bar)
        self.header.setContentsMargins(8, 2, 8, 2)
        self.title = QtWidgets.QLabel("")
        self.title.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        self.header.addWidget(self.title, 1)

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self.header_bar)
        layout.addWidget(self.plot, 1)

        self._x = np.empty(0)
        self._y = np.empty(0)
        self._user_zoomed = False
        self._full_x = (0.0, 1.0)
        self.set_active(False)

        self.vb.sigXRangeChanged.connect(self._on_x_range_changed)
        # pyqtgraph's own double-click events are unreliable; catch them on the widget.
        self.plot.viewport().installEventFilter(self)
        self.vb.sigContextMenu.connect(self._show_context_menu)
        self.vb.sigUserZoomed.connect(self._mark_zoomed)
        self.plot.scene().sigMouseMoved.connect(self._on_mouse_moved)

    # --- view management --------------------------------------------------

    def set_active(self, active: bool) -> None:
        self.header_bar.setStyleSheet(ACTIVE_HEADER_STYLE if active else INACTIVE_HEADER_STYLE)

    def eventFilter(self, obj, event):
        if event.type() == QtCore.QEvent.Type.MouseButtonPress:
            self.activated.emit()
        if (event.type() == QtCore.QEvent.Type.MouseButtonDblClick
                and event.button() == QtCore.Qt.MouseButton.LeftButton):
            self.reset_view()
            return True
        return super().eventFilter(obj, event)

    def _mark_zoomed(self):
        self._user_zoomed = True

    def set_full_x(self, lo: float, hi: float) -> None:
        if hi <= lo:
            hi = lo + 1.0
        pad = (hi - lo) * 0.02
        self._full_x = (lo - pad, hi + pad)
        self.vb.setLimits(xMin=self._full_x[0], xMax=self._full_x[1])

    def reset_view(self) -> None:
        self._user_zoomed = False
        self.vb.setXRange(*self._full_x, padding=0)
        self._on_x_range_changed()  # setXRange is silent when the range is unchanged

    def _visible_slice(self) -> slice:
        lo, hi = self.vb.viewRange()[0]
        return slice(int(np.searchsorted(self._x, lo)), int(np.searchsorted(self._x, hi, side="right")))

    def auto_y(self) -> None:
        y = self._y[self._visible_slice()]
        top = float(y.max()) if len(y) else 0.0
        self.vb.setYRange(0, top * Y_HEADROOM if top > 0 else 1.0, padding=0)

    def _on_x_range_changed(self, *_):
        self.auto_y()

    # --- hover ------------------------------------------------------------

    def _on_mouse_moved(self, pos):
        if not len(self._x) or not self.vb.sceneBoundingRect().contains(pos):
            return
        x = self.vb.mapSceneToView(pos).x()
        self.hover.emit(self.hover_text(x))

    def hover_text(self, x: float) -> str:
        return ""

    def nearest(self, x: float) -> int:
        i = int(np.searchsorted(self._x, x))
        if i <= 0:
            return 0
        if i >= len(self._x):
            return len(self._x) - 1
        return i if self._x[i] - x < x - self._x[i - 1] else i - 1

    # --- context menu -----------------------------------------------------

    def _show_context_menu(self, screen_pos: QtCore.QPoint, x: float) -> None:
        menu = QtWidgets.QMenu(self)
        self.add_menu_actions(menu, x)
        menu.addAction("Copy image", self.copy_image)
        export = menu.addAction(f"Export {self.export_label()} as CSV…", self.export_csv)
        export.setEnabled(len(self._x) > 0)
        menu.exec(screen_pos)

    def add_menu_actions(self, menu: QtWidgets.QMenu, x: float) -> None:
        """Hook for pane-specific context menu items (x is the click in data units)."""

    def export_label(self) -> str:
        return "data"

    def copy_image(self) -> None:
        QtWidgets.QApplication.clipboard().setPixmap(self.grab())

    def export_csv(self) -> None:
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export CSV", str(Path.home() / self.export_name), "CSV files (*.csv)"
        )
        if not path:
            return
        np.savetxt(path, np.column_stack([self._x, self._y]), delimiter=",",
                   header=self.export_header, comments="", fmt="%.6f")


class ChromatogramPane(PlotPane):
    clicked = QtCore.Signal(float)  # RT
    rangeSelected = QtCore.Signal(float, float)

    export_header = "rt_min,intensity"

    def __init__(self, parent=None):
        super().__init__("Retention time (min)", select_drag=True, parent=parent)
        self.curve = pg.PlotDataItem(pen=CHROM_PEN)
        self.curve.setDownsampling(auto=True, method="peak")
        self.curve.setClipToView(True)
        self.plot.addItem(self.curve)

        self.cursor = pg.InfiniteLine(angle=90, movable=False, pen=CURSOR_PEN)
        self.region = pg.LinearRegionItem(movable=False, brush=RANGE_BRUSH, pen=pg.mkPen(None))
        self.ms2_marker = pg.InfiniteLine(angle=90, movable=False, pen=MS2_MARKER_PEN)
        for item in (self.cursor, self.region, self.ms2_marker):
            item.setZValue(10)
            item.hide()
            self.plot.addItem(item, ignoreBounds=True)

        self._scans = np.empty(0, dtype=np.int32)
        self.kind = "TIC"
        self.vb.setLimits(minXRange=0.01)
        self.vb.sigClicked.connect(lambda x, _y: self.clicked.emit(x))
        self.vb.sigRangeDrag.connect(self._on_range_drag)

    def set_data(self, rt, intensity, scans, kind: str, keep_view: bool = False,
                 title: str | None = None, export_name: str | None = None) -> None:
        self._x, self._y, self._scans, self.kind = rt, intensity, scans, kind
        self.export_name = export_name or f"{kind.lower().replace(' ', '_')}.csv"
        self.curve.setData(rt, intensity)
        self.title.setText(title or f"{kind} chromatogram (MS1)")
        self.plot.getPlotItem().setLabel("left", kind)
        if len(rt):
            self.set_full_x(float(rt[0]), float(rt[-1]))
        if keep_view:
            self.auto_y()
        else:
            self.reset_view()

    def clear(self) -> None:
        self.set_data(np.empty(0), np.empty(0), np.empty(0, dtype=np.int32), self.kind)
        self.cursor.hide()
        self.region.hide()
        self.ms2_marker.hide()

    def show_ms2_marker(self, rt: float | None) -> None:
        """Dashed line at the RT of the MS2 currently shown."""
        if rt is None:
            self.ms2_marker.hide()
        else:
            self.ms2_marker.setValue(rt)
            self.ms2_marker.show()

    def show_cursor(self, rt: float) -> None:
        self.region.hide()
        self.cursor.setValue(rt)
        self.cursor.show()

    def show_range(self, rt0: float, rt1: float) -> None:
        self.cursor.hide()
        self.region.setRegion((rt0, rt1))
        self.region.show()

    def _on_range_drag(self, x0: float, x1: float, finished: bool) -> None:
        self.show_range(x0, x1)
        if finished:
            self.rangeSelected.emit(x0, x1)

    def export_label(self) -> str:
        return "chromatogram"

    def hover_text(self, x: float) -> str:
        i = self.nearest(x)
        return f"RT {self._x[i]:.3f} min   {self.kind} {self._y[i]:.4g}   MS1 scan {self._scans[i]}"


class SpectrumPane(PlotPane):
    mzClicked = QtCore.Signal(float)
    xicRequested = QtCore.Signal(float)

    export_header = "mz,intensity"

    def __init__(self, pen: QtGui.QPen, parent=None):
        super().__init__("m/z", parent=parent)
        self.plot.getPlotItem().setLabel("left", "Intensity")
        self.sticks = pg.PlotCurveItem(pen=pen, connect="pairs")
        self.profile = pg.PlotCurveItem(pen=pen)
        self.plot.addItem(self.sticks)
        self.plot.addItem(self.profile)

        self.precursor_line = pg.InfiniteLine(angle=90, movable=False, pen=PRECURSOR_PEN)
        self.isolation = pg.LinearRegionItem(movable=False, brush=ISOLATION_BRUSH, pen=pg.mkPen(None))
        for item in (self.precursor_line, self.isolation):
            item.setZValue(-5)
            item.hide()
            self.plot.addItem(item, ignoreBounds=True)

        font = QtGui.QFont()
        font.setPointSize(9)
        self.labels = []
        for _ in range(LABEL_COUNT):
            label = pg.TextItem(color="#333", anchor=(0.5, 1.0))
            label.setFont(font)
            label.hide()
            self.plot.addItem(label, ignoreBounds=True)
            self.labels.append(label)

        self.centroid = True
        self.vb.setLimits(minXRange=0.01)
        self.vb.sigClicked.connect(lambda x, _y: self.mzClicked.emit(x))

    def set_spectrum(self, spec: Spectrum | None, title: str, export_name: str = "spectrum.csv") -> None:
        self.title.setText(title)
        self.export_name = export_name
        if spec is None:
            spec = Spectrum.empty()
        self._x, self._y, self.centroid = spec.mz, spec.intensity, spec.centroid
        if spec.centroid:
            self.profile.setData([], [])
            y = np.column_stack([np.zeros_like(spec.intensity), spec.intensity]).ravel()
            self.sticks.setData(np.repeat(spec.mz, 2), y)
        else:
            self.sticks.setData([], [])
            self.profile.setData(spec.mz, spec.intensity)
        if len(spec):
            self.set_full_x(float(spec.mz[0]), float(spec.mz[-1]))
        if self._user_zoomed and len(spec):
            self.auto_y()
            self._update_labels()
        else:
            self.reset_view()

    def set_precursor(self, mz: float | None, window: tuple[float, float] | None) -> None:
        if mz is None or not np.isfinite(mz):
            self.precursor_line.hide()
        else:
            self.precursor_line.setValue(mz)
            self.precursor_line.show()
        if window is None:
            self.isolation.hide()
        else:
            self.isolation.setRegion(window)
            self.isolation.show()

    def _on_x_range_changed(self, *_):
        self.auto_y()
        self._update_labels()

    def _update_labels(self) -> None:
        sl = self._visible_slice()
        x, y = self._x[sl], self._y[sl]
        if not self.centroid and len(y) > 2:
            apex = np.flatnonzero((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:])) + 1
        else:
            apex = np.flatnonzero(y > 0)
        order = apex[np.argsort(y[apex])[::-1][: LABEL_COUNT * 8]]

        lo, hi = self.vb.viewRange()[0]
        px_per_x = self.vb.width() / (hi - lo) if hi > lo else 0
        placed: list[float] = []
        for i in order:
            if len(placed) == LABEL_COUNT:
                break
            mz = x[i] if self.centroid else self._refine_apex(x, y, i)
            if px_per_x and any(abs(mz - p) * px_per_x < LABEL_MIN_PX for p in placed):
                continue
            label = self.labels[len(placed)]
            # Keep labels near the edges inside the plot.
            px = (mz - lo) * px_per_x
            edge = LABEL_MIN_PX / 2
            label.setAnchor((0.0 if px < edge else 1.0 if px > self.vb.width() - edge else 0.5, 1.0))
            label.setText(f"{mz:.4f}")
            label.setPos(mz, y[i])
            label.show()
            placed.append(mz)
        for label in self.labels[len(placed):]:
            label.hide()

    @staticmethod
    def _refine_apex(x, y, i) -> float:
        """Parabolic interpolation of a profile peak's apex m/z."""
        if i <= 0 or i >= len(y) - 1:
            return float(x[i])
        y0, y1, y2 = y[i - 1], y[i], y[i + 1]
        denom = y0 - 2 * y1 + y2
        if denom == 0:
            return float(x[i])
        offset = 0.5 * (y0 - y2) / denom  # in units of sample spacing, within ±0.5
        step = (x[i + 1] - x[i - 1]) / 2
        return float(x[i] + offset * step)

    def export_label(self) -> str:
        return "spectrum"

    def add_menu_actions(self, menu: QtWidgets.QMenu, x: float) -> None:
        mz = self.peak_near(x)
        if mz is not None:
            menu.addAction(f"XIC of m/z {mz:.4f}", lambda: self.xicRequested.emit(mz))
            menu.addSeparator()

    def peak_near(self, x: float, px: float = 8.0) -> float | None:
        """m/z of the tallest peak within `px` pixels of x, or None."""
        lo, hi = self.vb.viewRange()[0]
        if not len(self._x) or hi <= lo:
            return None
        half = px * (hi - lo) / max(self.vb.width(), 1)
        i0 = int(np.searchsorted(self._x, x - half))
        i1 = int(np.searchsorted(self._x, x + half, side="right"))
        if i1 <= i0 or self._y[i0:i1].max() <= 0:
            return None
        i = i0 + int(np.argmax(self._y[i0:i1]))
        return float(self._x[i]) if self.centroid else self._refine_apex(self._x, self._y, i)

    def hover_text(self, x: float) -> str:
        i = self.nearest(x)
        return f"m/z {self._x[i]:.4f}   intensity {self._y[i]:.4g}"


class MS2Pane(SpectrumPane):
    prevRequested = QtCore.Signal()
    nextRequested = QtCore.Signal()

    def __init__(self, parent=None):
        super().__init__(MS2_PEN, parent=parent)
        self.prev_button = QtWidgets.QToolButton(text="◀")
        self.next_button = QtWidgets.QToolButton(text="▶")
        for button, tip in ((self.prev_button, "Previous MS2 (↑)"), (self.next_button, "Next MS2 (↓)")):
            button.setToolTip(tip)
            button.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
            button.setAutoRepeat(True)
            self.header.addWidget(button)
        self.prev_button.clicked.connect(self.prevRequested)
        self.next_button.clicked.connect(self.nextRequested)

    def set_nav_enabled(self, prev: bool, next_: bool) -> None:
        self.prev_button.setEnabled(prev)
        self.next_button.setEnabled(next_)
