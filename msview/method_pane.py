"""Method view: gradient plot on top, run/method parameters below.

Shown in place of the MS1/MS2 panes when View → Method is selected.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets

from msview.method import MethodInfo, Trace
from msview.panes import CURSOR_PEN, RANGE_BRUSH, PlotPane

PROGRAM_PEN = pg.mkPen("#1f4e9c", width=2, style=QtCore.Qt.PenStyle.DashLine)
LOGGED_PEN = pg.mkPen("#2ca02c", width=2.5)
RIGHT_PEN = pg.mkPen("#e67e00", width=1.2)
EXPAND_GROUP_MAX_ROWS = 30  # bigger groups (e.g. calibration constants) start collapsed


class GradientPane(PlotPane):
    """Programmed and logged %B on the left axis; any logged trace on the right."""

    export_header = "time_min,programmed_percent_b,logged_percent_b,right_axis"
    export_name = "gradient.csv"

    def __init__(self, parent=None):
        super().__init__("Time (min)", parent=parent)
        pi = self.plot.getPlotItem()
        pi.setLabel("left", "%B")
        self.program = pg.PlotDataItem(pen=PROGRAM_PEN, symbol="o", symbolSize=6,
                                       symbolBrush=PROGRAM_PEN.color(), symbolPen=None)
        self.logged = pg.PlotDataItem(pen=LOGGED_PEN)
        self.plot.addItem(self.logged)
        self.plot.addItem(self.program)

        # Second y-axis: its own ViewBox sharing the x-axis, behind the main one
        # so the main ViewBox still gets the mouse.
        self.right_vb = pg.ViewBox(enableMenu=False)
        self.right_vb.setMouseEnabled(False, False)
        self.right_vb.setZValue(-10)
        pi.showAxis("right")
        pi.scene().addItem(self.right_vb)
        pi.getAxis("right").linkToView(self.right_vb)
        pi.getAxis("right").enableAutoSIPrefix(False)
        self.right_vb.setXLink(self.vb)
        self.right_curve = pg.PlotCurveItem(pen=RIGHT_PEN)
        self.right_vb.addItem(self.right_curve)
        self.vb.sigResized.connect(self._sync_right_geometry)

        self.cursor = pg.InfiniteLine(angle=90, movable=False, pen=CURSOR_PEN)
        self.region = pg.LinearRegionItem(movable=False, brush=RANGE_BRUSH, pen=pg.mkPen(None))
        for item in (self.cursor, self.region):
            item.setZValue(10)
            item.hide()
            self.plot.addItem(item, ignoreBounds=True)

        self.trace_box = QtWidgets.QComboBox()
        self.trace_box.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        self.trace_box.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToContents)
        self.trace_box.currentIndexChanged.connect(self._on_trace_changed)
        self.header.addWidget(QtWidgets.QLabel("Right axis:"))
        self.header.addWidget(self.trace_box)

        self.vb.setLimits(minXRange=0.01)
        self._program_xy = (np.empty(0), np.empty(0))
        self._logged: Trace | None = None
        self._right: Trace | None = None
        self._traces: list[Trace] = []
        self._title = "Gradient"
        self.set_info(None, "Gradient")

    def _sync_right_geometry(self):
        self.right_vb.setGeometry(self.vb.sceneBoundingRect())
        self.right_vb.linkedViewChanged(self.vb, self.right_vb.XAxis)

    # --- data ---------------------------------------------------------------

    def set_info(self, info: MethodInfo | None, message: str = "") -> None:
        gradient = info.gradient if info is not None else None
        if gradient is not None:
            self._program_xy = (gradient["time"], gradient["percent_b"])
        else:
            self._program_xy = (np.empty(0), np.empty(0))
        self.program.setData(*self._program_xy)
        self._logged = info.percent_b_trace() if info is not None else None
        if self._logged is not None:
            self.logged.setData(self._logged.rt, self._logged.values)
        else:
            self.logged.setData([], [])

        # Right-axis choices: everything except the %B trace already on the left.
        self._traces = [t for t in (info.traces if info else []) if t is not self._logged]
        self.trace_box.blockSignals(True)
        self.trace_box.clear()
        for trace in self._traces:
            self.trace_box.addItem(trace.label)
        default = next((i for i, t in enumerate(self._traces) if "pressure" in t.name.lower()), 0)
        self.trace_box.setCurrentIndex(default if self._traces else -1)
        self.trace_box.setEnabled(bool(self._traces))
        self.trace_box.blockSignals(False)
        self._on_trace_changed(self.trace_box.currentIndex())

        times = [x for x in (self._program_xy[0], self._logged.rt if self._logged else None) if x is not None and len(x)]
        times += [t.rt for t in self._traces]
        if times:
            self.set_full_x(min(float(t.min()) for t in times), max(float(t.max()) for t in times))
        self._x = np.unique(np.concatenate(times)) if times else np.empty(0)
        self._y = np.zeros_like(self._x)  # unused; hover/export are overridden

        has_b = gradient is not None or self._logged is not None
        self.plot.getPlotItem().showAxis("left", has_b)
        if message:
            title = message
        elif info is None:
            title = "Gradient"
        elif gradient is None and self._logged is None:
            title = "Gradient: not recorded in this file" + ("; showing logged traces" if self._traces else "")
        else:
            parts = ["programmed %B" if gradient is not None else "",
                     "logged %B" if self._logged is not None else ""]
            title = "Gradient: " + " and ".join(p for p in parts if p)
        self._title = title
        self._rebuild_legend()
        self.reset_view()

    def _on_trace_changed(self, index: int) -> None:
        self._right = self._traces[index] if 0 <= index < len(self._traces) else None
        if self._right is None:
            self.right_curve.setData([], [])
            self.plot.getPlotItem().setLabel("right", "")
        else:
            self.right_curve.setData(self._right.rt, self._right.values)
            self.plot.getPlotItem().setLabel("right", self._right.name)
        self._rebuild_legend()
        self.auto_y()

    def _rebuild_legend(self) -> None:
        """Colour key in the header (an in-plot legend would cover the traces)."""
        key = []
        if len(self._program_xy[0]):
            key.append((PROGRAM_PEN, "programmed %B (dashed)"))
        if self._logged is not None:
            key.append((LOGGED_PEN, "logged %B"))
        if self._right is not None:
            key.append((RIGHT_PEN, f"{self._right.name} (right axis)"))
        swatches = "&nbsp;&nbsp;&nbsp;".join(
            f'<span style="color:{pen.color().name()}">&#9632;</span>&nbsp;{label}' for pen, label in key)
        self.title.setText(f"{self._title}&nbsp;&nbsp;&nbsp;&nbsp;{swatches}" if swatches else self._title)

    # --- selection markers (follow the chromatogram) --------------------------

    def show_cursor(self, rt: float) -> None:
        self.region.hide()
        self.cursor.setValue(rt)
        self.cursor.show()

    def show_range(self, rt0: float, rt1: float) -> None:
        self.cursor.hide()
        self.region.setRegion((rt0, rt1))
        self.region.show()

    # --- view -----------------------------------------------------------------

    def auto_y(self) -> None:
        self.vb.setYRange(0, 105, padding=0)
        if self._right is None:
            return
        lo, hi = self.vb.viewRange()[0]
        mask = (self._right.rt >= lo) & (self._right.rt <= hi)
        values = self._right.values[mask] if mask.any() else self._right.values
        values = values[np.isfinite(values)]
        if not len(values):
            return
        vmin, vmax = float(values.min()), float(values.max())
        pad = (vmax - vmin) * 0.05 or abs(vmax) * 0.05 or 1.0
        self.right_vb.setYRange(vmin - pad, vmax + pad, padding=0)

    def _at(self, x: float) -> dict[str, float]:
        values: dict[str, float] = {}
        t, b = self._program_xy
        if len(t) and t[0] <= x <= t[-1]:
            values["programmed %B"] = float(np.interp(x, t, b))
        for name, trace in (("logged %B", self._logged), (None, self._right)):
            if trace is not None and trace.rt[0] <= x <= trace.rt[-1]:
                values[name or trace.name] = float(np.interp(x, trace.rt, trace.values))
        return values

    def hover_text(self, x: float) -> str:
        parts = [f"t {x:.2f} min"] + [f"{k} {v:.4g}" for k, v in self._at(x).items()]
        return "   ".join(parts)

    def export_label(self) -> str:
        return "gradient"

    def export_csv(self) -> None:
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export CSV", str(Path.home() / self.export_name), "CSV files (*.csv)"
        )
        if not path:
            return
        base = self._logged or self._right
        t = base.rt if base is not None else self._program_xy[0]

        def column(xs, ys):
            if not len(xs):
                return np.full(len(t), np.nan)
            return np.where((t >= xs[0]) & (t <= xs[-1]), np.interp(t, xs, ys), np.nan)

        columns = [t, column(*self._program_xy)]
        columns.append(column(self._logged.rt, self._logged.values) if self._logged else np.full(len(t), np.nan))
        columns.append(column(self._right.rt, self._right.values) if self._right else np.full(len(t), np.nan))
        right_name = self._right.name.replace(",", ";") if self._right else "right_axis"
        header = self.export_header.replace("right_axis", right_name)
        np.savetxt(path, np.column_stack(columns), delimiter=",", header=header, comments="", fmt="%.6g")


class ParameterTree(QtWidgets.QWidget):
    """Searchable two-column tree: section › group › parameter = value."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.filter = QtWidgets.QLineEdit()
        self.filter.setPlaceholderText("Filter parameters (e.g. resolution, AGC, temperature)…")
        self.filter.setClearButtonEnabled(True)
        self.filter.textChanged.connect(self._apply_filter)
        self.tree = QtWidgets.QTreeWidget()
        self.tree.setHeaderLabels(["Parameter", "Value"])
        self.tree.setAlternatingRowColors(True)
        self.tree.setUniformRowHeights(True)
        self.tree.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.tree.setContextMenuPolicy(QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self._context_menu)
        copy = QtGui.QShortcut(QtGui.QKeySequence.StandardKey.Copy, self.tree)
        copy.setContext(QtCore.Qt.ShortcutContext.WidgetShortcut)
        copy.activated.connect(self.copy_selected)

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(self.filter)
        layout.addWidget(self.tree)

    def set_rows(self, rows: list[tuple[str, str, str, str]], placeholder: str = "") -> None:
        self.tree.clear()
        if not rows:
            if placeholder:
                QtWidgets.QTreeWidgetItem(self.tree, [placeholder, ""])
            return
        bold = QtGui.QFont()
        bold.setBold(True)
        sections: dict[str, QtWidgets.QTreeWidgetItem] = {}
        groups: dict[tuple[str, str], QtWidgets.QTreeWidgetItem] = {}
        for section, group, key, value in rows:
            parent = sections.get(section)
            if parent is None:
                parent = sections[section] = QtWidgets.QTreeWidgetItem(self.tree, [section, ""])
                parent.setFont(0, bold)
                parent.setFirstColumnSpanned(True)
            if group:
                group_item = groups.get((section, group))
                if group_item is None:
                    group_item = groups[(section, group)] = QtWidgets.QTreeWidgetItem(parent, [group, ""])
                    group_item.setFirstColumnSpanned(True)
                parent = group_item
            item = QtWidgets.QTreeWidgetItem(parent, [key, value])
            item.setToolTip(1, value)
        self._default_expansion()
        self.tree.resizeColumnToContents(0)
        self.tree.setColumnWidth(0, min(self.tree.columnWidth(0) + 20, 420))
        self._apply_filter(self.filter.text())

    def _default_expansion(self) -> None:
        for i in range(self.tree.topLevelItemCount()):
            section = self.tree.topLevelItem(i)
            section.setExpanded(True)
            for j in range(section.childCount()):
                child = section.child(j)
                child.setExpanded(child.childCount() <= EXPAND_GROUP_MAX_ROWS)

    def _apply_filter(self, text: str) -> None:
        needle = text.strip().lower()

        def visit(item: QtWidgets.QTreeWidgetItem, context: str) -> bool:
            label = f"{context} {item.text(0)} {item.text(1)}".lower()
            if item.childCount() == 0:
                visible = not needle or needle in label
            else:
                # A matching section/group shows all of its rows.
                if needle and needle in item.text(0).lower():
                    visible = True
                    for k in range(item.childCount()):
                        self._show_all(item.child(k))
                else:
                    visible = False
                    for k in range(item.childCount()):
                        visible |= visit(item.child(k), label)
                if needle and visible:
                    item.setExpanded(True)
            item.setHidden(not visible)
            return visible

        for i in range(self.tree.topLevelItemCount()):
            visit(self.tree.topLevelItem(i), "")
        if not needle:
            self._default_expansion()

    def _show_all(self, item: QtWidgets.QTreeWidgetItem) -> None:
        item.setHidden(False)
        for k in range(item.childCount()):
            self._show_all(item.child(k))

    def _context_menu(self, pos: QtCore.QPoint) -> None:
        menu = QtWidgets.QMenu(self)
        menu.addAction("Copy selected", self.copy_selected)
        menu.addAction("Copy all (visible)", self.copy_all)
        menu.exec(self.tree.viewport().mapToGlobal(pos))

    @staticmethod
    def _line(item: QtWidgets.QTreeWidgetItem) -> str:
        return f"{item.text(0)}\t{item.text(1)}" if item.text(1) else item.text(0)

    def copy_selected(self) -> None:
        lines = [self._line(item) for item in self.tree.selectedItems()]
        QtWidgets.QApplication.clipboard().setText("\n".join(lines))

    def copy_all(self) -> None:
        lines: list[str] = []

        def walk(item: QtWidgets.QTreeWidgetItem, depth: int) -> None:
            if item.isHidden():
                return
            lines.append("    " * depth + self._line(item))
            for k in range(item.childCount()):
                walk(item.child(k), depth + 1)

        for i in range(self.tree.topLevelItemCount()):
            walk(self.tree.topLevelItem(i), 0)
        QtWidgets.QApplication.clipboard().setText("\n".join(lines))


class MethodPane(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.gradient = GradientPane()
        self.parameters = ParameterTree()
        self.text = QtWidgets.QPlainTextEdit()
        self.text.setReadOnly(True)
        self.text.setLineWrapMode(QtWidgets.QPlainTextEdit.LineWrapMode.NoWrap)
        self.text.setFont(QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.SystemFont.FixedFont))

        self.tabs = QtWidgets.QTabWidget()
        self.tabs.addTab(self.parameters, "Parameters")
        self.tabs.addTab(self.text, "Method text")

        splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        splitter.setChildrenCollapsible(False)
        splitter.addWidget(self.gradient)
        splitter.addWidget(self.tabs)
        splitter.setSizes([350, 450])
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(splitter)
        self.set_info(None, "Open a file to see its method and run information")

    def set_info(self, info: MethodInfo | None, message: str = "") -> None:
        self.gradient.set_info(info, message if info is None else "")
        self.parameters.set_rows(info.rows if info else [], placeholder=message)
        if info is None:
            self.text.setPlainText(message)
        elif info.method_text:
            blocks = [f"{'=' * 78}\n{name}\n{'=' * 78}\n{text.rstrip()}\n" for name, text in info.method_text]
            self.text.setPlainText("\n".join(blocks))
        else:
            self.text.setPlainText("This file doesn't include the instrument method text.\n"
                                   "(mzML conversions don't carry it; open the original .raw file to see it.)")
