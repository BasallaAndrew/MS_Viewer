"""Method and run information: data classes plus vendor-text parsers.

Backends fill a MethodInfo; everything here is plain Python/NumPy so the
parsers can be tested against real method text without a .raw file.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

GRADIENT_DTYPE = np.dtype([("time", np.float64), ("percent_b", np.float64), ("flow", np.float64)])


@dataclass
class Trace:
    """A logged quantity over run time, e.g. pump pressure."""

    name: str  # e.g. "Pressure (bar)"
    device: str
    rt: np.ndarray  # minutes
    values: np.ndarray

    @property
    def label(self) -> str:
        return f"{self.name}  [{self.device}]" if self.device else self.name


@dataclass
class MethodInfo:
    # Flat rows of (section, group, parameter, value); the GUI builds the tree.
    rows: list[tuple[str, str, str, str]] = field(default_factory=list)
    method_text: list[tuple[str, str]] = field(default_factory=list)  # (device, text)
    gradient: np.ndarray | None = None  # programmed gradient, GRADIENT_DTYPE
    traces: list[Trace] = field(default_factory=list)

    def add(self, section: str, group: str, items) -> None:
        for key, value in items:
            value = "" if value is None else str(value).strip()
            if value:
                self.rows.append((section, group, str(key).strip(), value))

    def percent_b_trace(self) -> Trace | None:
        """The logged %B trace, if a pump recorded one."""
        for trace in self.traces:
            if is_percent_b(trace.name):
                return trace
        return None


def is_percent_b(name: str) -> bool:
    return bool(re.search(r"composition\s*b\b|%\s*b\b|\bb\s*\(%\)", name, re.I))


# --- method text -------------------------------------------------------------

_WINDOWS_PATH = re.compile(r"^(.*?\S)\s+([A-Za-z]:\\.*)$")
_GAP = re.compile(r"\s{2,}")
_KEY_GAP_VALUE = re.compile(r"^(\S.*?)\s{2,}(\S.*)$")
_KEY_COLON_VALUE = re.compile(r"^([^:]{1,80}?)\s*:\s*(.*)$")
# Thermo pads "Label   value", but long labels can leave a single space before a flag.
_KEY_FLAG = re.compile(r"^(.*\S)\s+(True|False|On|Off|Yes|No)$")


def _is_table_row(text: str) -> bool:
    return text[0].isdigit() or text[0] in "[(-" or len(_GAP.findall(text)) >= 2


def _is_major_header(text: str) -> bool:
    letters = [c for c in text if c.isalpha()]
    return bool(letters) and all(c.isupper() for c in letters)


def parse_method_text(text: str) -> list[tuple[str, str, str]]:
    """Best-effort (group, parameter, value) rows from instrument-method text.

    Handles "Key: value" (Agilent, Dionex) and "Label      value" (Thermo
    MS methods) layouts. Table rows are skipped (they stay readable in the
    Method text tab); headings become groups.
    """
    rows: list[tuple[str, str, str]] = []
    major = group = ""
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        m = _WINDOWS_PATH.match(s)
        if m:
            rows.append((group, m.group(1), m.group(2)))
            continue
        if _is_table_row(s):
            continue
        m = _KEY_GAP_VALUE.match(s)
        if m:
            value = m.group(2).strip()
            if value not in ("―", "-"):
                rows.append((group, m.group(1), value))
            continue
        m = _KEY_COLON_VALUE.match(s)
        if m:
            # "Key:" alone is an empty setting or a table title; either way no row.
            if m.group(2).strip():
                rows.append((group, m.group(1), m.group(2).strip()))
            continue
        m = _KEY_FLAG.match(s)
        if m:
            rows.append((group, m.group(1), m.group(2)))
        elif s.endswith("."):
            rows.append((group, "Note", s))
        elif _is_major_header(s):
            major = group = s
        else:
            group = f"{major} › {s}" if major else s
    return rows


# --- programmed gradient -----------------------------------------------------

# Agilent:  "  30.00(min)    0.300(ml/min)      A=86.0% B=14.0% C=0.0% D=0.0%"
_AGILENT_ROW = re.compile(
    r"^\s*(\d+(?:\.\d+)?)\s*\(min\)\s+(\d+(?:\.\d+)?)\s*\((ml|ul|µl|nl)/min\)\s+(.*)$", re.I
)
_AGILENT_B = re.compile(r"\bB\s*=\s*(\d+(?:\.\d+)?)\s*%")
# Dionex / Vanquish instrument scripts: "  5.000 [min]  Pump.%B.Value:  10.0 [%]"
_SCRIPT_B = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*(?:\[min\])?\s+\S*%B(?:\.Value)?\s*[:=]\s*(\d+(?:\.\d+)?)", re.I)
_FLOW_SCALE = {"ml": 1.0, "ul": 1e-3, "µl": 1e-3, "nl": 1e-6}


def parse_gradient(text: str) -> np.ndarray | None:
    """Programmed %B vs time from a pump method, or None if not recognised.
    Flow is in mL/min (NaN when the format doesn't give it)."""
    points = []
    for line in text.splitlines():
        m = _AGILENT_ROW.match(line)
        if m:
            b = _AGILENT_B.search(m.group(4))
            if b:
                flow = float(m.group(2)) * _FLOW_SCALE[m.group(3).lower()]
                points.append((float(m.group(1)), float(b.group(1)), flow))
            continue
        m = _SCRIPT_B.match(line)
        if m:
            points.append((float(m.group(1)), float(m.group(2)), np.nan))
    if len(points) < 2:
        return None
    gradient = np.array(points, dtype=GRADIENT_DTYPE)
    return gradient[np.argsort(gradient["time"], kind="stable")]


# --- summary from the scan index ---------------------------------------------

def _fmt_counts(values, labels=None, limit=8) -> str:
    uniq, counts = np.unique(values, return_counts=True)
    order = np.argsort(counts)[::-1][:limit]
    parts = [f"{labels(uniq[i]) if labels else uniq[i]} ({counts[i]:,})" for i in order]
    return ", ".join(parts)


def summary_rows(idx: np.ndarray) -> list[tuple[str, str, str]]:
    """(group, parameter, value) rows describing the acquisition, from the index."""
    rows: list[tuple[str, str, str]] = []
    if not len(idx):
        return rows
    ms1 = idx[idx["ms_level"] == 1]
    msn = idx[idx["ms_level"] >= 2]
    rows += [
        ("Scans", "Total scans", f"{len(idx):,}"),
        ("Scans", "MS1 scans", f"{len(ms1):,}"),
        ("Scans", "MSn scans", f"{len(msn):,}"),
        ("Scans", "Retention time", f"{idx['rt'][0]:.2f} – {idx['rt'][-1]:.2f} min"),
    ]
    for name, sub in (("MS1", ms1), ("MSn", msn)):
        if len(sub):
            n_centroid = int(sub["centroid"].sum())
            kind = "centroid" if n_centroid == len(sub) else "profile" if n_centroid == 0 else "mixed"
            rows.append(("Scans", f"{name} data type (as stored)", kind))
    if len(ms1) > 1:
        cycle = np.median(np.diff(ms1["rt"])) * 60
        rows.append(("Duty cycle", "Median MS1 cycle time", f"{cycle:.2f} s"))
    if len(ms1) and len(msn):
        # Bin k counts MSn scans after the k-th MS1 (bin 0: before the first MS1).
        parent = np.searchsorted(ms1["scan"], msn["scan"])
        per_cycle = np.bincount(parent, minlength=len(ms1) + 1)[1:]
        rows.append(("Duty cycle", "MSn per MS1 cycle (mean / max)", f"{per_cycle.mean():.1f} / {per_cycle.max()}"))
    if len(msn):
        energies = [f"{a} {e:g}" if np.isfinite(e) else str(a)
                    for a, e in zip(msn["activation"], msn["collision_energy"])]
        rows.append(("Fragmentation", "Activation / energy", _fmt_counts(np.array(energies))))
        widths = msn["isolation_width"][np.isfinite(msn["isolation_width"])]
        if len(widths):
            rows.append(("Fragmentation", "Isolation width (m/z)", _fmt_counts(np.round(widths, 3), lambda w: f"{w:g}")))
        rows.append(("Fragmentation", "Precursor charge", _fmt_counts(msn["charge"], lambda z: str(z) if z else "unknown")))
        mz = msn["precursor_mz"][np.isfinite(msn["precursor_mz"])]
        if len(mz):
            rows.append(("Fragmentation", "Precursor m/z range", f"{mz.min():.2f} – {mz.max():.2f}"))
    return rows
