"""Scan index: one row per scan, built once when a file is opened.

All selection logic (which MS1 is nearest a click, which MS2s belong to a
duty cycle or an RT range) lives here as pure functions over the index so it
can be tested without a GUI or a data file.
"""

from __future__ import annotations

import numpy as np

INDEX_DTYPE = np.dtype(
    [
        ("scan", np.int32),
        ("rt", np.float64),  # minutes
        ("ms_level", np.int8),
        ("precursor_mz", np.float64),  # monoisotopic when known; NaN for MS1
        ("charge", np.int16),  # 0 if unknown
        ("isolation_mz", np.float64),  # centre of the isolation window; NaN for MS1
        ("isolation_width", np.float32),  # full width in m/z; NaN if unknown
        ("activation", "U8"),  # HCD / CID / ETD ...
        ("collision_energy", np.float32),
        ("parent_ms1", np.int32),  # scan number of the preceding MS1, -1 if none
        ("tic", np.float64),
        ("base_peak", np.float64),
        ("centroid", np.bool_),
    ]
)


def empty_index(n: int) -> np.ndarray:
    idx = np.zeros(n, dtype=INDEX_DTYPE)
    for field in ("precursor_mz", "isolation_mz", "isolation_width", "collision_energy"):
        idx[field] = np.nan
    idx["parent_ms1"] = -1
    return idx


def finalize(idx: np.ndarray) -> np.ndarray:
    """Sort by scan number and fill in parent_ms1."""
    idx = np.sort(idx, order="scan")
    is_ms1 = idx["ms_level"] == 1
    # Position of the most recent MS1 at or before each row.
    last = np.where(is_ms1, np.arange(len(idx)), -1)
    last = np.maximum.accumulate(last) if len(idx) else last
    idx["parent_ms1"] = np.where(last >= 0, idx["scan"][np.maximum(last, 0)], -1)
    # Fall back to the isolation centre when the precursor m/z is missing, and
    # vice versa, so markers always have something to draw.
    msn = ~is_ms1
    p, c = idx["precursor_mz"], idx["isolation_mz"]
    idx["precursor_mz"] = np.where(msn & np.isnan(p), c, p)
    idx["isolation_mz"] = np.where(msn & np.isnan(c), idx["precursor_mz"], c)
    return idx


def row_of(idx: np.ndarray, scan: int) -> int:
    i = int(np.searchsorted(idx["scan"], scan))
    if i >= len(idx) or idx["scan"][i] != scan:
        raise KeyError(f"scan {scan} not in index")
    return i


def ms1_rows(idx: np.ndarray) -> np.ndarray:
    return np.flatnonzero(idx["ms_level"] == 1)


def nearest_ms1(idx: np.ndarray, rt: float) -> int | None:
    rows = ms1_rows(idx)
    if not len(rows):
        return None
    rts = idx["rt"][rows]
    j = int(np.searchsorted(rts, rt))
    candidates = [k for k in (j - 1, j) if 0 <= k < len(rows)]
    best = min(candidates, key=lambda k: abs(rts[k] - rt))
    return int(idx["scan"][rows[best]])


def step_ms1(idx: np.ndarray, scan: int, step: int) -> int | None:
    """MS1 scan `step` MS1s away from `scan` (clamped to the run)."""
    ms1_scans = idx["scan"][ms1_rows(idx)]
    if not len(ms1_scans):
        return None
    j = int(np.searchsorted(ms1_scans, scan))
    if j >= len(ms1_scans) or ms1_scans[j] != scan:
        # `scan` isn't an MS1; step from where it would sit.
        j = j - 1 if step > 0 else j
    j = int(np.clip(j + step, 0, len(ms1_scans) - 1))
    return int(ms1_scans[j])


def ms2_for_cycle(idx: np.ndarray, ms1_scan: int) -> np.ndarray:
    """Scan numbers of MSn scans triggered by `ms1_scan` (up to the next MS1)."""
    mask = (idx["ms_level"] >= 2) & (idx["parent_ms1"] == ms1_scan)
    return idx["scan"][mask]


def ms1_in_range(idx: np.ndarray, rt0: float, rt1: float) -> np.ndarray:
    mask = (idx["ms_level"] == 1) & (idx["rt"] >= rt0) & (idx["rt"] <= rt1)
    return idx["scan"][mask]


def ms2_in_range(idx: np.ndarray, rt0: float, rt1: float) -> np.ndarray:
    mask = (idx["ms_level"] >= 2) & (idx["rt"] >= rt0) & (idx["rt"] <= rt1)
    return idx["scan"][mask]


def evenly_subsample(scans: np.ndarray, cap: int) -> np.ndarray:
    if len(scans) <= cap:
        return scans
    picks = np.linspace(0, len(scans) - 1, cap).round().astype(int)
    return scans[np.unique(picks)]


def isolation_window(row) -> tuple[float, float] | None:
    """(low, high) m/z of an MSn row's isolation window, if known."""
    centre, width = float(row["isolation_mz"]), float(row["isolation_width"])
    if np.isnan(centre) or np.isnan(width) or width <= 0:
        return None
    return centre - width / 2, centre + width / 2


ISOTOPE_SPACING = 1.003355  # 13C - 12C


def precursor_matches(idx: np.ndarray, rows: np.ndarray, mz: float, ppm: float = 20.0) -> tuple[np.ndarray, np.ndarray]:
    """For MSn rows, whether `mz` is (a) the precursor or one of its first few
    13C isotopes within ppm, and (b) inside the isolation window."""
    r = idx[rows]
    tol = mz * ppm * 1e-6
    prec, z = r["precursor_mz"], r["charge"].astype(np.float64)
    k = np.where(z > 0, np.round((mz - prec) * z / ISOTOPE_SPACING), 0)
    k = np.clip(k, 0, 4)
    expected = prec + np.where(z > 0, k * ISOTOPE_SPACING / np.where(z > 0, z, 1), 0)
    exact = np.abs(mz - expected) <= tol
    half = np.where(np.isfinite(r["isolation_width"]), r["isolation_width"] / 2, 0.5)
    in_window = np.abs(r["isolation_mz"] - mz) <= half
    return exact, in_window


def ms2_for_precursor(idx: np.ndarray, mz: float, rt: float, rt_window: float,
                      ppm: float = 20.0) -> tuple[np.ndarray, int]:
    """MSn scans within rt ± rt_window that fragmented `mz`.

    Scans whose precursor (or an isotope of it) matches within ppm are
    preferred; if there are none, any scan whose isolation window contained
    mz is used. Returns (scans sorted by scan number, position of the scan
    closest in RT), or (empty, -1).
    """
    rows = np.flatnonzero((idx["ms_level"] >= 2) & (np.abs(idx["rt"] - rt) <= rt_window))
    if not len(rows):
        return idx["scan"][:0], -1
    exact, in_window = precursor_matches(idx, rows, mz, ppm)
    keep = exact if exact.any() else in_window
    rows = rows[keep]
    if not len(rows):
        return idx["scan"][:0], -1
    return idx["scan"][rows], int(np.argmin(np.abs(idx["rt"][rows] - rt)))
