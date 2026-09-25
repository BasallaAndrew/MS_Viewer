"""Scan averaging in NumPy.

Centroid data: pool all peaks, cluster by a ppm tolerance, sum intensities per
cluster and divide by the number of scans; report the intensity-weighted mean
m/z per cluster.

Profile data: build a common m/z grid from the pooled sample positions
(near-duplicates merged), linearly interpolate each scan onto it and average.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

from msview.data.base import Spectrum


def average_centroid(spectra: Sequence[Spectrum], ppm: float = 5.0) -> Spectrum:
    n = len(spectra)
    if n == 0:
        return Spectrum.empty(centroid=True)
    mz = np.concatenate([s.mz for s in spectra])
    inten = np.concatenate([s.intensity for s in spectra]).astype(np.float64)
    keep = inten > 0
    mz, inten = mz[keep], inten[keep]
    if not len(mz):
        return Spectrum.empty(centroid=True)
    order = np.argsort(mz, kind="stable")
    mz, inten = mz[order], inten[order]
    # Start a new cluster wherever the gap to the previous peak exceeds tolerance.
    gap_ppm = np.diff(mz) / mz[:-1] * 1e6
    starts = np.concatenate(([0], np.flatnonzero(gap_ppm > ppm) + 1))
    summed = np.add.reduceat(inten, starts)
    weighted_mz = np.add.reduceat(mz * inten, starts) / summed
    return Spectrum(weighted_mz, summed / n, centroid=True)


def _median_rel_spacing(mz: np.ndarray) -> float:
    if len(mz) < 3:
        return np.nan
    rel = np.diff(mz) / mz[:-1]
    rel = rel[rel > 0]
    return float(np.median(rel)) if len(rel) else np.nan


def average_profile(spectra: Sequence[Spectrum], gap_factor: float = 8.0) -> Spectrum:
    """Average profile spectra on a shared grid.

    Thermo profile scans only store points around peaks (bracketed by zeros),
    and each scan's sampling positions differ slightly. The grid is therefore
    the union of all sample positions with near-duplicates (closer than a
    quarter of the typical spacing) merged. Gaps much wider than the typical
    spacing are treated as zero rather than interpolated across, in case the
    writer dropped the bracketing zeros.
    """
    n = len(spectra)
    nonempty = [s for s in spectra if len(s.mz) > 1]
    if not nonempty:
        return Spectrum.empty(centroid=False)
    spacings = [_median_rel_spacing(s.mz) for s in nonempty]
    spacings = [x for x in spacings if np.isfinite(x) and x > 0]
    spacing = float(np.median(spacings)) if spacings else 1e-6

    pooled = np.sort(np.concatenate([s.mz for s in nonempty]))
    keep = np.concatenate(([True], np.diff(pooled) / pooled[:-1] > spacing / 4))
    grid = pooled[keep]

    acc = np.zeros(len(grid))
    for s in nonempty:
        vals = np.interp(grid, s.mz, s.intensity, left=0.0, right=0.0)
        # Wide gaps only matter when an edge is nonzero; zero-bracketed gaps
        # already interpolate to zero.
        rel = np.diff(s.mz) / s.mz[:-1]
        edge_nonzero = (s.intensity[:-1] > 0) | (s.intensity[1:] > 0)
        gaps = np.flatnonzero((rel > gap_factor * spacing) & edge_nonzero)
        if len(gaps):
            lo = np.searchsorted(grid, s.mz[gaps], side="right")
            hi = np.searchsorted(grid, s.mz[gaps + 1], side="left")
            mask = np.zeros(len(grid) + 1, dtype=np.int32)
            np.add.at(mask, lo, 1)
            np.add.at(mask, hi, -1)
            vals[np.cumsum(mask[:-1]) > 0] = 0.0
        acc += vals
    acc /= n

    # Drop long zero runs, but keep the zeros that bracket each peak.
    nz = acc > 0
    keep = nz.copy()
    keep[1:] |= nz[:-1]
    keep[:-1] |= nz[1:]
    return Spectrum(grid[keep], acc[keep], centroid=False)


def average_spectra(spectra: Sequence[Spectrum], ppm: float = 5.0) -> Spectrum:
    """Average scans, choosing the method from the data type.

    Mixed centroid/profile input (unusual) falls back to centroid binning.
    """
    if spectra and all(not s.centroid for s in spectra):
        return average_profile(spectra)
    return average_centroid(spectra, ppm=ppm)
