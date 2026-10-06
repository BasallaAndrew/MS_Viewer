"""The RawSource interface. The GUI talks only to this, never to a backend."""

from __future__ import annotations

import threading
from abc import ABC, abstractmethod
from collections import OrderedDict
from dataclasses import dataclass
from typing import Sequence

import numpy as np


class BackendUnavailable(RuntimeError):
    """Raised when a backend's dependencies (DLLs, runtimes) are missing."""


@dataclass(frozen=True)
class Spectrum:
    mz: np.ndarray  # float64, ascending
    intensity: np.ndarray  # float64
    centroid: bool

    @classmethod
    def empty(cls, centroid: bool = True) -> "Spectrum":
        return cls(np.empty(0), np.empty(0), centroid)

    def __len__(self) -> int:
        return len(self.mz)


class RawSource(ABC):
    """A readable LC-MS run.

    Subclasses implement `_build_index` and `_read_spectrum`; this class adds
    the LRU spectrum cache, chromatogram extraction and averaging. Spectrum
    reads must be safe to call from a worker thread.
    """

    cache_size = 500

    def __init__(self, path: str):
        self.path = path
        self._index: np.ndarray | None = None
        self._cache: OrderedDict[int, Spectrum] = OrderedDict()
        self._cache_lock = threading.Lock()

    # --- to implement -----------------------------------------------------

    @abstractmethod
    def _build_index(self) -> np.ndarray:
        """Return a finalized array of msview.index.INDEX_DTYPE."""

    @abstractmethod
    def _read_spectrum(self, scan: int) -> Spectrum: ...

    def close(self) -> None:
        pass

    # --- public API -------------------------------------------------------

    def scan_index(self) -> np.ndarray:
        if self._index is None:
            self._index = self._build_index()
        return self._index

    def chromatogram(self, kind: str = "tic") -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(rt, intensity, scan) over MS1 scans; kind is 'tic' or 'bpc'."""
        idx = self.scan_index()
        ms1 = idx[idx["ms_level"] == 1]
        field = "tic" if kind == "tic" else "base_peak"
        return ms1["rt"], ms1[field], ms1["scan"]

    def xic(self, mz: float, tolerance: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Extracted ion chromatogram over MS1 scans: summed intensity within
        mz ± tolerance (m/z units). Returns (rt, intensity, scan).

        Reads spectra without the LRU cache so a full pass doesn't evict the
        scans the user is looking at.
        """
        idx = self.scan_index()
        ms1 = idx[idx["ms_level"] == 1]
        out = np.zeros(len(ms1))
        lo_mz, hi_mz = mz - tolerance, mz + tolerance
        for i, scan in enumerate(ms1["scan"]):
            spec = self._read_spectrum(int(scan))
            lo = np.searchsorted(spec.mz, lo_mz, side="left")
            hi = np.searchsorted(spec.mz, hi_mz, side="right")
            out[i] = spec.intensity[lo:hi].sum()
        return ms1["rt"], out, ms1["scan"]

    def method_info(self):
        """Run and method information (msview.method.MethodInfo).

        The base version summarises the scan index; backends add whatever
        their format stores (sample, instrument, method text, pump logs).
        """
        from msview.method import MethodInfo, summary_rows

        info = MethodInfo()
        for group, key, value in summary_rows(self.scan_index()):
            info.rows.append(("Acquisition summary", group, key, value))
        return info

    def spectrum(self, scan: int) -> Spectrum:
        with self._cache_lock:
            hit = self._cache.get(scan)
            if hit is not None:
                self._cache.move_to_end(scan)
                return hit
        spec = self._read_spectrum(scan)
        with self._cache_lock:
            self._cache[scan] = spec
            while len(self._cache) > self.cache_size:
                self._cache.popitem(last=False)
        return spec

    def average(self, scans: Sequence[int], ppm: float = 5.0) -> Spectrum:
        from msview.averaging import average_spectra

        return average_spectra([self.spectrum(int(s)) for s in scans], ppm=ppm)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
