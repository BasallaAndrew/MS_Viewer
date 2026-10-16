"""Intact-mass deconvolution with UniDec (https://github.com/michaelmarty/UniDec).

UniDec is optional (`pip install unidec`). Its engine runs a native program
that reads and writes files, so every run gets its own temporary folder.
Runs are serialized: the engine isn't meant for concurrent use.
"""

from __future__ import annotations

import contextlib
import io
import os
import tempfile
import threading
import time
import warnings
from dataclasses import dataclass, field

import numpy as np

_run_lock = threading.Lock()


class UniDecUnavailable(RuntimeError):
    pass


@dataclass
class DeconSettings:
    z_min: int = 1
    z_max: int = 50
    mass_min: float = 5_000.0
    mass_max: float = 200_000.0
    mass_bin: float = 1.0  # Da
    peak_fwhm: float | None = None  # m/z; None = let UniDec estimate it
    mz_min: float | None = None  # restrict the input m/z range
    mz_max: float | None = None
    peak_window: float = 20.0  # Da; peaks closer than this are merged
    peak_threshold: float = 0.05  # fraction of the tallest mass peak


@dataclass
class DeconPeak:
    mass: float
    height: float  # % of the tallest peak


@dataclass
class DeconResult:
    mass: np.ndarray  # zero-charge mass axis (Da)
    intensity: np.ndarray  # % of max
    peaks: list[DeconPeak]
    data_mz: np.ndarray  # the m/z data UniDec fitted (after its processing)
    fit: np.ndarray  # UniDec's reconstructed spectrum on data_mz (same scale as data)
    data_intensity: np.ndarray
    peak_fwhm: float
    r_squared: float
    seconds: float
    log: str = field(repr=False, default="")


def _import_engine():
    # UniDec pulls in matplotlib; keep it off the GUI backends (we may be on a worker thread).
    os.environ.setdefault("MPLBACKEND", "Agg")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from unidec.engine import UniDec
    except ImportError as e:
        raise UniDecUnavailable(
            "UniDec is not installed in this environment. Install it with:\n"
            "    .venv/bin/pip install unidec"
        ) from e
    return UniDec


def available() -> bool:
    try:
        import importlib.util

        return importlib.util.find_spec("unidec") is not None
    except Exception:
        return False


def deconvolve(mz: np.ndarray, intensity: np.ndarray, settings: DeconSettings) -> DeconResult:
    """Run UniDec on one spectrum and return the zero-charge mass spectrum."""
    mz = np.asarray(mz, dtype=np.float64)
    intensity = np.asarray(intensity, dtype=np.float64)
    if settings.mz_min is not None or settings.mz_max is not None:
        lo = -np.inf if settings.mz_min is None else settings.mz_min
        hi = np.inf if settings.mz_max is None else settings.mz_max
        keep = (mz >= lo) & (mz <= hi)
        mz, intensity = mz[keep], intensity[keep]
    if len(mz) < 10 or not np.any(intensity > 0):
        raise ValueError("Not enough signal in the selected m/z range to deconvolve")
    if settings.z_min > settings.z_max or settings.mass_min >= settings.mass_max:
        raise ValueError("Check the charge and mass ranges (minimum must be below maximum)")

    UniDec = _import_engine()
    log = io.StringIO()
    start = time.perf_counter()
    with _run_lock, tempfile.TemporaryDirectory(prefix="msview_unidec_") as workdir, \
            contextlib.redirect_stdout(log), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        engine = UniDec(ignore_args=True)
        engine.pass_data_in(np.column_stack([mz, intensity]), dirname=workdir, fname="spectrum.txt", silent=True)
        config = engine.config
        config.startz, config.endz = int(settings.z_min), int(settings.z_max)
        config.masslb, config.massub = float(settings.mass_min), float(settings.mass_max)
        config.massbins = float(settings.mass_bin)
        config.peakwindow = float(settings.peak_window)
        config.peakthresh = float(settings.peak_threshold)
        engine.process_data(silent=True)
        if settings.peak_fwhm is None:
            engine.get_auto_peak_width()
        else:
            config.mzsig = float(settings.peak_fwhm)
        status = engine.run_unidec(silent=True)
        if status != 0:
            raise RuntimeError(f"UniDec failed (code {status}).\n\n{log.getvalue()[-2000:]}")
        engine.pick_peaks(calc_dscore=False)

        massdat = np.asarray(engine.data.massdat, dtype=np.float64)
        data2 = np.asarray(engine.data.data2, dtype=np.float64)
        fit = np.asarray(engine.data.fitdat, dtype=np.float64)
        peaks = sorted((DeconPeak(float(p.mass), float(p.height)) for p in engine.pks.peaks),
                       key=lambda p: -p.height)
        fwhm, r2 = float(config.mzsig), float(config.error)

    top = massdat[:, 1].max() if len(massdat) else 1.0
    return DeconResult(
        mass=massdat[:, 0], intensity=massdat[:, 1] / (top or 1.0) * 100.0, peaks=peaks,
        data_mz=data2[:, 0], data_intensity=data2[:, 1],
        fit=fit if len(fit) == len(data2) else np.full(len(data2), np.nan),
        peak_fwhm=fwhm, r_squared=r2, seconds=time.perf_counter() - start, log=log.getvalue(),
    )


def charge_ladder(mass: float, z_min: int, z_max: int, mz_lo: float, mz_hi: float,
                  adduct: float = 1.007276467) -> tuple[np.ndarray, np.ndarray]:
    """(z, m/z) of the charge states of `mass` that fall inside [mz_lo, mz_hi]."""
    z = np.arange(max(1, z_min), z_max + 1)
    mz = (mass + z * adduct) / z
    keep = (mz >= mz_lo) & (mz <= mz_hi)
    return z[keep], mz[keep]
