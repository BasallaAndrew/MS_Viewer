"""Data backends. Use `open_source(path)`; the GUI never imports a backend directly."""

from __future__ import annotations

from pathlib import Path

from msview.data.base import BackendUnavailable, RawSource, Spectrum

SUPPORTED_SUFFIXES = (".raw", ".mzml")

__all__ = ["BackendUnavailable", "RawSource", "Spectrum", "open_source", "SUPPORTED_SUFFIXES"]


def open_source(path: str | Path) -> RawSource:
    path = Path(path).expanduser()
    suffix = path.suffix.lower()
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    if suffix == ".raw":
        if path.is_dir():
            raise ValueError(
                f"{path.name} is a folder; that's a Waters/Agilent .raw, not Thermo. "
                "Only Thermo .raw files are supported."
            )
        from msview.data.thermo import ThermoSource

        return ThermoSource(str(path))
    if suffix == ".mzml":
        from msview.data.mzml import MzMLSource

        return MzMLSource(str(path))
    raise ValueError(f"Unsupported file type: {path.suffix} (expected .raw or .mzML)")
