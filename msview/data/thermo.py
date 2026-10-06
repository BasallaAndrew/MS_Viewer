"""Thermo .raw backend: Thermo's RawFileReader .NET assemblies via pythonnet.

The DLLs can't be redistributed, so they're loaded from a user folder
(default ~/.msview/rawfilereader, override with MSVIEW_RAWFILEREADER_DIR).
"""

from __future__ import annotations

import ctypes
import os
import re
import threading
from pathlib import Path

import numpy as np

from msview import index as ix
from msview.data.base import BackendUnavailable, RawSource, Spectrum

DOWNLOAD_URL = "https://github.com/thermofisherlsms/RawFileReader"
DEFAULT_DLL_DIR = Path.home() / ".msview" / "rawfilereader"
# Load order matters: dependencies first.
ASSEMBLIES = ("OpenMcdf", "ThermoFisher.CommonCore.Data", "ThermoFisher.CommonCore.RawFileReader")
REQUIRED = ("ThermoFisher.CommonCore.Data", "ThermoFisher.CommonCore.RawFileReader")
DOTNET_ROOTS = (Path.home() / ".dotnet", Path("/usr/local/share/dotnet"), Path("/opt/homebrew/share/dotnet"))

ACTIVATION_NAMES = {
    "CollisionInducedDissociation": "CID",
    "HigherEnergyCollisionalDissociation": "HCD",
    "ElectronTransferDissociation": "ETD",
    "ElectronCaptureDissociation": "ECD",
    "MultiPhotonDissociation": "MPD",
    "UltraVioletPhotoDissociation": "UVPD",
    "NegativeElectronTransferDissociation": "NETD",
    "PQD": "PQD",
}

# Status-log fields that are codes/states rather than measurements.
_STATUS_NOISE = re.compile(r"not ready|buffer size|\b(state|mode|code)\b", re.I)

_clr_lock = threading.Lock()
_clr_ready = False


def dll_dir() -> Path:
    return Path(os.environ.get("MSVIEW_RAWFILEREADER_DIR", DEFAULT_DLL_DIR)).expanduser()


def ensure_runtime() -> None:
    """Load .NET and the RawFileReader assemblies once (idempotent)."""
    global _clr_ready
    with _clr_lock:
        if _clr_ready:
            return
        folder = dll_dir()
        missing = [n for n in REQUIRED if not (folder / f"{n}.dll").exists()]
        if missing:
            raise BackendUnavailable(
                f"Thermo RawFileReader DLLs not found in {folder}\n"
                f"(missing: {', '.join(m + '.dll' for m in missing)}).\n\n"
                f"Download them from {DOWNLOAD_URL} and copy the .NET 8 DLLs "
                "into that folder. See the README for step-by-step instructions."
            )
        try:
            import pythonnet
        except ImportError as e:
            raise BackendUnavailable(
                "pythonnet is not installed. Run: pip install 'msview[thermo]'"
            ) from e
        if "DOTNET_ROOT" not in os.environ:
            for root in DOTNET_ROOTS:
                if (root / "shared" / "Microsoft.NETCore.App").is_dir():
                    os.environ["DOTNET_ROOT"] = str(root)
                    break
        try:
            pythonnet.load("coreclr")
        except Exception as e:
            raise BackendUnavailable(
                "Could not start the .NET 8 runtime. Install it (see README) or set "
                f"DOTNET_ROOT to its folder.\n\nDetails: {e}"
            ) from e
        import clr

        for name in ASSEMBLIES:
            path = folder / f"{name}.dll"
            if path.exists():
                clr.AddReference(str(path))
        _clr_ready = True


def _to_numpy(net_array) -> np.ndarray:
    """Copy a .NET double[] into a NumPy array without per-element calls."""
    if net_array is None or net_array.Length == 0:
        return np.empty(0)
    from System.Runtime.InteropServices import GCHandle, GCHandleType

    n = net_array.Length
    handle = GCHandle.Alloc(net_array, GCHandleType.Pinned)
    try:
        address = handle.AddrOfPinnedObject().ToInt64()
        return np.frombuffer((ctypes.c_double * n).from_address(address), np.float64).copy()
    finally:
        handle.Free()


class ThermoSource(RawSource):
    def __init__(self, path: str):
        super().__init__(path)
        ensure_runtime()
        from ThermoFisher.CommonCore.Data.Business import Device
        from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter

        # RawFileReader objects aren't thread-safe; serialize all access.
        self._lock = threading.Lock()
        raw = RawFileReaderAdapter.FileFactory(str(path))
        if raw is None or raw.IsError:
            msg = raw.FileError.ErrorMessage if raw is not None else "unknown error"
            raise OSError(f"RawFileReader could not open {path}: {msg}")
        if raw.InAcquisition:
            raw.Dispose()
            raise OSError(f"{path} is still being acquired")
        raw.SelectInstrument(Device.MS, 1)
        self._raw = raw

    def close(self) -> None:
        with self._lock:
            if self._raw is not None:
                self._raw.Dispose()
                self._raw = None

    def _trailer_slots(self) -> tuple[int | None, int | None]:
        mono = charge = None
        for i, header in enumerate(self._raw.GetTrailerExtraHeaderInformation()):
            label = str(header.Label).strip().rstrip(":").lower()
            if label == "monoisotopic m/z":
                mono = i
            elif label == "charge state":
                charge = i
        return mono, charge

    def _build_index(self) -> np.ndarray:
        with self._lock:
            raw = self._raw
            first, last = raw.RunHeaderEx.FirstSpectrum, raw.RunHeaderEx.LastSpectrum
            mono_slot, charge_slot = self._trailer_slots()
            idx = ix.empty_index(max(0, last - first + 1))
            for i, scan in enumerate(range(first, last + 1)):
                row = idx[i]
                row["scan"] = scan
                row["rt"] = raw.RetentionTimeFromScanNumber(scan)
                order = int(raw.GetFilterForScanNumber(scan).MSOrder)
                # Negative orders are neutral loss/gain/parent scans; treat as MSn.
                level = order if order >= 1 else 2
                row["ms_level"] = level
                stats = raw.GetScanStatsForScanNumber(scan)
                row["tic"] = stats.TIC
                row["base_peak"] = stats.BasePeakIntensity
                row["centroid"] = stats.IsCentroidScan
                if level >= 2:
                    self._fill_precursor(row, scan, level, mono_slot, charge_slot)
        return ix.finalize(idx)

    def _fill_precursor(self, row, scan, level, mono_slot, charge_slot) -> None:
        raw = self._raw
        event = raw.GetScanEventForScanNumber(scan)
        n_reactions = event.MassCount
        if n_reactions == 0:
            return
        reaction = event.GetReaction(min(level - 2, n_reactions - 1))
        target = reaction.PrecursorMass
        row["isolation_mz"] = target + reaction.IsolationWidthOffset
        row["isolation_width"] = reaction.IsolationWidth
        row["collision_energy"] = reaction.CollisionEnergy
        activation = str(reaction.ActivationType)
        row["activation"] = ACTIVATION_NAMES.get(activation, activation[:8])
        precursor = target
        if mono_slot is not None:
            try:
                mono = float(raw.GetTrailerExtraValue(scan, mono_slot))
                if mono > 0:
                    precursor = mono
            except Exception:
                pass
        row["precursor_mz"] = precursor
        if charge_slot is not None:
            try:
                row["charge"] = int(raw.GetTrailerExtraValue(scan, charge_slot))
            except Exception:
                pass

    def _read_spectrum(self, scan: int) -> Spectrum:
        with self._lock:
            raw = self._raw
            stream = raw.GetCentroidStream(scan, False)
            if stream is not None and stream.Length > 0:
                return Spectrum(_to_numpy(stream.Masses), _to_numpy(stream.Intensities), True)
            stats = raw.GetScanStatsForScanNumber(scan)
            seg = raw.GetSegmentedScanFromScanNumber(scan, stats)
            return Spectrum(
                _to_numpy(seg.Positions), _to_numpy(seg.Intensities), bool(stats.IsCentroidScan)
            )

    # --- method / run information ------------------------------------------

    def method_info(self):
        from ThermoFisher.CommonCore.Data.Business import Device

        info = super().method_info()
        with self._lock:
            try:
                for name, part in (
                    ("sample", self._add_sample_info),
                    ("instrument", self._add_instrument_info),
                    ("tune data", self._add_tune_data),
                    ("instrument method", self._add_methods),
                    ("status logs", self._add_status_logs),
                ):
                    try:
                        part(info)
                    except Exception as e:  # one vendor quirk shouldn't hide the rest
                        info.rows.append(("Problems reading the file", "", name, str(e)))
            finally:
                self._raw.SelectInstrument(Device.MS, 1)
        return info

    def _add_sample_info(self, info) -> None:
        si, header = self._raw.SampleInformation, self._raw.FileHeader
        info.add("Sample", "Sample", [
            ("Sample name", si.SampleName), ("Sample ID", si.SampleId), ("Sample type", si.SampleType),
            ("Vial", si.Vial), ("Injection volume (µL)", f"{si.InjectionVolume:g}"),
            ("Dilution factor", f"{si.DilutionFactor:g}"), ("Comment", si.Comment),
            ("Sequence row", si.RowNumber),
        ])
        info.add("Sample", "User text", [(f"User text {i + 1}", t) for i, t in enumerate(si.UserText or [])])
        info.add("Sample", "File", [
            ("Acquired", header.CreationDate), ("Operator", header.WhoCreatedId),
            ("Description", header.FileDescription), ("Original path", si.RawFileName),
            ("Instrument method file", si.InstrumentMethodFile),
            ("Processing method file", si.ProcessingMethodFile),
        ])

    def _add_instrument_info(self, info) -> None:
        data, run = self._raw.GetInstrumentData(), self._raw.RunHeaderEx
        info.add("Instrument", "Mass spectrometer", [
            ("Model", data.Model or data.Name), ("Serial number", data.SerialNumber),
            ("Software version", data.SoftwareVersion), ("Hardware version", data.HardwareVersion),
            ("Mass range (all scans)", f"{run.LowMass:g} – {run.HighMass:g} m/z"),
            ("Acquisition time", f"{run.StartTime:.2f} – {run.EndTime:.2f} min"),
        ])
        try:
            names = list(self._raw.GetAllInstrumentNamesFromInstrumentMethod())
        except Exception:
            names = []
        info.add("Instrument", "Devices in method", [(f"Device {i + 1}", n) for i, n in enumerate(names)])

    def _add_tune_data(self, info) -> None:
        raw = self._raw
        count = raw.GetTuneDataCount()
        labels = [str(h.Label) for h in raw.GetTuneDataHeaderInformation()] if count else []
        for t in range(count):
            section = "Tune" if count == 1 else f"Tune {t + 1}"
            group = ""
            for label, value in zip(labels, raw.GetTuneDataValues(t, True).Values):
                label = label.strip().rstrip(":").strip()
                if label.startswith("==="):
                    group = label.strip("= ").rstrip(":")
                    continue
                info.add(section, group, [(label, value)])

    def _add_methods(self, info) -> None:
        from msview.method import parse_gradient, parse_method_text

        raw = self._raw
        try:
            names = list(raw.GetAllInstrumentNamesFromInstrumentMethod())
        except Exception:
            names = []
        for i in range(raw.InstrumentMethodsCount):
            text = str(raw.GetInstrumentMethod(i) or "")
            if not text.strip():
                continue
            name = str(names[i]) if i < len(names) else text.strip().splitlines()[0].strip()
            info.method_text.append((name, text))
            for group, key, value in parse_method_text(text):
                info.rows.append((f"Method: {name}", group, key, value))
            if info.gradient is None:
                info.gradient = parse_gradient(text)

    def _add_status_logs(self, info) -> None:
        """Numeric, non-constant status-log fields of every device as traces
        (LC devices first: pump %B, pressure and flow, oven temperature...)."""
        from ThermoFisher.CommonCore.Data.Business import Device

        raw = self._raw
        devices = [(Device.Other, i) for i in range(1, raw.GetInstrumentCountOfType(Device.Other) + 1)]
        devices += [(Device.MS, 1)]
        for device_type, number in devices:
            raw.SelectInstrument(device_type, number)
            name = str(raw.GetInstrumentData().Name or device_type)
            info.traces += self._status_log_traces(name)

    def _status_log_traces(self, device_name: str) -> list:
        from msview.method import Trace

        raw = self._raw
        n = raw.GetStatusLogEntriesCount()
        if n < 2:
            return []
        labels = [str(h.Label).strip().rstrip(":").strip() for h in raw.GetStatusLogHeaderInformation()]
        rt = np.empty(n)
        table = []
        for i in range(n):
            entry = raw.GetStatusLogValues(i, False)
            rt[i] = entry.RetentionTime
            table.append(list(entry.Values))
        order = np.argsort(rt, kind="stable")
        traces = []
        for j, label in enumerate(labels):
            values = np.full(n, np.nan)
            for i, row in enumerate(table):
                try:
                    values[i] = float(row[j])
                except (IndexError, ValueError, TypeError):
                    pass
            finite = np.isfinite(values)
            if (finite.mean() < 0.9 or np.nanmax(values) == np.nanmin(values)
                    or _STATUS_NOISE.search(label)):
                continue  # text fields, constants and status codes make dull plots
            traces.append(Trace(label, device_name, rt[order], values[order]))
        return traces
