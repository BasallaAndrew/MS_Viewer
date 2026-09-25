"""mzML backend (fallback for files converted with ThermoRawFileParser/msconvert).

A small purpose-built reader rather than a general XML parser: it finds each
<spectrum> by byte offset (from the indexedmzML index, or by scanning the file
if there is none), parses only the metadata before <binaryDataArrayList> to
build the scan index, and decodes arrays on demand. os.pread keeps reads
thread-safe without a shared file position.
"""

from __future__ import annotations

import base64
import mmap
import os
import re
import zlib

import numpy as np

from msview import index as ix
from msview.data.base import RawSource, Spectrum

_CV = re.compile(rb"<cvParam\b[^>]*?>")
_ATTR = re.compile(rb'(\w+)="([^"]*)"')
_USER = re.compile(rb'<userParam\b[^>]*?name="([^"]*)"[^>]*?value="([^"]*)"')
_SPECTRUM_TAG = re.compile(rb"<spectrum\b[^>]*>")
_OFFSET = re.compile(rb'<offset\s+idRef="([^"]*)"[^>]*>(\d+)</offset>')
_ARRAY = re.compile(rb"<binaryDataArray\b.*?</binaryDataArray>", re.S)
_BINARY = re.compile(rb"<binary>([^<]*)</binary>")
_SCAN_NO = re.compile(rb"scan=(\d+)")

MS_LEVEL = b"MS:1000511"
CENTROID = b"MS:1000127"
SCAN_START = b"MS:1000016"
TIC = b"MS:1000285"
BASE_PEAK_INT = b"MS:1000505"
SELECTED_MZ = b"MS:1000744"
CHARGE = b"MS:1000041"
ISO_TARGET = b"MS:1000827"
ISO_LOWER = b"MS:1000828"
ISO_UPPER = b"MS:1000829"
COLLISION_ENERGY = b"MS:1000045"
SECONDS = (b"UO:0000010",)

ACTIVATIONS = {
    b"MS:1000133": "CID",
    b"MS:1000422": "HCD",
    b"MS:1002481": "HCD",
    b"MS:1000598": "ETD",
    b"MS:1000250": "ECD",
    b"MS:1000262": "IRMPD",
    b"MS:1000435": "PD",
}

DTYPES = {b"MS:1000523": "<f8", b"MS:1000521": "<f4", b"MS:1000522": "<i8", b"MS:1000519": "<i4"}
ZLIB = b"MS:1000574"
NUMPRESS = (b"MS:1002312", b"MS:1002313", b"MS:1002314")
MZ_ARRAY = b"MS:1000514"
INTENSITY_ARRAY = b"MS:1000515"


def _cv_params(block: bytes) -> dict[bytes, tuple[bytes, bytes]]:
    """accession -> (value, unitAccession), first occurrence wins."""
    out: dict[bytes, tuple[bytes, bytes]] = {}
    for tag in _CV.findall(block):
        attrs = dict(_ATTR.findall(tag))
        acc = attrs.get(b"accession")
        if acc and acc not in out:
            out[acc] = (attrs.get(b"value", b""), attrs.get(b"unitAccession", b""))
    return out


def _float(cv, acc, default=np.nan) -> float:
    try:
        return float(cv[acc][0])
    except (KeyError, ValueError):
        return default


class MzMLSource(RawSource):
    def __init__(self, path: str):
        super().__init__(path)
        self._fd = os.open(path, os.O_RDONLY)
        self._size = os.fstat(self._fd).st_size
        self._offsets = self._spectrum_offsets()
        self._scan_to_pos: dict[int, int] = {}

    def close(self) -> None:
        if self._fd is not None:
            os.close(self._fd)
            self._fd = None

    # --- locating spectra ---------------------------------------------------

    def _spectrum_offsets(self) -> np.ndarray:
        offsets = self._offsets_from_index()
        if offsets is None:
            offsets = self._offsets_by_scanning()
        return offsets

    def _offsets_from_index(self) -> np.ndarray | None:
        tail = os.pread(self._fd, 4096, max(0, self._size - 4096))
        m = re.search(rb"<indexListOffset>(\d+)</indexListOffset>", tail)
        if not m:
            return None
        start = int(m.group(1))
        block = os.pread(self._fd, self._size - start, start)
        s = block.find(b'<index name="spectrum"')
        if s < 0:
            return None
        e = block.find(b"</index>", s)
        offsets = np.array([int(o) for _, o in _OFFSET.findall(block[s:e])], dtype=np.int64)
        # Guard against a stale or broken index.
        if not len(offsets) or not os.pread(self._fd, 9, int(offsets[0])).startswith(b"<spectrum"):
            return None
        return offsets

    def _offsets_by_scanning(self) -> np.ndarray:
        with open(self.path, "rb") as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            return np.array([m.start() for m in re.finditer(rb"<spectrum\s", mm)], dtype=np.int64)

    def _read_until(self, offset: int, marker: bytes, chunk: int = 16384) -> bytes:
        """Bytes from offset up to (not including) marker, or to </spectrum>."""
        data = b""
        while True:
            data += os.pread(self._fd, chunk, offset + len(data))
            for end_marker in (marker, b"</spectrum>"):
                end = data.find(end_marker)
                if end >= 0:
                    return data[:end]
            if offset + len(data) >= self._size:
                return data
            chunk *= 2

    # --- index --------------------------------------------------------------

    def _build_index(self) -> np.ndarray:
        idx = ix.empty_index(len(self._offsets))
        scan_numbers = []
        for i, offset in enumerate(self._offsets):
            header = self._read_until(int(offset), b"<binaryDataArrayList")
            tag = _SPECTRUM_TAG.search(header)
            spectrum_attrs = dict(_ATTR.findall(tag.group(0))) if tag else {}
            m = _SCAN_NO.search(spectrum_attrs.get(b"id", b""))
            scan_numbers.append(int(m.group(1)) if m else i + 1)
            self._fill_row(idx[i], header)
        if len(set(scan_numbers)) != len(scan_numbers):
            scan_numbers = list(range(1, len(self._offsets) + 1))
        idx["scan"] = scan_numbers
        self._scan_to_pos = {s: p for p, s in enumerate(scan_numbers)}
        self._fill_missing_totals(idx)
        return ix.finalize(idx)

    @staticmethod
    def _fill_row(row, header: bytes) -> None:
        cv = _cv_params(header)
        level = int(_float(cv, MS_LEVEL, 1))
        row["ms_level"] = level
        row["centroid"] = CENTROID in cv
        rt = _float(cv, SCAN_START, 0.0)
        if SCAN_START in cv and cv[SCAN_START][1] in SECONDS:
            rt /= 60.0
        row["rt"] = rt
        row["tic"] = _float(cv, TIC)
        row["base_peak"] = _float(cv, BASE_PEAK_INT)
        if level < 2:
            return
        precursor = _float(cv, SELECTED_MZ)
        for name, value in _USER.findall(header):
            if b"Monoisotopic M/Z" in name:
                try:
                    mono = float(value)
                    if mono > 0:
                        precursor = mono
                except ValueError:
                    pass
        row["precursor_mz"] = precursor
        row["charge"] = int(_float(cv, CHARGE, 0))
        target = _float(cv, ISO_TARGET)
        lower, upper = _float(cv, ISO_LOWER), _float(cv, ISO_UPPER)
        if not np.isnan(lower) and not np.isnan(upper):
            row["isolation_width"] = lower + upper
            if not np.isnan(target):
                row["isolation_mz"] = target + (upper - lower) / 2
        elif not np.isnan(target):
            row["isolation_mz"] = target
        row["collision_energy"] = _float(cv, COLLISION_ENERGY)
        names = [ACTIVATIONS[a] for a in ACTIVATIONS if a in cv]
        if "ETD" in names and ("HCD" in names or "CID" in names):
            row["activation"] = "EThcD" if "HCD" in names else "ETciD"
        elif names:
            row["activation"] = names[0]

    def _fill_missing_totals(self, idx: np.ndarray) -> None:
        """Compute TIC/base peak from the data for writers that omit them."""
        missing = np.flatnonzero(np.isnan(idx["tic"]) | np.isnan(idx["base_peak"]))
        for i in missing:
            spec = self._decode(int(self._offsets[i]))
            total = float(spec.intensity.sum()) if len(spec) else 0.0
            peak = float(spec.intensity.max()) if len(spec) else 0.0
            if np.isnan(idx["tic"][i]):
                idx["tic"][i] = total
            if np.isnan(idx["base_peak"][i]):
                idx["base_peak"][i] = peak

    # --- spectra ------------------------------------------------------------

    def _read_spectrum(self, scan: int) -> Spectrum:
        if not self._scan_to_pos:
            self.scan_index()
        return self._decode(int(self._offsets[self._scan_to_pos[scan]]))

    def _decode(self, offset: int) -> Spectrum:
        element = self._read_until(offset, b"</spectrum>", chunk=65536)
        centroid = CENTROID in _cv_params(element[: element.find(b"<binaryDataArrayList")])
        mz = intensity = None
        for block in _ARRAY.findall(element):
            header, _, _ = block.partition(b"<binary>")
            cv = _cv_params(header)
            if any(n in cv for n in NUMPRESS):
                raise ValueError("MS-Numpress compressed mzML is not supported; re-convert without it")
            dtype = next((DTYPES[a] for a in DTYPES if a in cv), "<f8")
            m = _BINARY.search(block)
            raw = base64.b64decode(m.group(1)) if m else b""
            if ZLIB in cv and raw:
                raw = zlib.decompress(raw)
            values = np.frombuffer(raw, dtype=dtype).astype(np.float64)
            if MZ_ARRAY in cv:
                mz = values
            elif INTENSITY_ARRAY in cv:
                intensity = values
        if mz is None or intensity is None or len(mz) != len(intensity):
            return Spectrum.empty(centroid)
        if len(mz) > 1 and np.any(np.diff(mz) < 0):
            order = np.argsort(mz, kind="stable")
            mz, intensity = mz[order], intensity[order]
        return Spectrum(mz, intensity, centroid)
