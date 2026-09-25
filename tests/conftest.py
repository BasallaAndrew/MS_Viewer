import base64
import zlib
from pathlib import Path

import numpy as np
import pytest

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def _array_xml(values, accession, dtype, compress):
    raw = np.asarray(values, dtype=dtype).tobytes()
    if compress:
        raw = zlib.compress(raw)
    encoded = base64.b64encode(raw).decode()
    type_acc = "MS:1000523" if dtype == "<f8" else "MS:1000521"
    comp = '<cvParam cvRef="MS" accession="MS:1000574" name="zlib compression" value=""/>' if compress else (
        '<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>'
    )
    return (
        f'<binaryDataArray encodedLength="{len(encoded)}">'
        f'<cvParam cvRef="MS" accession="{type_acc}" name="float" value=""/>{comp}'
        f'<cvParam cvRef="MS" accession="{accession}" name="array" value=""/>'
        f"<binary>{encoded}</binary></binaryDataArray>"
    )


def write_mzml(path: Path, scans: list[dict], indexed: bool = True, rt_unit: str = "minute") -> None:
    """Write a minimal mzML. Each scan dict: level, rt, mz, intensity, and for
    MS2: precursor, charge, iso (half-width), activation accession, ce."""
    unit = "UO:0000031" if rt_unit == "minute" else "UO:0000010"
    parts = ['<?xml version="1.0" encoding="utf-8"?>\n<indexedmzML><mzML><run id="t">'
             f'<spectrumList count="{len(scans)}">\n']
    offsets = []
    for i, s in enumerate(scans):
        body = []
        mode = "MS:1000127" if s.get("centroid", True) else "MS:1000128"
        body.append(f'<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="{s["level"]}"/>')
        body.append(f'<cvParam cvRef="MS" accession="{mode}" name="mode" value=""/>')
        if s.get("write_totals", True):
            inten = np.asarray(s["intensity"])
            body.append(f'<cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="{inten.sum()}"/>')
            body.append(f'<cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="{inten.max() if len(inten) else 0}"/>')
        body.append(
            f'<scanList count="1"><scan><cvParam cvRef="MS" accession="MS:1000016" name="scan start time" '
            f'value="{s["rt"]}" unitCvRef="UO" unitAccession="{unit}" unitName="t"/></scan></scanList>'
        )
        if s["level"] >= 2:
            hw = s.get("iso", 1.0)
            body.append(
                '<precursorList count="1"><precursor><isolationWindow>'
                f'<cvParam cvRef="MS" accession="MS:1000827" name="target" value="{s["precursor"]}"/>'
                f'<cvParam cvRef="MS" accession="MS:1000828" name="lower" value="{hw}"/>'
                f'<cvParam cvRef="MS" accession="MS:1000829" name="upper" value="{hw}"/>'
                '</isolationWindow><selectedIonList count="1"><selectedIon>'
                f'<cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="{s["precursor"]}"/>'
                f'<cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="{s.get("charge", 2)}"/>'
                "</selectedIon></selectedIonList><activation>"
                f'<cvParam cvRef="MS" accession="{s.get("activation", "MS:1000422")}" name="act" value=""/>'
                f'<cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="{s.get("ce", 30)}"/>'
                "</activation></precursor></precursorList>"
            )
        body.append('<binaryDataArrayList count="2">')
        body.append(_array_xml(s["mz"], "MS:1000514", "<f8", compress=True))
        body.append(_array_xml(s["intensity"], "MS:1000515", "<f4", compress=i % 2 == 0))
        body.append("</binaryDataArrayList>")
        scan_no = s.get("scan", i + 1)
        head = "".join(parts).encode()
        offsets.append(len(head))
        parts.append(
            f'<spectrum index="{i}" id="controllerType=0 controllerNumber=1 scan={scan_no}" '
            f'defaultArrayLength="{len(s["mz"])}">' + "".join(body) + "</spectrum>\n"
        )
    parts.append("</spectrumList></run></mzML>\n")
    if indexed:
        index_offset = len("".join(parts).encode())
        parts.append('<indexList count="1"><index name="spectrum">')
        for i, (s, off) in enumerate(zip(scans, offsets)):
            parts.append(f'<offset idRef="scan={s.get("scan", i + 1)}">{off}</offset>')
        parts.append(f"</index></indexList><indexListOffset>{index_offset}</indexListOffset>")
    parts.append("</indexedmzML>\n")
    path.write_bytes("".join(parts).encode())


def dda_scans(n_cycles=5, n_ms2=3):
    """A small DDA run: each cycle is one MS1 followed by n_ms2 MS2s."""
    scans = []
    rt = 0.0
    for c in range(n_cycles):
        scans.append(dict(level=1, rt=rt, mz=[400.0, 500.0 + c, 600.0], intensity=[10.0, 100.0 + c, 50.0]))
        for k in range(n_ms2):
            rt += 0.01
            scans.append(dict(level=2, rt=rt, precursor=500.0 + c + k * 0.1, charge=2 + k, iso=0.8,
                              mz=[150.0, 250.0], intensity=[5.0, 7.0]))
        rt += 0.05
    return scans


@pytest.fixture
def dda_mzml(tmp_path):
    path = tmp_path / "dda.mzML"
    write_mzml(path, dda_scans())
    return path
