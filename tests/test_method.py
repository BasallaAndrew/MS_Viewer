import numpy as np
import pytest

from conftest import dda_scans, write_mzml
from msview import index as ix
from msview.data import open_source
from msview.method import is_percent_b, parse_gradient, parse_method_text, summary_rows

# Excerpts of real method text (Agilent 1260 pump, Q Exactive HF) from a .raw file.
AGILENT_PUMP = """Agilent1260 G1311 Quat Pump 1

Solvent A: 0.1% FA in water
Solvent B: 0.1% FA in ACN
Solvent C: 
Minimum pressure limit (bar): 0.0
Maximum pressure limit (bar): 600.0

Gradient program:
    Time           Flow Rate          Composition
      0.00(min)    0.300(ml/min)      A=98.0% B=2.0% C=0.0% D=0.0%
     30.00(min)    0.300(ml/min)      A=86.0% B=14.0% C=0.0% D=0.0%
     70.00(min)    0.300(ml/min)      A=60.0% B=40.0% C=0.0% D=0.0%
     74.00(min)    300.0(ul/min)      A=98.0% B=2.0% C=0.0% D=0.0%
"""

QEXACTIVE = """                            Method of Q Exactive HF

FULL MS / DD-MS² (TOPN)
General
Runtime                                                              0 to 80 min
Inclusion                                                                  ― 
Full MS
Resolution                                                            60,000 
AGC target                                                               3e6 
dd-MS² / dd-SIM
Resolution                                                            15,000 

TUNEFILES
General
Base Tunefile C:\\Xcalibur\\methods\\tune.mstune

CONTACT CLOSURE
General
Start in Closed True 

EXCLUSION LIST
    2 entries
     Mass Formula Species  CS Polarity Start   End Comment
    [m/z]     [M]         [z]          [min] [min] 
195.12274                     Positive             
"""


def test_parse_gradient_agilent():
    g = parse_gradient(AGILENT_PUMP)
    assert g["time"].tolist() == [0, 30, 70, 74]
    assert g["percent_b"].tolist() == [2, 14, 40, 2]
    assert g["flow"].tolist() == pytest.approx([0.3, 0.3, 0.3, 0.3])  # µL/min converted


def test_parse_gradient_script_style_and_none():
    text = "  0.000 [min]  Pump.%B.Value:  4.0 [%]\n 30.000 [min]  Pump.%B.Value: 35.0 [%]\n"
    g = parse_gradient(text)
    assert g["percent_b"].tolist() == [4, 35] and np.isnan(g["flow"]).all()
    assert parse_gradient("Injection volume (ul): 10.0") is None


def test_parse_agilent_method_text():
    rows = parse_method_text(AGILENT_PUMP)
    params = {k: (g, v) for g, k, v in rows}
    assert params["Solvent B"] == ("Agilent1260 G1311 Quat Pump 1", "0.1% FA in ACN")
    # An empty "Solvent C:" must not become a heading for the next rows.
    assert params["Maximum pressure limit (bar)"][0] == "Agilent1260 G1311 Quat Pump 1"
    assert "Solvent C" not in params and "Time" not in params  # empty value, table header


def test_parse_thermo_method_text():
    rows = parse_method_text(QEXACTIVE)
    found = {(g, k): v for g, k, v in rows}
    assert found[("FULL MS / DD-MS² (TOPN) › Full MS", "Resolution")] == "60,000"
    assert found[("FULL MS / DD-MS² (TOPN) › dd-MS² / dd-SIM", "Resolution")] == "15,000"
    assert found[("TUNEFILES › General", "Base Tunefile")] == "C:\\Xcalibur\\methods\\tune.mstune"
    assert found[("CONTACT CLOSURE › General", "Start in Closed")] == "True"
    assert not any(k == "Inclusion" for _, k, _ in rows)  # "―" means unset
    assert not any(k.startswith(("195", "Mass Formula")) for _, k, _ in rows)  # table rows


def test_is_percent_b():
    assert is_percent_b("Composition B (%)")
    assert is_percent_b("%B")
    assert not is_percent_b("Composition A (%)")
    assert not is_percent_b("Pressure (bar)")


def test_summary_rows():
    idx = ix.empty_index(5)
    idx["scan"] = [1, 2, 3, 4, 5]
    idx["ms_level"] = [1, 2, 2, 1, 2]
    idx["rt"] = [0.0, 0.01, 0.02, 0.05, 0.06]
    idx["activation"] = ["", "HCD", "HCD", "", "CID"]
    idx["collision_energy"] = [np.nan, 30, 30, np.nan, 35]
    idx["charge"] = [0, 2, 3, 0, 2]
    rows = {k: v for _, k, v in summary_rows(ix.finalize(idx))}
    assert rows["MS1 scans"] == "2" and rows["MSn scans"] == "3"
    assert rows["Median MS1 cycle time"] == "3.00 s"
    assert rows["MSn per MS1 cycle (mean / max)"] == "1.5 / 2"
    assert rows["Activation / energy"] == "HCD 30 (2), CID 35 (1)"
    assert rows["Precursor charge"] == "2 (2), 3 (1)"


def test_mzml_method_info(tmp_path):
    path = tmp_path / "with_pressure.mzML"
    chroms = [
        dict(id="TIC", time=[0.0, 1.0], values=[1.0, 2.0], accession="MS:1000235"),
        dict(id="Pump Pressure 1", time=[0.0, 0.5, 1.0], values=[3e7, 3.5e7, 4e7], unit_name="pascal"),
    ]
    write_mzml(path, dda_scans(2, 1), chromatograms=chroms)
    with open_source(path) as src:
        info = src.method_info()
        assert len(src.spectrum(1)) == 3  # spectra unaffected by the chromatogram list
    rows = {(s, k): v for s, _, k, v in info.rows}
    assert rows[("Instrument", "Model")] == "Orbitrap Exploris 480"
    assert rows[("Instrument", "Instrument serial number")] == "SN123"
    assert rows[("Instrument", "Analyzer")] == "orbitrap"
    assert rows[("Instrument", "Xcalibur")] == "4.5"
    assert rows[("Sample", "Acquired")] == "2026-01-02T03:04:05Z"
    assert rows[("Acquisition summary", "MS1 scans")] == "2"
    assert [t.name for t in info.traces] == ["Pump Pressure 1 (bar)"]  # TIC skipped
    assert info.traces[0].values.tolist() == pytest.approx([300, 350, 400])  # Pa -> bar
    assert info.gradient is None and info.method_text == []


SAMPLE_RAW = __import__("conftest").DATA_DIR / "Waters_mAB_10ug.raw"


@pytest.mark.skipif(not SAMPLE_RAW.exists(), reason="sample .raw not present")
def test_sample_raw_method_info():
    from msview.data.base import BackendUnavailable

    try:
        src = open_source(SAMPLE_RAW)
    except BackendUnavailable as e:
        pytest.skip(str(e))
    with src:
        info = src.method_info()
        spec = src.spectrum(100)  # device selection restored to MS afterwards
    assert len(spec) > 0
    assert info.gradient is not None and info.gradient["percent_b"].max() == 85
    assert info.percent_b_trace() is not None
    assert any("Pressure" in t.name for t in info.traces)
    assert any(name.startswith("Agilent1260 G1311") for name, _ in info.method_text)
    sections = {s for s, *_ in info.rows}
    assert {"Acquisition summary", "Sample", "Instrument", "Tune"} <= sections
    assert not any(s == "Problems reading the file" for s in sections)
