import numpy as np
import pytest

from conftest import DATA_DIR, dda_scans, write_mzml
from msview import index as ix
from msview.data import open_source


def test_index_from_dda_file(dda_mzml):
    with open_source(dda_mzml) as src:
        idx = src.scan_index()
    assert len(idx) == 20
    assert (idx["ms_level"] == 1).sum() == 5
    ms2 = idx[idx["ms_level"] == 2]
    assert ms2["activation"].tolist() == ["HCD"] * 15
    assert ms2["charge"][:3].tolist() == [2, 3, 4]
    assert ms2["isolation_width"][0] == pytest.approx(1.6)
    assert ms2["precursor_mz"][1] == pytest.approx(500.1)
    assert ix.ms2_for_cycle(idx, 5).tolist() == [6, 7, 8]


def test_spectrum_and_chromatogram(dda_mzml):
    with open_source(dda_mzml) as src:
        spec = src.spectrum(5)
        assert spec.mz.tolist() == [400.0, 501.0, 600.0]
        assert spec.intensity.tolist() == [10.0, 101.0, 50.0]
        rt, tic, scans = src.chromatogram("tic")
        assert scans.tolist() == [1, 5, 9, 13, 17]
        assert tic[0] == pytest.approx(160.0)
        _, bpc, _ = src.chromatogram("bpc")
        assert bpc.tolist() == [100, 101, 102, 103, 104]
        assert src.spectrum(5) is spec  # cached


def test_unindexed_file_and_seconds(tmp_path):
    path = tmp_path / "noindex.mzML"
    scans = dda_scans(2, 1)
    for s in scans:
        s["rt"] *= 60
    write_mzml(path, scans, indexed=False, rt_unit="second")
    with open_source(path) as src:
        idx = src.scan_index()
        assert len(idx) == 4
        assert idx["rt"][2] == pytest.approx(scans[2]["rt"] / 60)
        assert src.spectrum(4).mz.tolist() == [150.0, 250.0]


def test_missing_totals_are_computed(tmp_path):
    path = tmp_path / "nototals.mzML"
    scans = dda_scans(2, 1)
    for s in scans:
        s["write_totals"] = False
    write_mzml(path, scans)
    with open_source(path) as src:
        _, tic, _ = src.chromatogram()
        assert tic.tolist() == pytest.approx([160.0, 161.0])


def test_average_range(dda_mzml):
    with open_source(dda_mzml) as src:
        avg = src.average([1, 5])
    assert avg.centroid
    assert len(avg) == 4  # 400, 500, 501, 600
    assert avg.intensity[0] == pytest.approx(10.0)


def test_unsupported_suffix(tmp_path):
    p = tmp_path / "x.txt"
    p.write_text("")
    with pytest.raises(ValueError):
        open_source(p)


# --- optional checks against real files in data/ ---------------------------

SAMPLE_MZML = DATA_DIR / "022824_A_BSA_01.mzML"
SAMPLE_RAW = DATA_DIR / "Waters_mAB_10ug.raw"


@pytest.mark.skipif(not SAMPLE_MZML.exists(), reason="sample mzML not present")
def test_sample_mzml():
    with open_source(SAMPLE_MZML) as src:
        idx = src.scan_index()
        assert np.all(np.diff(idx["rt"]) >= 0)
        ms1 = ix.nearest_ms1(idx, 20.0)
        assert len(src.spectrum(ms1)) > 0
        assert len(ix.ms2_for_cycle(idx, ms1)) >= 0


@pytest.mark.skipif(not SAMPLE_RAW.exists(), reason="sample .raw not present")
def test_sample_raw():
    from msview.data.base import BackendUnavailable

    try:
        src = open_source(SAMPLE_RAW)
    except BackendUnavailable as e:
        pytest.skip(str(e))
    with src:
        idx = src.scan_index()
        ms2 = idx[idx["ms_level"] == 2]
        assert len(ms2) and np.all(ms2["parent_ms1"] > 0)
        assert np.isfinite(ms2["precursor_mz"]).all()
        spec = src.spectrum(int(ms2["scan"][0]))
        assert len(spec) > 0 and np.all(np.diff(spec.mz) >= 0)


def test_xic(dda_mzml):
    # MS1 cycle c has peaks at 400, 500 + c, 600 with intensities 10, 100 + c, 50.
    with open_source(dda_mzml) as src:
        rt, inten, scans = src.xic(502.0, 0.01)
        assert scans.tolist() == [1, 5, 9, 13, 17]
        assert inten.tolist() == [0, 0, 102, 0, 0]
        _, wide, _ = src.xic(502.0, 2.0)  # 500..504 inclusive
        assert wide.tolist() == [100, 101, 102, 103, 104]
        _, edge, _ = src.xic(400.5, 0.5)  # window edge exactly on the 400.0 peak
        assert edge.tolist() == [10] * 5
