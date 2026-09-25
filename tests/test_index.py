import numpy as np

from msview import index as ix


def make_index():
    # scans: 1 MS1, 2-3 MS2, 4 MS1, 5-7 MS2, 8 MS1
    levels = [1, 2, 2, 1, 2, 2, 2, 1]
    idx = ix.empty_index(len(levels))
    idx["scan"] = np.arange(1, len(levels) + 1)[::-1]  # deliberately unsorted
    idx["ms_level"] = levels[::-1]
    idx["rt"] = np.linspace(0, 7, len(levels))[::-1]
    return ix.finalize(idx)


def test_finalize_sorts_and_sets_parents():
    idx = make_index()
    assert idx["scan"].tolist() == list(range(1, 9))
    assert idx["parent_ms1"].tolist() == [1, 1, 1, 4, 4, 4, 4, 8]


def test_msn_before_first_ms1_has_no_parent():
    idx = ix.empty_index(3)
    idx["scan"] = [1, 2, 3]
    idx["ms_level"] = [2, 1, 2]
    idx = ix.finalize(idx)
    assert idx["parent_ms1"].tolist() == [-1, 2, 2]


def test_ms2_for_cycle():
    idx = make_index()
    assert ix.ms2_for_cycle(idx, 1).tolist() == [2, 3]
    assert ix.ms2_for_cycle(idx, 4).tolist() == [5, 6, 7]
    assert ix.ms2_for_cycle(idx, 8).tolist() == []


def test_nearest_ms1():
    idx = make_index()  # rt == scan - 1
    assert ix.nearest_ms1(idx, 0.2) == 1
    assert ix.nearest_ms1(idx, 1.4) == 1
    assert ix.nearest_ms1(idx, 1.6) == 4
    assert ix.nearest_ms1(idx, 100) == 8


def test_step_ms1():
    idx = make_index()
    assert ix.step_ms1(idx, 1, 1) == 4
    assert ix.step_ms1(idx, 4, -1) == 1
    assert ix.step_ms1(idx, 8, 1) == 8
    assert ix.step_ms1(idx, 1, -1) == 1
    # From an MS2 (cycle of scan 4), step to the nearest MS1 in that direction.
    assert ix.step_ms1(idx, 5, 1) == 8
    assert ix.step_ms1(idx, 5, -1) == 4


def test_ranges():
    idx = make_index()
    assert ix.ms1_in_range(idx, 0.5, 6.5).tolist() == [4]
    assert ix.ms2_in_range(idx, 0.5, 6.5).tolist() == [2, 3, 5, 6, 7]


def test_evenly_subsample():
    scans = np.arange(1000)
    sub = ix.evenly_subsample(scans, 500)
    assert len(sub) == 500 and sub[0] == 0 and sub[-1] == 999
    assert ix.evenly_subsample(scans[:10], 500).tolist() == list(range(10))


def test_isolation_window_and_precursor_fallback():
    idx = ix.empty_index(2)
    idx["scan"] = [1, 2]
    idx["ms_level"] = [1, 2]
    idx["isolation_mz"][1] = 500.0
    idx["isolation_width"][1] = 2.0
    idx = ix.finalize(idx)
    assert idx["precursor_mz"][1] == 500.0
    assert ix.isolation_window(idx[1]) == (499.0, 501.0)
    assert ix.isolation_window(idx[0]) is None


def precursor_index():
    # scan: 1 MS1 | 2 MS2 of 500.0 z2 | 3 MS2 of 500.6 z2 | 4 MS1 | 5 MS2 of 500.0 z2 | 6 MS2 of 800 z1
    idx = ix.empty_index(6)
    idx["scan"] = np.arange(1, 7)
    idx["ms_level"] = [1, 2, 2, 1, 2, 2]
    idx["rt"] = [10.0, 10.01, 10.02, 10.5, 10.51, 30.0]
    idx["precursor_mz"][[1, 2, 4, 5]] = [500.0, 500.6, 500.0, 800.0]
    idx["charge"][[1, 2, 4, 5]] = [2, 2, 2, 1]
    idx["isolation_width"][[1, 2, 4, 5]] = 2.0
    return ix.finalize(idx)


def test_ms2_for_precursor_prefers_exact_and_nearest_rt():
    idx = precursor_index()
    scans, pos = ix.ms2_for_precursor(idx, 500.0, rt=10.45, rt_window=1.0)
    assert scans.tolist() == [2, 5]  # 500.6 is inside the window but not a match
    assert scans[pos] == 5


def test_ms2_for_precursor_isotope_and_window_fallback():
    idx = precursor_index()
    scans, _ = ix.ms2_for_precursor(idx, 500.0 + 1.003355 / 2, rt=10.0, rt_window=1.0)
    assert scans.tolist() == [2, 5]  # clicked the M+1 isotope of the z=2 precursor
    scans, pos = ix.ms2_for_precursor(idx, 500.3, rt=10.0, rt_window=1.0)
    assert scans.tolist() == [2, 3, 5] and scans[pos] == 2  # only window containment


def test_ms2_for_precursor_respects_rt_window():
    idx = precursor_index()
    scans, pos = ix.ms2_for_precursor(idx, 800.0, rt=10.0, rt_window=1.0)
    assert len(scans) == 0 and pos == -1
    scans, pos = ix.ms2_for_precursor(idx, 800.0, rt=29.5, rt_window=1.0)
    assert scans.tolist() == [6]
