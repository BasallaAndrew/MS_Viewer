import numpy as np
import pytest

from msview.averaging import average_centroid, average_profile, average_spectra
from msview.data.base import Spectrum


def S(mz, inten, centroid=True):
    return Spectrum(np.asarray(mz, float), np.asarray(inten, float), centroid)


def test_centroid_bins_within_ppm_and_weights_mz():
    a = S([500.0000, 600.0], [100, 10])
    b = S([500.0010, 700.0], [300, 20])  # 2 ppm from 500.0000
    out = average_centroid([a, b], ppm=5)
    assert len(out) == 3
    assert out.mz[0] == pytest.approx((500.0 * 100 + 500.001 * 300) / 400)
    assert out.intensity.tolist() == pytest.approx([200.0, 5.0, 10.0])


def test_centroid_separates_beyond_ppm():
    out = average_centroid([S([500.0], [1]), S([500.01], [1])], ppm=5)  # 20 ppm apart
    assert len(out) == 2


def test_centroid_empty_inputs():
    assert len(average_centroid([])) == 0
    out = average_centroid([S([], []), S([400.0], [8])])
    assert out.intensity.tolist() == [4.0]  # divided by both scans


def test_profile_identical_scans_average_to_same():
    mz = np.linspace(400, 401, 101)
    inten = np.exp(-((mz - 400.5) ** 2) / 0.001)
    out = average_profile([S(mz, inten, False)] * 3)
    assert not out.centroid
    peak = np.argmax(out.intensity)
    assert out.mz[peak] == pytest.approx(400.5, abs=0.011)
    assert out.intensity.max() == pytest.approx(inten.max(), rel=1e-6)


def test_profile_does_not_bridge_wide_gaps():
    # Two peak regions with no zeros between them: the gap must not fill in.
    mz1 = np.r_[np.linspace(400, 400.02, 5), np.linspace(410, 410.02, 5)]
    inten = np.full(10, 100.0)
    shifted = mz1 + 0.001
    out = average_profile([S(mz1, inten, False), S(shifted, inten, False)])
    between = (out.mz > 401) & (out.mz < 409)
    assert not np.any(out.intensity[between] > 0)


def test_dispatch():
    assert average_spectra([S([1.0], [1.0], True)]).centroid
    assert not average_spectra([S([1.0, 2.0], [1.0, 1.0], False)] * 2).centroid
