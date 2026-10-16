import numpy as np
import pytest

from msview.deconv import DeconSettings, available, charge_ladder, deconvolve

PROTON = 1.007276467


def test_charge_ladder():
    z, mz = charge_ladder(12358.0, 1, 50, 600, 1300)
    assert z.tolist() == list(range(10, 21))
    assert mz[0] == pytest.approx((12358.0 + 10 * PROTON) / 10)
    assert np.all((mz >= 600) & (mz <= 1300))


def synthetic_spectrum(mass, charges, width=0.3):
    mz = np.arange(500.0, 3000.0, 0.02)
    intensity = np.zeros_like(mz)
    for z in charges:
        centre = (mass + z * PROTON) / z
        intensity += np.exp(-((z - np.mean(charges)) ** 2) / 8) * np.exp(-((mz - centre) ** 2) / (2 * (width / 2.355) ** 2))
    return mz, intensity


@pytest.mark.skipif(not available(), reason="UniDec not installed")
def test_deconvolve_synthetic_protein():
    mz, intensity = synthetic_spectrum(25_000.0, range(10, 22))
    result = deconvolve(mz, intensity, DeconSettings(z_max=40, mass_max=60_000, peak_fwhm=0.3))
    assert result.peaks, "no mass found"
    assert result.peaks[0].mass == pytest.approx(25_000.0, abs=2.0)
    assert result.intensity.max() == pytest.approx(100.0)
    assert len(result.fit) == len(result.data_mz)


@pytest.mark.skipif(not available(), reason="UniDec not installed")
def test_deconvolve_rejects_empty_range():
    mz, intensity = synthetic_spectrum(25_000.0, range(10, 22))
    with pytest.raises(ValueError):
        deconvolve(mz, intensity, DeconSettings(mz_min=4000, mz_max=5000))
