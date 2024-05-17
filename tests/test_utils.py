import mimsim
import galsim
import pytest


def test_latitude():
    from lsst.obs.lsst.translators.lsst import SIMONYI_LOCATION as RUBIN_LOC
    lat = mimsim.utils.get_latitude()

    expected = RUBIN_LOC.lat.deg * galsim.degrees

    assert lat == expected


def test_extract_boresight():
    obsdata = mimsim.simtools.load_example_obsdata()

    boresight = mimsim.utils.extract_boresight(obsdata)
    expected = galsim.CelestialCoord(
        obsdata['rightascension'] * galsim.degrees,
        obsdata['declination'] * galsim.degrees,
    )

    assert boresight == expected


@pytest.mark.parametrize('flux', [100, 1.e7])
def test_initial_draw_method(flux):
    from mimsim.defaults import FFT_FLUX_THRESH, FFT_SB_THRESH
    draw_method = mimsim.utils.get_initial_draw_method(flux)

    if flux < FFT_FLUX_THRESH or flux < FFT_SB_THRESH:
        assert draw_method == 'phot'
    else:
        assert draw_method == 'fft'
