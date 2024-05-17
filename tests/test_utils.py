import mimsim
import mimsim.simtools
import galsim


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
