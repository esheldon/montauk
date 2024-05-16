def extract_boresight(obsdata):
    """
    Get the boresight from the opsim data as a galsim.CelestialCoord
    """
    import galsim
    return galsim.CelestialCoord(
        obsdata['rightascension'] * galsim.degrees,
        obsdata['declination'] * galsim.degrees,
    )


def get_latitude():
    """
    Get the latitude of the site as a galsim.Angle
    """
    import galsim
    from lsst.obs.lsst.translators.lsst import SIMONYI_LOCATION as RUBIN_LOC

    return RUBIN_LOC.lat.deg * galsim.degrees
