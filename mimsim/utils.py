def extract_boresight(obsdata):
    """
    Get the boresight from the opsim data as a galsim.CelestialCoord
    """
    import galsim
    return galsim.CelestialCoord(
        obsdata['rightascension'] * galsim.degrees,
        obsdata['declination'] * galsim.degrees,
    )
