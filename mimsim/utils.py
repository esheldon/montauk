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


def get_initial_draw_method(flux):
    """
    Determine if we should use an fft for a bright object.  We may override
    this in get_fft_psf_maybe

    Returns
    -------
    draw method: str
    """
    from .defaults import FFT_FLUX_THRESH, FFT_SB_THRESH

    # this is from imsim.stamp, and the default fft_sb_thresh
    # in the example configs is 2.0e5 while FFT_FLUX_THRESH was
    # always 1.e6
    if flux < FFT_FLUX_THRESH or flux < FFT_SB_THRESH:
        method = 'phot'
    else:
        method = 'fft'

    return method
