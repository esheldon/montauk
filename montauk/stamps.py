from .defaults import PIXEL_SCALE


def get_stamp_size(obj, flux, noise_var, obsdata, pixel_scale=PIXEL_SCALE):
    """
    Get a good stamp size for the object

    Parameters
    ----------
    obj: GSObject or Sum
        The object of interest, before PSF
    flux: float
        The flux of the object
    noise_var: float
        Variance of background
    obsdata: dict
        Information about this observation
    pixel_scale: float, optional
        Pixel scale, default mimsim.defaults.PIXEL_SCALE

    Returns
    -------
    stamp_size: int
    """
    import imsim
    from .defaults import MIN_STAMP_SIZE, MAX_STAMP_SIZE

    obj_achrom = obj.evaluateAtWavelength(
        obsdata['bandpass'].effective_wavelength
    )
    stamp_size = imsim.stamp_utils.get_stamp_size(
        obj_achrom=obj_achrom,
        nominal_flux=flux,
        noise_var=noise_var,
        airmass=obsdata['airmass'],
        rawSeeing=obsdata['rawSeeing'],
        band=obsdata['band'],
        Nmax=MAX_STAMP_SIZE,
        pixel_scale=pixel_scale,
    )

    if stamp_size < MIN_STAMP_SIZE:
        stamp_size = MIN_STAMP_SIZE

    return stamp_size


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
