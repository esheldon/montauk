from .defaults import MAX_STAMP_SIZE, PIXEL_SCALE


def get_stamp_size(obj, flux, noise_var, obsdata):
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

    Returns
    -------
    stamp_size: int
    """
    import imsim

    obj_achrom = obj.evaluateAtWavelength(
        obsdata['bandpass'].effective_wavelength
    )
    return imsim.stamp_utils.get_stamp_size(
        obj_achrom=obj_achrom,
        nominal_flux=flux,
        noise_var=noise_var,
        airmass=obsdata['airmass'],
        rawSeeing=obsdata['rawSeeing'],
        band=obsdata['band'],
        Nmax=MAX_STAMP_SIZE,
        pixel_scale=PIXEL_SCALE,
    )
