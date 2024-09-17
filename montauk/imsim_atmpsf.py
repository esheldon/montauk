def make_imsim_atmpsf(obsdata, gs_rng, psf_config={}):
    """
    Make an imsim AtmosphericPSF

    Parameters
    ----------
    obsdata: dict
        Information about this observation, must have keys
            altitude, azimuth, band, boresight, rawSeeing, exptime
    gs_rng: galsim rng
        The galsim random number generator
    psf_config: dict, optional
        Config for psf, can have t0, nlayers, screen_size, exponent,
        kcrit, screen_scale, save_file, field_x, field_y

    Returns
    --------
    psfws.simulate_atm_psf.AtmosphericPSF
    """
    import imsim.atmPSF

    return imsim.atmPSF.AtmosphericPSF(
        airmass=obsdata['airmass'],
        rawSeeing=obsdata['rawSeeing'],
        band=obsdata['band'],
        boresight=obsdata['boresight'],
        rng=gs_rng,
        exptime=obsdata['exptime'],
        **psf_config
    )
