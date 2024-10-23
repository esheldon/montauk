def make_psfws_psf(obsdata, gs_rng, psf_config={}):
    """
    Make a psfws AtmosphericPSF

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
    import psfws.simulate_atm_psf

    return psfws.simulate_atm_psf.AtmosphericPSF(
        rng=gs_rng,
        # alt=obsdata['altitude'].deg,
        # az=obsdata['azimuth'].deg,
        # new psfws wants units
        alt=obsdata['altitude'],
        az=obsdata['azimuth'],
        band=obsdata['band'],
        boresight=obsdata['boresight'],
        rawSeeing=obsdata['rawSeeing'],
        exptime=obsdata['exptime'],
        **psf_config
    )
