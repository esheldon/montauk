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
        Config for psf, can have t0, kkcrit, screen_size, screen_scale,
        exponent, nproc

    Returns
    --------
    ImSimAtmosphericPSFWrapper
    """
    import imsim.atmPSF

    atm_psf = imsim.atmPSF.AtmosphericPSF(
        airmass=obsdata['airmass'],
        rawSeeing=obsdata['rawSeeing'],
        band=obsdata['band'],
        boresight=obsdata['boresight'],
        rng=gs_rng,
        exptime=obsdata['exptime'],
        **psf_config
    )
    return ImSimAtmosphericPSFWrapper(atm_psf)


class ImSimAtmosphericPSFWrapper:
    """
    Wrapper that always does second kick
    """
    def __init__(self, atm_psf):
        self.atm_psf = atm_psf

    def getPSF(self, field_pos, gsparams=None):  # noqa
        import galsim

        atm_psf = self.atm_psf

        theta = (field_pos.x * galsim.arcsec, field_pos.y * galsim.arcsec)
        return galsim.ChromaticAtmosphere(
            atm_psf.atm.makePSF(
                atm_psf.wlen_eff,
                aper=atm_psf.aper,
                theta=theta,
                t0=atm_psf.t0,
                exptime=atm_psf.exptime,
                gsparams=gsparams,
            ),
            alpha=atm_psf.exponent,
            base_wavelength=atm_psf.wlen_eff,
            # Turns off DCR, since we apply that later using PhotonDCR.
            zenith_angle=0 * galsim.degrees
        )
