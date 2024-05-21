def load_obsdata_from_instcat(instcat):
    """
    load observational data from an instcat header, applying data
    conversions

    Parameters
    ----------
    instcat: str
        The instcat

    Returns
    -------
    OpsimDataLoader
    """
    import imsim
    import galsim
    from .utils import extract_boresight

    obsdata = imsim.opsim_data.OpsimDataLoader(instcat).meta
    obsdata['boresight'] = extract_boresight(obsdata)
    obsdata['bandpass'] = imsim.bandpass.RubinBandpass(obsdata['band'])
    obsdata['rotTelPos'] = obsdata['rotTelPos'] * galsim.degrees
    obsdata['altitude'] = obsdata['altitude'] * galsim.degrees
    obsdata['azimuth'] = obsdata['azimuth'] * galsim.degrees
    obsdata['HA'] = obsdata['HA'] * galsim.degrees
    return obsdata
