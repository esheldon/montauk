def make_batoid_wcs(
    obsdata,
    dm_detector,
    order=3,
    temperature=280,
    pressure=72.7,
    H2O_pressure=1.0,

):
    """
    Parameters
    ----------
    obsdata: dict
        Information about this observation
    dm_detector: lsst.afw.cameraGeom.Detector
        Data management detector object.  Use make_dm_detector(detnum)
    order: int, optional
        Order of the SIP polynomial [default: 3]
    temperature : float, optional
        Ambient temperature in Kelvin [default: 280 K]
    pressure : float, optional
        Ambient pressure in kPa [default: based on LSST heigh of 2715 meters]
    H2O_pressure : float, optional
        Water vapor pressure in kPa [default: 1.0 kPa]

    Returns
    --------
    wcs, icrf_to_field
        Both are galsim.FittedSIPWCS
    """
    import astropy.time
    import imsim
    from .telescope import make_detector_telescope

    obstime = astropy.time.Time(
        obsdata['mjd'], format='mjd', scale='tai',
    )
    wavelength = obsdata['bandpass'].effective_wavelength

    detector_telescope = make_detector_telescope(
        band=obsdata['band'],
        rot_tel_pos=obsdata['rotTelPos'],
        dm_detector=dm_detector,
    )

    factory = imsim.batoid_wcs.BatoidWCSFactory(
        boresight=obsdata['boresight'],
        obstime=obstime,
        telescope=detector_telescope,
        wavelength=wavelength,
        camera='LsstCam',
        temperature=temperature,
        pressure=pressure,
        H2O_pressure=H2O_pressure,
    )
    return (
        factory.getWCS(det=dm_detector, order=order),
        factory.get_icrf_to_field(det=dm_detector, order=order),
    )


def get_pixel_scale(wcs, bbox):
    """
    get the pixel scale at the center of the image

    Parameters
    ----------
    wcs: galsim.FitsWCS
        The galsim wcs
    bbox: lsst.geom.Box2I
        The image bounding box
    """
    import numpy as np
    import galsim

    x = bbox.getWidth() / 2
    y = bbox.getHeight() / 2
    return np.sqrt(wcs.pixelArea(galsim.PositionD(x, y)))
