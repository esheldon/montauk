class OpticsMaker(object):
    """
    A class functor to create a RubinDiffractionOptics

    Parameters
    ----------
    altitude: galsim.Angle
        Altitude of telescope
    azimuth: galsim.Angle
        Azimuth of telescope
    boresight: galsim.CelestialCoord
        Boresight for observation
    rot_tel_pos: galsim.Angle
        Angle of telescope rotation
    band: str
        Band of observation, e.g. u, g, r, i, z
    dm_detector: lsst.afw.cameraGeom.Detector
        Data management detector object.  Use make_dm_detector(detnum)
    wcs: galsim.GSFitsWCS
        The wcs, probably crude
    icrf_to_field: galsim.GSFitsWCS
        WCS to convert to field coords

    Calling the object returns
    --------------------------
    RubinDiffractionOptics
    """

    def __init__(
        self,
        altitude,
        azimuth,
        boresight,
        rot_tel_pos,
        band,
        dm_detector,
        wcs,
        icrf_to_field,
    ):
        from .telescope import make_detector_telescope
        from .camera import make_camera

        self.altitude = altitude
        self.azimuth = azimuth
        self.boresight = boresight
        self.rot_tel_pos = rot_tel_pos
        self.band = band
        self.dm_detector = dm_detector
        self.wcs = wcs
        self.icrf_to_field = icrf_to_field

        self.detector_telescope = make_detector_telescope(
            band=band,
            rot_tel_pos=rot_tel_pos,
            dm_detector=dm_detector,
        )
        self.camera = make_camera()

    def __call__(self, obj_coord):
        import imsim
        from .utils import get_latitude

        rubin_diffraction = imsim.photon_ops.RubinDiffraction(
            telescope=self.detector_telescope,
            latitude=get_latitude().rad,
            altitude=self.altitude.rad,
            azimuth=self.azimuth.rad,
            sky_pos=obj_coord,
            icrf_to_field=self.icrf_to_field,
        )

        obj_image_pos = self.wcs.toImage(obj_coord)

        return imsim.photon_ops.RubinDiffractionOptics(
            telescope=self.detector_telescope,
            boresight=self.boresight,
            image_pos=obj_image_pos,
            icrf_to_field=self.icrf_to_field,
            det_name=self.dm_detector.getName(),
            camera=self.camera,
            rubin_diffraction=rubin_diffraction,
            sky_pos=obj_coord,
        )


def get_focus_depth_value(band):
    """
    get the focus depth value for the input band

    Parameters
    -----------
    band: str
        band, case insensitive

    Returns
    -------
    the focus depth
    """
    from .defaults import FOCUS_DEPTH_DICT
    return FOCUS_DEPTH_DICT[band.lower()]


def make_focus_depth(band):
    """
    Make a FocusDepth object

    Parameters
    ----------
    band: str
        The band

    Returns
    -------
    galsim.FocusDepth
    """
    import galsim

    focus_depth_val = get_focus_depth_value(band)
    return galsim.FocusDepth(focus_depth_val)


def make_refraction():
    """
    Create a galsim Refraction object

    TODO index_ratio was for Y (Josh) need for other bands

    Returns
    -------
    galsim.Refraction(
    """
    import galsim
    return galsim.Refraction(index_ratio=3.9)
