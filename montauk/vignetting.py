class Vignetter(object):
    """
    A class to apply vignetting to a detector

    Parameters
    ---------
    dm_detector: lsst.afw.cameraGeom.Detector
        Data management detector object.  Use make_dm_detector(detnum)
    """
    def __init__(self, dm_detector):
        from imsim.vignetting import Vignetting
        from lsst.afw import cameraGeom
        from .camera import make_camera

        self.dm_detector = dm_detector
        detname = self.dm_detector.getName()

        camera = make_camera()
        radii = Vignetting.get_pixel_radii(camera[detname])

        self.vignetting = Vignetting('LSSTCam_vignetting_data.json')
        self.vals = self.vignetting.apply_to_radii(radii)

        self.pix_to_fp = dm_detector.getTransform(
            cameraGeom.PIXELS,
            cameraGeom.FOCAL_PLANE,
        )

    def apply(self, image):
        """
        Apply vignetting to a image

        Parameters
        ----------
        image: galsim.Image
            The sky image
        """
        image.array[:] *= self.vals

    def at_sky_coord(self, wcs, sky_pos):
        """
        evaluate vignetting at the input location

        Parameters
        ----------
        wcs: galsim wcs
            The wcs to convert sky position to image position
        sky_pos: galsim.CelestialCoord
            The sky location

        Returns
        -------
        vignetting: float
        """
        return self.vignetting.at_sky_coord(
            sky_pos, wcs, self.pix_to_fp,
        )
