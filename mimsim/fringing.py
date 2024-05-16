def should_apply_fringing(band, dm_detector):
    """
    Determine if we should apply fringing

    Parameters
    ----------
    band: str
        If not set to Y, do not apply fringing
    dm_detector: lsst.afw.cameraGeom.Detector
        Data management detector object.  Use make_dm_detector(detnum)
        If it is not an E2V, do not apply fringing
    """
    from .camera import get_ccd_vendor

    if band == 'Y' and get_ccd_vendor(dm_detector) == 'E2V':
        return True
    else:
        return False


class Fringer(object):
    """
    A class to apply fringing.  Only makes sens for Y band E2V sensors

    Parameters
    ----------
    boresight: galsim.CelestialCoord
        Boresight for observation
    dm_detector: lsst.afw.cameraGeom.Detector
        Data management detector object.  Use make_dm_detector(detnum)
    """
    def __init__(self, boresight, dm_detector):
        import hashlib

        self.boresight = boresight
        self.dm_detector = dm_detector

        # Use the hash value of the serial number as random seed number to
        # make sure the height map of the same sensor remains unchanged for
        # different exposures.
        serial_number = dm_detector.getSerial()

        # Note: the regular Python hash function is non-deterministic, which is
        # not good.
        # Instead we use hashlib.sha256, which is deterministic and convert
        # that to an integer.
        # https://stackoverflow.com/questions/27954892/deterministic-hashing-in-python-3
        # Only apply fringing to e2v sensors.

        self.seed = int(
            hashlib.sha256(serial_number.encode('UTF-8')).hexdigest(), 16
        ) & 0xFFFFFFFF

    def apply(self, image):
        """
        Apply fringing to the input image

        Parameters
        ----------
        image: galsim.Image
            The image to which we apply fringing
        """
        import numpy as np
        from imsim.sky_model import CCD_Fringing

        skycen = image.wcs.toWorld(image.true_center)
        ccd_fringing = CCD_Fringing(
            true_center=skycen,
            boresight=self.boresight,
            seed=self.seed,
            spatial_vary=True,
        )
        ny, nx = image.array.shape
        xarr, yarr = np.meshgrid(range(nx), range(ny))
        fringing_map = ccd_fringing.calculate_fringe_amplitude(xarr, yarr)
        image.array[:] *= fringing_map
