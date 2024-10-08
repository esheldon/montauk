def make_camera():
    """
    Create a lsst.obs.lsst.LsstCam

    Returns
    -------
    camera: lsst.obs.lsst.LsstCam
    """
    from imsim.camera import get_camera
    return get_camera('LsstCam')


def make_dm_detector(detnum):
    """
    Get the DM detector object

    Parameters
    ----------
    detnum: int
        Detector number

    Returns
    ---------
    dm_detector: lsst.afw.cameraGeom.Detector
    """

    dm_camera = make_camera()
    return dm_camera[detnum]


def get_ccd_vendor(dm_detector):
    """
    Get the vendor, E2V or ITL

    Parameters
    ----------
    dm_detector: lsst.afw.cameraGeom.Detector
        Data management detector object.  Use make_dm_detector(detnum)
        If it is not an E2V, do not apply fringing

    Returns
    -------
    vendor: str
    """
    serial_number = dm_detector.getSerial()
    return serial_number[:3]
