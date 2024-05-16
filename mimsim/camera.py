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
