def make_sensor(dm_detector, gs_rng, tree_rings=None):
    """
    Make a galsim.SiliconSensor object

    Parameters
    ----------
    dm_detector: lsst.afw.cameraGeom.Detector
        Data management detector object.  Use make_dm_detector(detnum)
    gs_rng: galsim rng
        The galsim rng
    tree_rings: imsim.treerings.TreeRings, optional
        Optional tree rings

    Returns
    -------
    galsim.SiliconSensor
    """
    import galsim

    detname = dm_detector.getName()

    kw = {
        'strength': 1.0e-6,
        'rng': gs_rng,
    }
    if tree_rings is not None:
        center, func = tree_rings.info[detname]
        kw['treering_func'] = func
        kw['treering_center'] = center

    serial_number = dm_detector.getSerial()

    if serial_number[:3] == 'E2V':
        kw['name'] = 'lsst_e2v_50_8'
    else:
        kw['name'] = 'lsst_itl_50_8'

    return galsim.SiliconSensor(**kw)


def make_tree_rings(dm_detectors):
    """
    Make a imsim.treerinsg.TreeRings instance

    Parameters
    -----------
    dm_detectors: [lsst.afw.cameraGeom.Detector]
        List of data management detector objects.  Use make_dm_detector(detnum)

    Returns
    -------
    imsim.treerings.TreeRings
    """
    import imsim

    names = [det.getName() for det in dm_detectors]
    return imsim.treerings.TreeRings(
        'tree_ring_parameters_2018-04-26.txt',
        only_dets=names,
        defer_load=False,
    )
