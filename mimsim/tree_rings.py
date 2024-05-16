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
