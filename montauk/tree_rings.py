def make_tree_rings(detnums):
    """
    Make a imsim.treerinsg.TreeRings instance

    Parameters
    -----------
    detnums: sequence
        The detector numbers, e.g. [3, 25]

    Returns
    -------
    imsim.treerings.TreeRings
    """
    import imsim
    from .camera import make_dm_detector
    names = []
    for detnum in detnums:
        det = make_dm_detector(detnum)
        names.append(det.getName())

    return imsim.treerings.TreeRings(
        'tree_ring_parameters_2018-04-26.txt',
        only_dets=names,
        defer_load=False,
    )
