def make_detector_telescope(band, rot_tel_pos, dm_detector):
    """
    get "detector telescope" as a batoid.Optic (batoid.CompoundOptic)

    Parameters
    ----------
    band: str
        ugriz
    rot_tel_pos: galsim.Angle
        Rotation angle for telescope
    dm_detector: lsst.afw.cameraGeom.Detector
        Data management detector object.  Use make_dm_detector(detnum)

    Returns
    -------
    batoid.Optic
        Most likely the batoid.CompountOptic subclass
    """
    import imsim

    telname = f'LSST_{band}.yaml'
    det_tel_obj = imsim.telescope_loader.DetectorTelescope(
        telname,
        rotTelPos=rot_tel_pos,
    )

    detname = dm_detector.getName()

    z_offset = det_tel_obj.calculate_z_offset(detname)
    return det_tel_obj.get_telescope(z_offset)


def make_pupil_sampler():
    """
    get a galsim PupulAnnulusSampler for Rubin
    """
    import galsim

    return galsim.PupilAnnulusSampler(
        R_outer=4.18,
        # M1 inner diameter is 2.558, but we need a bit of slack for
        # off-axis rays
        R_inner=2.55,
    )
