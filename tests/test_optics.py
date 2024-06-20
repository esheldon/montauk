import galsim
import montauk


def test_focus_depth():
    for band in 'ugrizYy':
        focus_depth = montauk.optics.make_focus_depth(band)

        expected = galsim.FocusDepth(
            montauk.optics.get_focus_depth_value(band)
        )
        assert focus_depth == expected


def test_optics_maker():

    seed = 55
    gs_rng = galsim.BaseDeviate(seed)

    detnum = 91
    dm_detector = montauk.camera.make_dm_detector(detnum)
    obsdata = montauk.simtools.load_example_obsdata()

    wcs, icrf_to_field = montauk.wcs.make_batoid_wcs(
        obsdata=obsdata, dm_detector=dm_detector,
    )

    optics = montauk.optics.OpticsMaker(
        altitude=obsdata['altitude'],
        azimuth=obsdata['azimuth'],
        boresight=obsdata['boresight'],
        rot_tel_pos=obsdata['rotTelPos'],
        band=obsdata['band'],
        dm_detector=dm_detector,
        wcs=wcs,
        icrf_to_field=icrf_to_field,
    )

    image_pos = galsim.PositionD(250.1, 1091.3)
    sky_pos = wcs.toWorld(image_pos)

    rubin_diffraction_optics = optics(sky_pos)

    pa = galsim.PhotonArray(
        N=100, wavelength=obsdata['bandpass'].effective_wavelength,
    )
    local_wcs = wcs.local(image_pos=image_pos)

    time_sampler = galsim.TimeSampler(exptime=obsdata['vistime'])
    pupil_sampler = montauk.telescope.make_pupil_sampler()

    time_sampler.applyTo(pa, rng=gs_rng)
    assert pa.x.std() == 0

    # this modifies pupil_u, pupil_v but not x, y
    pupil_sampler.applyTo(pa, local_wcs=local_wcs, rng=gs_rng)
    assert pa.x.std() == 0

    ustd = pa.pupil_u.std()
    vstd = pa.pupil_v.std()
    assert ustd > 2 and vstd > 2

    # modifies x, y
    # note rng can be sent but is ignored
    rubin_diffraction_optics.applyTo(pa, local_wcs=local_wcs)

    # finite spread from optics
    assert pa.x.std() > 0.2

    # There is a bug in imsim.photon_ops.RubinOptics.applyTo
    # that modifies the pupil values in place
    # https://github.com/LSSTDESC/imSim/issues/475
    # assert pa.pupil_u.std() == ustd
    # assert pa.pupil_v.std() == vstd


if __name__ == '__main__':
    test_optics_maker()
