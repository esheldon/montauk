import galsim
import mimsim
import pytest


@pytest.mark.parametrize('use_dcr', [False, True])
def test_photon_ops_maker(use_dcr):

    detnum = 91
    dm_detector = mimsim.camera.make_dm_detector(detnum)
    obsdata = mimsim.simtools.load_example_obsdata()

    wcs, icrf_to_field = mimsim.wcs.make_batoid_wcs(
        obsdata=obsdata, dm_detector=dm_detector,
    )

    image_pos = galsim.PositionD(250.1, 1091.3)
    sky_pos = wcs.toWorld(image_pos)

    if use_dcr:
        dcr = mimsim.dcr.DCRMaker(
            bandpass=obsdata['bandpass'],
            hour_angle=obsdata['HA'],
        )
    else:
        dcr = None

    optics = mimsim.optics.OpticsMaker(
        altitude=obsdata['altitude'],
        azimuth=obsdata['azimuth'],
        boresight=obsdata['boresight'],
        rot_tel_pos=obsdata['rotTelPos'],
        band=obsdata['band'],
        dm_detector=dm_detector,
        wcs=wcs,
        icrf_to_field=icrf_to_field,
    )

    photon_ops_maker = mimsim.photon_ops.PhotonOpsMaker(
        exptime=obsdata['exptime'],
        band=obsdata['band'],
        dcr=dcr,
        optics=optics,
    )
    ops = photon_ops_maker(sky_pos)

    expected = [
        galsim.TimeSampler(t0=0.0, exptime=obsdata['exptime']),
        mimsim.telescope.make_pupil_sampler(),
    ]

    if use_dcr:
        expected += [dcr(sky_pos)]

    expected += [
        optics(sky_pos),
        mimsim.optics.make_focus_depth(obsdata['band']),
        mimsim.optics.make_refraction(),
    ]

    assert len(ops) == len(expected)
    for op, expected_op in zip(ops, expected):
        assert op == expected_op


if __name__ == '__main__':
    test_photon_ops_maker()
