import galsim
import mimsim


def test_psfws():
    import psfws.simulate_atm_psf

    seed = 1991

    # default is 800, use 100 for speed
    psf_config = {'screen_size': 100}

    obsdata = mimsim.simtools.load_example_obsdata()

    psf = mimsim.psfws.make_psfws_psf(
        obsdata=obsdata,
        gs_rng=galsim.BaseDeviate(seed),
        psf_config=psf_config,
    )

    expected_psf = psfws.simulate_atm_psf.AtmosphericPSF(
        rng=galsim.BaseDeviate(seed),
        alt=obsdata['altitude'].deg,
        az=obsdata['azimuth'].deg,
        band=obsdata['band'],
        boresight=obsdata['boresight'],
        rawSeeing=obsdata['rawSeeing'],
        exptime=obsdata['exptime'],
        **psf_config
    )

    image_pos = galsim.PositionD(25.2, 100.8)

    psf_at_pos = psf.getPSF(image_pos)
    expected_psf_at_pos = expected_psf.getPSF(image_pos)
    assert psf_at_pos == expected_psf_at_pos

    eval_psf = mimsim.psf.eval_psf(psf, image_pos)
    assert eval_psf == expected_psf_at_pos
