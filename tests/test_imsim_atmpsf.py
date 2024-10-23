import galsim
import montauk


def test_imsim_atmpsf():
    import imsim.atmPSF

    seed = 1991

    # default is 800, use 100 for speed
    psf_config = {'screen_size': 100}

    obsdata = montauk.simtools.load_example_obsdata()

    psf = montauk.imsim_atmpsf.make_imsim_atmpsf(
        obsdata=obsdata,
        gs_rng=galsim.BaseDeviate(seed),
        psf_config=psf_config,
    )

    epsf = imsim.atmPSF.AtmosphericPSF(
        airmass=obsdata['airmass'],
        rawSeeing=obsdata['rawSeeing'],
        band=obsdata['band'],
        boresight=obsdata['boresight'],
        rng=galsim.BaseDeviate(seed),
        exptime=obsdata['exptime'],
        **psf_config
    )
    expected_psf = montauk.imsim_atmpsf.ImSimAtmosphericPSFWrapper(epsf)

    image_pos = galsim.PositionD(25.2, 100.8)

    psf_at_pos = psf.getPSF(image_pos)
    expected_psf_at_pos = expected_psf.getPSF(image_pos)

    assert psf_at_pos == expected_psf_at_pos

    eval_psf = montauk.psf.eval_psf(psf, image_pos)
    assert eval_psf == expected_psf_at_pos


if __name__ == '__main__':
    test_imsim_atmpsf()
