import galsim
import mimsim
import pytest


@pytest.mark.parametrize('psf_type', ['psfws', galsim.Gaussian])
def test_eval_psf(psf_type):

    seed = 1991

    # default is 800, use 100 for speed
    psf_config = {'screen_size': 100}

    obsdata = mimsim.simtools.load_example_obsdata()
    image_pos = galsim.PositionD(25.2, 100.8)

    if psf_type == 'psfws':
        psf = mimsim.psfws.make_psfws_psf(
            obsdata=obsdata,
            gs_rng=galsim.BaseDeviate(seed),
            psf_config=psf_config,
        )
        expected_psf_at_pos = psf.getPSF(image_pos)
    else:
        psf = psf_type(half_light_radius=1)
        expected_psf_at_pos = psf

    psf_at_pos = mimsim.psf.eval_psf(psf, image_pos)
    assert psf_at_pos == expected_psf_at_pos
