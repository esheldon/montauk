import galsim
import montauk
import pytest


@pytest.mark.parametrize('psf_type', ['psfws', galsim.Gaussian])
def test_eval_psf(psf_type):

    seed = 1991

    # default is 800, use 100 for speed
    psf_config = {'screen_size': 100}

    obsdata = montauk.simtools.load_example_obsdata()
    image_pos = galsim.PositionD(25.2, 100.8)

    if psf_type == 'psfws':
        psf = montauk.psfws.make_psfws_psf(
            obsdata=obsdata,
            gs_rng=galsim.BaseDeviate(seed),
            psf_config=psf_config,
        )
        expected_psf_at_pos = psf.getPSF(image_pos)
    else:
        psf = psf_type(half_light_radius=1)
        expected_psf_at_pos = psf

    psf_at_pos = montauk.psf.eval_psf(psf, image_pos)
    assert psf_at_pos == expected_psf_at_pos


@pytest.mark.parametrize('psf_type', ['psfws', 'gauss', 'moffat'])
def test_fixed_psf(psf_type):

    seed = 99

    # default is 800, use 100 for speed
    if psf_type == 'psfws':
        psf_options = {'screen_size': 100}
    else:
        psf_options = {'fwhm': 1.0}
        if psf_type == 'moffat':
            psf_options['beta'] = 3.0

    obsdata = montauk.simtools.load_example_obsdata()
    image_pos = galsim.PositionD(25.2, 100.8)

    if psf_type == 'psfws':
        psf = montauk.psfws.make_psfws_psf(
            obsdata=obsdata,
            gs_rng=galsim.BaseDeviate(seed),
            psf_config=psf_options,
        )
        expected_psf_at_pos = psf.getPSF(image_pos)
    else:
        psf = montauk.fixed_psf.make_fixed_psf(
            type=psf_type, options=psf_options,
        )
        expected_psf_at_pos = psf

    psf_at_pos = montauk.psf.eval_psf(psf, image_pos)
    assert psf_at_pos == expected_psf_at_pos

    if psf_type != 'psfws':
        if psf_type == 'gauss':
            assert isinstance(psf, galsim.Gaussian)
        elif psf_type == 'moffat':
            assert isinstance(psf, galsim.Moffat)
        else:
            raise RuntimeError(f'unexpected psf type {psf_type}')
