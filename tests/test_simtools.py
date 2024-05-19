import mimsim
import numpy as np
import pytest


@pytest.mark.parametrize('band', ['i', 'Y'])
def test_load_example_obsdata(band):
    obsdata = mimsim.simtools.load_example_obsdata(band)
    assert obsdata['band'] == band or obsdata['band'] == band.lower()


@pytest.mark.parametrize('band', ['i', 'Y'])
def test_load_example_instcat(band):
    seed = 8991
    rng = np.random.default_rng(seed)

    cat = mimsim.simtools.load_example_instcat(rng=rng, band=band)
    # this means we go the random ra, dec within the detector region
    assert cat.getNObjects() == 5


@pytest.mark.parametrize('n', [None, 3])
def test_ccd_radec_generator(n):
    seed = 3113
    rng = np.random.default_rng(seed)

    dm_detector = mimsim.camera.make_dm_detector(35)
    obsdata = mimsim.simtools.load_example_obsdata(band='i')

    wcs, _ = mimsim.wcs.make_batoid_wcs(
        obsdata=obsdata, dm_detector=dm_detector,
    )

    radec_gen = mimsim.simtools.CCDRadecGenerator(rng=rng, wcs=wcs)

    ra, dec = radec_gen(n)
    if n is None:
        assert isinstance(ra, np.float64)
    else:
        assert hasattr(ra, 'size')
        assert ra.size == n


if __name__ == '__main__':
    test_load_example_obsdata()
    test_load_example_instcat()
    test_ccd_radec_generator(None)
