import galsim
import montauk
import numpy as np
import pytest


@pytest.mark.parametrize('cosmic_ray_rate', [0, 10])
def test_cosmic_rays(cosmic_ray_rate):

    seed = 1982

    detnum = 79
    dm_detector = montauk.camera.make_dm_detector(detnum)
    obsdata = montauk.simtools.load_example_obsdata()

    cosmics = montauk.cosmic_rays.CosmicRays(
        cosmic_ray_rate=cosmic_ray_rate,
        exptime=obsdata['vistime'],
        gs_rng=galsim.BaseDeviate(seed),
    )

    wcs, _ = montauk.wcs.make_batoid_wcs(
        obsdata=obsdata, dm_detector=dm_detector,
    )

    bbox = dm_detector.getBBox()
    image = galsim.ImageD(bbox.width, bbox.height, wcs=wcs)
    image.array[:, :] = 0

    cosmics.add(image)
    w = np.where(image.array != 0)
    if cosmic_ray_rate == 10:
        assert w[0].size == 4268
    elif cosmic_ray_rate == 0:
        assert w[0].size == 0
    else:
        raise ValueError(f'got cosmic_ray_rate {cosmic_ray_rate} != 0 or 10')


if __name__ == '__main__':
    test_cosmic_rays(10)
