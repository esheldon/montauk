import galsim
import imsim
import mimsim
import mimsim.simtools
import numpy as np


def test_gradient():

    detnum = 88
    dm_detector = mimsim.camera.make_dm_detector(detnum)
    obsdata = mimsim.simtools.load_example_obsdata()

    wcs, _ = mimsim.wcs.make_batoid_wcs(
        obsdata=obsdata, dm_detector=dm_detector,
    )

    sky_model = imsim.sky_model.SkyModel(
        exptime=obsdata['exptime'],
        mjd=obsdata['mjd'],
        bandpass=obsdata['bandpass'],
    )

    gradient = mimsim.sky.FixedSkyGradient(sky_model)

    bbox = dm_detector.getBBox()
    image = galsim.ImageD(bbox.width, bbox.height, wcs=wcs)
    image.array[:, :] = 1

    gradient.apply(image)

    ny, nx = image.array.shape

    world_center = wcs.toWorld(image.true_center)

    sky_gradient = imsim.sky_model.SkyGradient(
        sky_model=sky_model,
        wcs=image.wcs,
        world_center=world_center,
        image_xsize=nx,
    )
    xarr, yarr = np.meshgrid(range(nx), range(ny))
    assert np.all(image.array[:, :] == sky_gradient(xarr, yarr))


if __name__ == '__main__':
    test_gradient()
