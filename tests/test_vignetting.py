import galsim
import imsim
import mimsim
import numpy as np


def test_vignetting():

    detnum = 88
    dm_detector = mimsim.camera.make_dm_detector(detnum)

    obsdata = mimsim.simtools.load_example_obsdata()

    wcs, _ = mimsim.wcs.make_batoid_wcs(
        obsdata=obsdata, dm_detector=dm_detector,
    )

    vignetter = mimsim.vignetting.Vignetter(dm_detector)
    assert isinstance(vignetter.vignetting, imsim.vignetting.Vignetting)

    bbox = dm_detector.getBBox()
    image = galsim.ImageD(bbox.width, bbox.height, wcs=wcs)
    image.array[:, :] = 1
    vignetter.apply(image)
    assert np.all(image.array[:, :] == vignetter.vals)

    px = 50
    py = 100
    image_pos = galsim.PositionD(px, py)
    sky_pos = wcs.toWorld(image_pos)
    val = vignetter.at_sky_coord(wcs=wcs, sky_pos=sky_pos)
    assert np.allclose(val, image.array[py, px])


if __name__ == '__main__':
    test_vignetting()
