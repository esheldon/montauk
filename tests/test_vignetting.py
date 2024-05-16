import galsim
import imsim
import mimsim
import numpy as np


def test_vignetting():

    detnum = 88
    dm_detector = mimsim.camera.make_dm_detector(detnum)

    vignetter = mimsim.vignetting.Vignetter(dm_detector)
    assert isinstance(vignetter.vignetting, imsim.vignetting.Vignetting)

    bbox = dm_detector.getBBox()
    image = galsim.ImageD(bbox.width, bbox.height)
    image.array[:, :] = 1
    vignetter.apply(image)
    assert np.all(image.array[:, :] == vignetter.vals)


if __name__ == '__main__':
    test_vignetting()
