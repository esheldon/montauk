import galsim
import montauk
from montauk.fringing import should_apply_fringing
from imsim.sky_model import CCD_Fringing
import numpy as np
import pytest


@pytest.mark.parametrize('band', ['r', 'Y'])
@pytest.mark.parametrize('detnum', [35, 88])
def test_fringing(band, detnum):

    obsdata = montauk.simtools.load_example_obsdata(band=band)

    dm_detector = montauk.camera.make_dm_detector(detnum)

    wcs, _ = montauk.wcs.make_batoid_wcs(
        obsdata=obsdata, dm_detector=dm_detector,
    )

    serial_number = dm_detector.getSerial()
    if serial_number[:3] == 'E2V' and band == 'Y':
        assert should_apply_fringing(band=band, dm_detector=dm_detector)

        fringer = montauk.fringing.Fringer(
            boresight=obsdata['boresight'], dm_detector=dm_detector,
        )

        bbox = dm_detector.getBBox()
        image = galsim.ImageD(bbox.width, bbox.height, wcs=wcs)
        image.array[:, :] = 1

        fringer.apply(image)

        # make sure we did it right
        skycen = image.wcs.toWorld(image.true_center)
        ccd_fringing = CCD_Fringing(
            true_center=skycen,
            boresight=obsdata['boresight'],
            seed=fringer.seed,
            spatial_vary=True,
        )
        ny, nx = image.array.shape
        xarr, yarr = np.meshgrid(range(nx), range(ny))
        fringing_map = ccd_fringing.calculate_fringe_amplitude(xarr, yarr)

        assert np.all(image.array[:, :] == fringing_map)

    else:
        assert not should_apply_fringing(band=band, dm_detector=dm_detector)


if __name__ == '__main__':
    test_fringing()
