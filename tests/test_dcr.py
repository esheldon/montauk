import numpy as np
import galsim
import mimsim


def test_dcr():

    detnum = 130
    seed = 9191

    gs_rng = galsim.BaseDeviate(seed)

    obsdata = mimsim.simtools.load_example_obsdata(band='g')

    dm_detector = mimsim.camera.make_dm_detector(detnum)

    wcs, icrf_to_field = mimsim.wcs.make_batoid_wcs(
        obsdata=obsdata, dm_detector=dm_detector,
    )

    dcr_maker = mimsim.dcr.DCRMaker(
        bandpass=obsdata['bandpass'],
        hour_angle=obsdata['HA'],
    )

    image_pos = galsim.PositionD(81.98, 55.23)
    sky_pos = wcs.toWorld(image_pos)

    dcr = dcr_maker(sky_pos)

    expected = galsim.PhotonDCR(
        base_wavelength=obsdata['bandpass'].effective_wavelength,
        latitude=mimsim.utils.get_latitude(),
        HA=obsdata['HA'],
        obj_coord=sky_pos,
    )
    assert dcr == expected

    # DCR is differential to effective_wavelength, so offset
    pa = galsim.PhotonArray(
        N=100, wavelength=obsdata['bandpass'].effective_wavelength * 1.2,
    )
    local_wcs = wcs.local(image_pos=image_pos)

    # no need for samplers here

    xold = pa.x.copy()
    yold = pa.y.copy()
    dcr.applyTo(pa, local_wcs=local_wcs)

    assert np.allclose(pa.x - xold, 1.0143603223428197)
    assert np.allclose(pa.y - yold, 0.23943545583254222)


if __name__ == '__main__':
    test_dcr()
