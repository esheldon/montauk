import galsim
import imsim
import mimsim
import numpy as np
import pytest
import logging


@pytest.mark.parametrize(
    'options',
    [{'gradient': False,
      'fringing': False,
      'vignetting': False,
      'sensor_areas': False},
     {'gradient': True,
      'fringing': False,
      'vignetting': False,
      'sensor_areas': False},
     {'gradient': False,
      'fringing': True,
      'vignetting': False,
      'sensor_areas': False},
     {'gradient': False,
      'fringing': False,
      'vignetting': True,
      'sensor_areas': False},
     {'gradient': False,
      'fringing': False,
      'vignetting': False,
      'sensor_areas': True},
     {'gradient': True,
      'fringing': True,
      'vignetting': True,
      'sensor_areas': True}]
)
def test_make_sky(options):
    # 88 is E2V
    band = 'Y'
    detnum = 88

    seed = 31415
    gs_rng = galsim.BaseDeviate(seed)

    logger = logging.getLogger('imsim-runner')

    dm_detector = mimsim.camera.make_dm_detector(detnum)
    obsdata = mimsim.simtools.load_example_obsdata()

    wcs, _ = mimsim.wcs.make_batoid_wcs(
        obsdata=obsdata, dm_detector=dm_detector,
    )

    sky_model = imsim.SkyModel(
        exptime=obsdata['exptime'],
        mjd=obsdata['mjd'],
        bandpass=obsdata['bandpass'],
    )

    if options['gradient']:
        gradient = mimsim.sky.FixedSkyGradient(sky_model)
    else:
        gradient = None

    if options['vignetting']:
        vignetter = mimsim.vignetting.Vignetter(dm_detector)
    else:
        vignetter = None

    if options['fringing']:
        assert mimsim.fringing.should_apply_fringing(
            band=band, dm_detector=dm_detector,
        )
        fringer = mimsim.fringing.Fringer(
            boresight=obsdata['boresight'], dm_detector=dm_detector,
        )
    else:
        fringer = None

    if options['sensor_areas']:
        sensor = mimsim.sensor.make_sensor(
            dm_detector=dm_detector, gs_rng=gs_rng,
        )
    else:
        sensor = None

    bbox = dm_detector.getBBox()
    nx = bbox.width
    ny = bbox.height
    image = mimsim.sky.make_sky_image(
        sky_model=sky_model,
        wcs=wcs,
        nx=nx,
        ny=ny,
        gradient=gradient,
        vignetting=vignetter,
        fringing=fringer,
        sensor=sensor,
        logger=logger,
    )

    test_image = image.copy()
    test_image.array[:, :] = 1

    x = nx / 2
    y = ny / 2
    image_pos = galsim.PositionD(x=x, y=y)
    sky_pos = wcs.toWorld(image_pos)
    sky_level = sky_model.get_sky_level(sky_pos)

    # values will vary due to wcs pixel size variations
    wcs.makeSkyImage(test_image, sky_level)

    print('sky_level:', sky_level)
    print('med sky image:', np.median(test_image.array[:, :]))

    if options['gradient']:
        gradient.apply(test_image)

    if options['vignetting']:
        vignetter.apply(test_image)

    if options['fringing']:
        fringer.apply(test_image)

    if options['sensor_areas']:
        area = sensor.calculate_pixel_areas(test_image)
        test_image *= area

    assert np.all(test_image.array[:, :] == image.array[:, :])


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
    test_make_sky({
        'gradient': False,
        'fringing': False,
        'vignetting': False,
        'sensor_areas': False,
    })
