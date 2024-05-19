import galsim
import imsim
import mimsim
import numpy as np
import pytest


@pytest.mark.parametrize(
    'options',
    [
        {'band': 'i',
         'limit': False,
         'select': False},
        {'band': 'i',
         'limit': True,
         'select': False},
        {'band': 'i',
         'limit': False,
         'select': True},
        {'band': 'Y',
         'limit': False,
         'select': False}
    ]
)
def test_runner_smoke(options):
    # 88 is E2V, which we want for fringing

    band = options['band']
    limit = options['limit']
    select = options['select']

    detnum = 88

    # default is 800, use 100 for speed
    psf_config = {'screen_size': 100}
    cosmic_ray_rate = mimsim.defaults.DEFAULT_COSMIC_RAY_RATE

    seed = 31415
    rng = np.random.default_rng(seed)
    gs_rng = galsim.BaseDeviate(rng.integers(0, 2**60))

    dm_detector = mimsim.camera.make_dm_detector(detnum)
    obsdata = mimsim.simtools.load_example_obsdata(band=band)

    wcs, icrf_to_field = mimsim.wcs.make_batoid_wcs(
        obsdata=obsdata, dm_detector=dm_detector,
    )

    cat = mimsim.simtools.load_example_instcat(
        rng=rng, band=band, detnum=detnum,
    )

    if limit:
        nobj = 2
        indices = np.arange(2)
    else:
        nobj = cat.getNObjects()
        indices = None

    def example_selector(d, i):
        return d.magnorm[i] > magmax

    def default_selector(d, i):
        return True

    if select:
        magmax = 21
        selector = example_selector
        select_num = np.sum(cat.magnorm > magmax)
    else:
        selector = default_selector
        select_num = nobj

    sky_model = imsim.SkyModel(
        exptime=obsdata['exptime'],
        mjd=obsdata['mjd'],
        bandpass=obsdata['bandpass'],
    )

    gradient = mimsim.sky.FixedSkyGradient(sky_model)
    vignetter = mimsim.vignetting.Vignetter(dm_detector)
    if band == 'Y':
        assert mimsim.fringing.should_apply_fringing(
            band=band, dm_detector=dm_detector,
        )
        fringer = mimsim.fringing.Fringer(
            boresight=obsdata['boresight'], dm_detector=dm_detector,
        )
    else:
        assert not mimsim.fringing.should_apply_fringing(
            band=band, dm_detector=dm_detector,
        )
        fringer = None

    sensor = mimsim.sensor.make_sensor(
        dm_detector=dm_detector, gs_rng=gs_rng,
    )

    dcr = mimsim.dcr.DCRMaker(
        bandpass=obsdata['bandpass'],
        hour_angle=obsdata['HA'],
    )

    optics = mimsim.optics.OpticsMaker(
        altitude=obsdata['altitude'],
        azimuth=obsdata['azimuth'],
        boresight=obsdata['boresight'],
        rot_tel_pos=obsdata['rotTelPos'],
        band=obsdata['band'],
        dm_detector=dm_detector,
        wcs=wcs,
        icrf_to_field=icrf_to_field,
    )

    photon_ops_maker = mimsim.photon_ops.PhotonOpsMaker(
        exptime=obsdata['exptime'],
        band=obsdata['band'],
        dcr=dcr,
        optics=optics,
    )

    diffraction_fft = imsim.stamp.DiffractionFFT(
        exptime=obsdata['exptime'],
        altitude=obsdata['altitude'],
        azimuth=obsdata['azimuth'],
        rotTelPos=obsdata['rotTelPos'],
    )

    artist = mimsim.artist.Artist(
        bandpass=obsdata['bandpass'],
        sensor=sensor,
        photon_ops_maker=photon_ops_maker,
        diffraction_fft=diffraction_fft,
        gs_rng=gs_rng,
    )

    psf = mimsim.psfws.make_psfws_psf(
        obsdata=obsdata,
        gs_rng=gs_rng,
        psf_config=psf_config,
    )

    cosmics = mimsim.cosmic_rays.CosmicRays(
        cosmic_ray_rate=cosmic_ray_rate,
        exptime=obsdata['exptime'],
        gs_rng=gs_rng,
    )

    image, sky_image, truth = mimsim.runner.run_sim(
        rng=rng,
        cat=cat,
        obsdata=obsdata,
        artist=artist,
        psf=psf,
        wcs=wcs,
        sky_model=sky_model,
        sensor=sensor,
        dm_detector=dm_detector,
        cosmics=cosmics,
        sky_gradient=gradient,
        vignetting=vignetter,
        fringing=fringer,
        indices=indices,
        apply_pixel_areas=False,  # for speed
        selector=selector,
    )
    assert image.array.shape == sky_image.array.shape
    if limit:
        assert truth.size == indices.size
    else:
        assert truth.size == cat.getNObjects()

    if select:
        assert np.sum(~truth['skipped']) == select_num


@pytest.mark.parametrize(
    'send_type',
    ['default_rng', 'RandomState', 'BaseDeviate']
)
def test_get_rngs(send_type):
    if send_type == 'default_rng':
        rng = np.random.default_rng()
    elif send_type == 'RandomState':
        rng = np.random.RandomState()
    else:
        rng = galsim.BaseDeviate()

    np_rng, gs_rng = mimsim.runner.get_rngs(rng)
    assert isinstance(np_rng, np.random.Generator)
    assert isinstance(gs_rng, galsim.BaseDeviate)

    with pytest.raises(ValueError):
        mimsim.runner.get_rngs(3)


if __name__ == '__main__':
    test_runner_smoke()
