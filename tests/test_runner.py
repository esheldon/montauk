import galsim
import imsim
import mimsim
import numpy as np
import pytest


def test_runner_smoke():
    mimsim.logging.setup_logging('info')

    seed = 919
    rng = np.random.default_rng(seed)
    gs_rng = galsim.BaseDeviate(rng.integers(0, 2**60))

    band = 'i'
    detnum = 35

    # default is 800, use 100 for speed
    psf_config = {'screen_size': 100}
    cosmic_ray_rate = mimsim.defaults.DEFAULT_COSMIC_RAY_RATE

    dm_detector = mimsim.camera.make_dm_detector(detnum)
    obsdata = mimsim.simtools.load_example_obsdata(band=band)

    wcs, icrf_to_field = mimsim.wcs.make_batoid_wcs(
        obsdata=obsdata, dm_detector=dm_detector,
    )

    cat = mimsim.simtools.load_example_instcat(
        rng=rng, band=band, detnum=detnum,
    )

    sky_model = imsim.SkyModel(
        exptime=obsdata['vistime'],
        mjd=obsdata['mjd'],
        bandpass=obsdata['bandpass'],
    )

    gradient = mimsim.sky.FixedSkyGradient(sky_model)
    vignetter = mimsim.vignetting.Vignetter(dm_detector)

    tree_rings = mimsim.tree_rings.make_tree_rings([detnum])
    sensor = mimsim.sensor.make_sensor(
        dm_detector=dm_detector,
        tree_rings=tree_rings,
        gs_rng=gs_rng,
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
        exptime=obsdata['vistime'],
        band=obsdata['band'],
        dcr=dcr,
        optics=optics,
    )

    diffraction_fft = imsim.stamp.DiffractionFFT(
        exptime=obsdata['vistime'],
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
        exptime=obsdata['vistime'],
        gs_rng=gs_rng,
    )

    mimsim.runner.run_sim(
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
        apply_pixel_areas=False,  # for speed
    )


@pytest.mark.parametrize(
    'options',
    [
        {'band': 'i',
         'limit': False,
         'select': False,
         'cosmics_and_diff_flux': True},
        {'band': 'i',
         'limit': False,
         'select': False,
         'cosmics_and_diff_flux': False},

        {'band': 'i',
         'limit': False,
         'select': True,
         'cosmics_and_diff_flux': True},

        {'band': 'i',
         'limit': True,
         'select': False,
         'cosmics_and_diff_flux': False},
        {'band': 'i',
         'limit': False,
         'select': True,
         'cosmics_and_diff_flux': False},
        {'band': 'Y',
         'limit': False,
         'select': False,
         'cosmics_and_diff_flux': False}
    ]
)
def test_runner(options):
    mimsim.logging.setup_logging('info')

    seed = 919
    rng = np.random.default_rng(seed)
    gs_rng = galsim.BaseDeviate(rng.integers(0, 2**60))

    band = options['band']
    limit = options['limit']
    select = options['select']

    # 88 is E2V, which we want for fringing
    detnum = 88

    # default is 800, use 100 for speed
    psf_config = {'screen_size': 100}
    cosmic_ray_rate = mimsim.defaults.DEFAULT_COSMIC_RAY_RATE

    dm_detector = mimsim.camera.make_dm_detector(detnum)
    obsdata = mimsim.simtools.load_example_obsdata(band=band)

    wcs, icrf_to_field = mimsim.wcs.make_batoid_wcs(
        obsdata=obsdata, dm_detector=dm_detector,
    )

    cat = mimsim.simtools.load_example_instcat(
        rng=rng, band=band, detnum=detnum,
    )
    # make all visible but not too brig
    if not options['cosmics_and_diff_flux']:
        cat.magnorm[:] = 19
    elif select:
        # this will be skipped
        cat.magnorm[-1] = 5

    if limit:
        # indices = np.arange(2)
        indices = np.array([1, 3])
        nobj = indices.size
    else:
        nobj = cat.getNObjects()
        indices = None

    def example_selector(d, i):
        return d.magnorm[i] > magmax

    def default_selector(d, i):
        return True

    if select:
        magmax = 18.5
        selector = example_selector
        select_num = np.sum(cat.magnorm > magmax)
    else:
        selector = default_selector
        select_num = nobj

    sky_model = imsim.SkyModel(
        exptime=obsdata['vistime'],
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

    tree_rings = mimsim.tree_rings.make_tree_rings([detnum])
    sensor = mimsim.sensor.make_sensor(
        dm_detector=dm_detector,
        tree_rings=tree_rings,
        gs_rng=gs_rng,
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
        exptime=obsdata['vistime'],
        band=obsdata['band'],
        dcr=dcr,
        optics=optics,
    )

    diffraction_fft = imsim.stamp.DiffractionFFT(
        exptime=obsdata['vistime'],
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

    if options['cosmics_and_diff_flux']:
        cosmics = mimsim.cosmic_rays.CosmicRays(
            cosmic_ray_rate=cosmic_ray_rate,
            exptime=obsdata['vistime'],
            gs_rng=gs_rng,
        )
    else:
        cosmics = None

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

    # cosmics not in truth
    if not options['cosmics_and_diff_flux']:
        exposure = _image_to_exposure(
            image=image, sky_image=sky_image, band=band,
        )
        sources = _detect(exposure=exposure)

        flux = sources['base_PsfFlux_instFlux']
        flux_err = sources['base_PsfFlux_instFluxErr']
        w, = np.where(flux / flux_err > 10)

        assert w.size == truth.size


def _doplot(image, sky_image, truth, sources):
    import matplotlib.pyplot as mplt

    fig, ax = mplt.subplots()
    ims = np.log10((image.array - sky_image.array).clip(min=0.1))
    ax.imshow(ims)
    flux = sources['base_PsfFlux_instFlux']
    flux_err = sources['base_PsfFlux_instFluxErr']
    w, = np.where(flux / flux_err > 10)
    ax.scatter(
        sources['base_SdssCentroid_x'][w],
        sources['base_SdssCentroid_y'][w],
        c='red',
    )
    mplt.show()


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


# Maybe move this to a module

def _detect(exposure, thresh=5):
    import lsst.afw.table as afw_table
    from lsst.meas.base import (
        SingleFrameMeasurementConfig,
        SingleFrameMeasurementTask,
    )
    from lsst.meas.algorithms import (
        SourceDetectionTask, SourceDetectionConfig,
    )
    from lsst.meas.deblender import SourceDeblendTask, SourceDeblendConfig
    # from metadetect.lsst.util import get_stats_mask

    schema = afw_table.SourceTable.makeMinimalSchema()

    # Setup algorithms to run
    meas_config = SingleFrameMeasurementConfig()
    meas_config.plugins.names = [
        'base_PixelFlags',
        'base_SdssCentroid',
        'base_PsfFlux',
        'base_SkyCoord',
        'base_SdssShape',
    ]

    # set these slots to none because we aren't running these algorithms
    meas_config.slots.apFlux = None
    meas_config.slots.gaussianFlux = None
    meas_config.slots.calibFlux = None
    meas_config.slots.modelFlux = None

    # fix odd issue where it things things are near the edge
    meas_config.plugins['base_SdssCentroid'].binmax = 1

    # sub-pixel offsets in the psf rendering
    # meas_config.plugins['ext_shapeHSM_HsmPsfMoments'].useSourceCentroidOffset = True  # noqa

    meas_task = SingleFrameMeasurementTask(
        config=meas_config,
        schema=schema,
    )

    detection_config = SourceDetectionConfig()
    detection_config.thresholdValue = thresh

    # these will be ignored when finding the image standard deviation
    # detection_config.statsMask = get_stats_mask(exposure)

    detection_task = SourceDetectionTask(config=detection_config)

    # these tasks must use the same schema and all be constructed before
    # any other tasks, and must use the same schema because the schema is
    # modified in place by tasks, and the constructor does a check that
    # fails if we do this afterward

    deblend_task = SourceDeblendTask(
        config=SourceDeblendConfig(),
        schema=schema,
    )

    table = afw_table.SourceTable.make(schema)

    result = detection_task.run(table, exposure)

    if result is not None:
        sources = result.sources
        deblend_task.run(exposure, sources)

        for source in sources:

            if source.get('deblend_nChild') != 0:
                continue

            meas_task.callMeasure(source, exposure)

    else:
        sources = []

    return sources


def _image_to_exposure(image, sky_image, band):
    import lsst.afw.image as afw_image

    ny, nx = image.array.shape
    masked_image = afw_image.MaskedImageF(nx, ny)
    masked_image.image.array[:, :] = image.array
    masked_image.image.array[:, :] -= sky_image.array
    masked_image.variance.array[:, :] = sky_image.array

    exposure = afw_image.ExposureF(masked_image)

    filter_label = afw_image.FilterLabel(band=band, physical=band)
    dm_wcs = _gs_wcs_to_dm_wcs(image.wcs, image.bounds)
    dm_psf = _make_fixed_psf(fwhm=0.8)

    exposure.setFilter(filter_label)
    exposure.setPsf(dm_psf)
    exposure.setWcs(dm_wcs)
    return exposure


def _gs_wcs_to_dm_wcs(gs_wcs, bounds):
    header = {}
    gs_wcs.writeToFitsHeader(header, bounds)
    return _header_to_wcs(header)


def _header_to_wcs(hdr):
    """
    convert a header to a WCS

    Note this will not capture tree rings
    """
    from lsst.daf.base import PropertyList
    from lsst.afw.geom import makeSkyWcs

    prop = PropertyList()
    for key in hdr:
        if key[:3] == 'GS_':
            continue
        prop.set(key, hdr[key])

    return makeSkyWcs(prop)


def _make_fixed_psf(fwhm, rng=None):
    """
    make a KernelPsf(FixedKernel()) for a gaussian with the input fwhm
    """
    import galsim
    from lsst.meas.algorithms import KernelPsf
    from lsst.afw.math import FixedKernel
    import lsst.afw.image as afw_image

    g = galsim.Gaussian(fwhm=fwhm)
    psf_image = g.drawImage(scale=0.2, nx=25, ny=25).array

    if rng is not None:
        noise = psf_image.max() / 1000
        psf_image += rng.normal(scale=noise, size=psf_image.shape)

    psf_image = psf_image.astype(float)

    return KernelPsf(
        FixedKernel(afw_image.ImageD(psf_image))
    )


if __name__ == '__main__':

    test_runner_smoke()
    # options = {
    #     'band': 'i',
    #     'limit': False,
    #     'select': True,
    #     'cosmics_and_diff_flux': True,
    # }
    # test_runner(options)
