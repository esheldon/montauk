def run_sim(
    rng,
    cat,
    obsdata,
    artist,
    psf,
    wcs,
    sky_model,
    sensor,
    dm_detector,
    cosmics=None,
    sky_gradient=None,
    vignetting=None,
    fringing=None,
    apply_pixel_areas=True,
    calc_xy_indices=None,
    skip_bright=False,
    selector=lambda d, i: True,
):
    """
    An example simulation runner

    Parameters
    ----------
    rng: random number generator
        One of
            - np.random.Generator (from np.random.default_rng)
            - np.random.RandomState
            - galsim.BaseDeviate
    cat: imsim.InstCatalog
        The input catalog
    obsdata: dict
        Information about this observation
    artist: artist object, e.g. Artist
        An artist object with draw_phot and draw_fft methods
    psf: the psf creator
        Must be a galsim object or have the getPSF(image_pos) method
    wcs: galsim.FITSWCS
        The wcs for the image
    sky_model: isim.SkyModel
        The imsim sky model
    sensor: galsim.SiliconSensor, optional
        The sensor to use for modifying pixel areas
    dm_detector: lsst.afw.cameraGeom.Detector
        Data management detector object.  Use make_dm_detector(detnum)
    cosmics: cosmic ray adder, optional
        Must have a add(image) method
    sky_gradient: gradient object, optional
        Must have an apply(image) method
    vignetting: vignetting object, optional
        Must have an apply(image) method
    fringing: fringing object, optional
        Must have an apply(image) method
    apply_pixel_areas: bool, optional
        If set to False, do not apply pixel areas.  This saves a large amount
        of time before drawing begins, good for testing.  Default True
    calc_xy_indices: array, optional
        If sent, objects corresponding to these indices will have their
        realized x,y positions calculated and stored in the truth output as
        realized_x and realized_y
    skip_bright: bool, optional
        Skip objects so bright we would use an FFT rather than photon shooting
    selector: function
        A way to select on the catalog, must be callable with
            selector(cat, iobj)
        and return True if the object is to be kept
    """
    import numpy as np
    import logging
    import galsim
    from esutil.pbar import prange
    import imsim
    from imsim.psf_utils import get_fft_psf_maybe
    from .sky import make_sky_image
    from .psf import eval_psf
    from .stamps import get_stamp_size, get_initial_draw_method
    from .wcs import get_pixel_scale
    from .defaults import FFT_SB_THRESH

    logger = logging.getLogger('imsim-runner')

    np_rng, gs_rng = get_rngs(rng)

    bbox = dm_detector.getBBox()

    image = galsim.ImageF(bbox.width, bbox.height, wcs=wcs)
    pixel_scale = get_pixel_scale(wcs=wcs, bbox=bbox)

    sky_image = make_sky_image(
        sky_model=sky_model,
        wcs=wcs,
        nx=bbox.width,
        ny=bbox.height,
        logger=logger,
        gradient=sky_gradient,
        vignetting=vignetting,
        fringing=fringing,
        sensor=None if not apply_pixel_areas else sensor,
    )
    med_noise_var = np.median(sky_image.array)

    nobj = cat.getNObjects()

    nskipped_low_flux = 0
    nskipped_select = 0
    nskipped_bounds = 0
    nskipped_bright = 0

    truth = make_truth(nobj, with_realized_pos=calc_xy_indices is not None)

    for iobj in prange(nobj):

        obj_coord = cat.world_pos[iobj]
        image_pos = cat.image_pos[iobj]

        truth['x'][iobj] = image_pos.x
        truth['y'][iobj] = image_pos.y
        truth['ra'][iobj] = obj_coord.ra.deg
        truth['dec'][iobj] = obj_coord.dec.deg

        if not selector(cat, iobj):
            nskipped_select += 1
            continue

        obj = cat.getObj(index=iobj, rng=gs_rng, exptime=obsdata['exptime'])

        flux = obj.calculateFlux(obsdata['bandpass'])
        truth['nominal_flux'][iobj] = flux

        if flux <= 0:  # pragma: no cover
            nskipped_low_flux += 1
            continue

        imsim.stamp.LSST_SiliconBuilder._fix_seds(
            prof=obj, bandpass=obsdata['bandpass'], logger=logger,
        )

        psf_at_pos = eval_psf(psf=psf, image_pos=image_pos)

        draw_method = get_initial_draw_method(flux)

        if draw_method != 'phot' and skip_bright:
            nskipped_bright += 1
            continue

        if draw_method != 'phot':
            psf_at_pos, draw_method, fft_flux = get_fft_psf_maybe(
                obj=obj,
                nominal_flux=flux,
                psf=psf_at_pos,
                bandpass=obsdata['bandpass'],
                wcs=wcs,
                fft_sb_thresh=FFT_SB_THRESH,
                pixel_scale=pixel_scale,
                dm_detector=dm_detector,
                vignetting=vignetting.vignetting,
                sky_pos=obj_coord,
            )

        if draw_method != 'phot' and skip_bright:
            nskipped_bright += 1
            continue

        stamp_size = get_stamp_size(
            obj=obj, flux=flux, noise_var=med_noise_var, obsdata=obsdata,
            pixel_scale=pixel_scale,
        )

        logger.debug('draw_method: %s stamp_size: %s flux: %s',
                     draw_method, stamp_size, flux)
        local_wcs = wcs.local(image_pos=image_pos)

        if draw_method == 'phot':
            stamp = artist.phot_draw(
                obj=obj,
                obj_coord=obj_coord,
                image_pos=image_pos,
                flux=flux,
                stamp_size=stamp_size,
                local_wcs=local_wcs,
                psf=psf_at_pos,
            )
        else:
            if fft_flux != flux:
                obj = obj.withFlux(fft_flux, obsdata['bandpass'])

            stamp = artist.fft_draw(
                obj=obj,
                image_pos=image_pos,
                stamp_size=stamp_size,
                local_wcs=local_wcs,
                psf=psf_at_pos,
            )

        bounds = stamp.bounds & image.bounds
        if not bounds.isDefined():  # pragma: no cover
            nskipped_bounds += 1
            continue

        # if wcs_fitter is not None and iobj in wcs_fitter:
        if calc_xy_indices is not None and iobj in calc_xy_indices:
            posdata = artist.get_pos(
                obj=obj, obj_coord=obj_coord, image_pos=image_pos,
                stamp_size=stamp_size, local_wcs=local_wcs,
                psf=psf_at_pos,
            )
            truth['realized_x'][iobj] = posdata['x']
            truth['realized_y'][iobj] = posdata['y']

        image[bounds] += stamp[bounds]
        truth['realized_flux'][iobj] = stamp.added_flux
        truth['skipped'][iobj] = False

    # variance contains poisson level from signal
    # this will itself be noisy. DM does this
    variance = sky_image.array + image.array

    # add noise from sky, noise from objects already in
    image.array[:, :] += np_rng.poisson(lam=sky_image.array)

    # should go in after poisson noise
    if cosmics is not None:
        logger.info('adding cosmic rays')
        cosmics.add(image)

    nskipped = (
        nskipped_select + nskipped_low_flux + nskipped_bounds + nskipped_bright
    )
    logger.info(f'skipped {nskipped}/{nobj}')
    logger.info(f'skipped {nskipped_low_flux}/{nobj} low flux')
    logger.info(f'skipped {nskipped_select}/{nobj} selector')
    logger.info(f'skipped {nskipped_bounds}/{nobj} bounds')
    logger.info(f'skipped {nskipped_bright}/{nobj} bright')

    return image, variance, sky_image, truth


def make_truth(nobj, with_realized_pos=False):
    """
    Make the truth for run_sim.  The catalog will have fields
        ('skipped', bool),
        ('ra', 'f8'),
        ('dec', 'f8'),
        ('x', 'f4'),
        ('y', 'f4'),
        ('nominal_flux', 'f4'),
        ('realized_flux', 'f4'),

    Parameters
    ----------
    nobj: int
        The number of rows

    Returns
    -------
    array
    """
    import numpy as np

    dtype = [
        ('skipped', bool),
        ('ra', 'f8'),
        ('dec', 'f8'),
        ('x', 'f4'),
        ('y', 'f4'),
        ('nominal_flux', 'f4'),
        ('realized_flux', 'f4'),
    ]
    if with_realized_pos:
        dtype += [
            ('realized_x', 'f4'),
            ('realized_y', 'f4'),
        ]
    st = np.zeros(nobj, dtype=dtype)
    st['skipped'] = True
    st['realized_flux'] = np.nan

    if with_realized_pos:
        st['realized_x'] = np.nan
        st['realized_y'] = np.nan

    return st


def get_rngs(rng):
    """
    extract numpy and galsim random number generators

    Parameters
    ----------
    rng:
        One of
            - np.random.Generator (from np.random.default_rng)
            - np.random.RandomState
            - galsim.BaseDeviate

    Returns
    -------
    np_rng, gs_rng:
        a np.random.default_rng and a galsim.BaseDeviate
    """
    import numpy as np
    import galsim

    if isinstance(rng, galsim.BaseDeviate):
        np_rng = np.random.default_rng(rng.raw())
        gs_rng = rng
    elif isinstance(rng, np.random.Generator):
        gs_rng = galsim.BaseDeviate(rng.integers(0, 2**60))
        np_rng = rng
    elif isinstance(rng, np.random.RandomState):
        np_rng = np.random.default_rng(rng.randint(0, 2**60))
        gs_rng = galsim.BaseDeviate(rng.randint(0, 2**60))
    else:
        raise ValueError(
            f'got rng {type(rng)}, expected one of np.random.default_rng, '
            'np.random.RandomState, or galsim.BaseDeviate'
        )

    return np_rng, gs_rng
