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
    limit=None,
    apply_pixel_areas=True,
    selector=lambda d, i: True,
):
    """
    An example simulation runner

    Parameters
    ----------
    rng: np.random.default_rng
        A numpy random number generator
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
    limit: int, optional
        Optionally limit to drawing this many objects.  We need this because
        the imsim.InstCatalog cannot read in a subset of the data, and there is
        no way to extract a subset from the InstCatalog.
    apply_pixel_areas: bool, optional
        If set to False, do not apply pixel areas.  This saves a large amount
        of time before drawing begins, good for testing.  Default True
    selector: function
        A way to select on the catalog, must be callable with selector(cat,
        iobj) and return True if the object is to be kept
    """
    import numpy as np
    import logging
    import galsim
    from tqdm import trange
    import imsim
    from imsim.psf_utils import get_fft_psf_maybe
    from .sky import make_sky_image
    from .psf import eval_psf
    from .stamps import get_stamp_size, get_initial_draw_method
    from .wcs import get_pixel_scale
    from .defaults import FFT_SB_THRESH

    logger = logging.getLogger('imsim-runner')

    gs_rng = galsim.BaseDeviate(rng.integers(0, 2**60))

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
    if limit is not None and limit < nobj:
        logger.info(f'will limit to {limit}/{nobj} objects')
        nobj = limit

    nskipped_low_flux = 0
    nskipped_select = 0
    nskipped_bounds = 0

    truth = make_truth(nobj)

    for iobj in trange(nobj):

        obj_coord = cat.world_pos[iobj]
        truth['ra'][iobj] = obj_coord.ra.deg
        truth['dec'][iobj] = obj_coord.dec.deg

        if not selector(cat, iobj):
            nskipped_select += 1
            continue

        obj = cat.getObj(index=iobj, rng=gs_rng, exptime=obsdata['exptime'])

        flux = obj.calculateFlux(obsdata['bandpass'])
        truth['nominal_flux'][iobj] = flux

        if flux <= 0:
            nskipped_low_flux += 1
            continue

        imsim.stamp.LSST_SiliconBuilder._fix_seds(
            prof=obj, bandpass=obsdata['bandpass'], logger=logger,
        )

        image_pos = cat.image_pos[iobj]
        truth['x'][iobj] = image_pos.x
        truth['y'][iobj] = image_pos.y

        psf_at_pos = eval_psf(psf=psf, image_pos=image_pos)

        draw_method = get_initial_draw_method(flux)
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

        stamp_size = get_stamp_size(
            obj=obj, flux=flux, noise_var=med_noise_var, obsdata=obsdata,
            pixel_scale=pixel_scale,
        )

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
        if not bounds.isDefined():
            nskipped_bounds += 1
            continue

        image[bounds] += stamp[bounds]
        truth['realized_flux'][iobj] = stamp.added_flux
        truth['skipped'][iobj] = False

    image.array[:, :] += rng.poisson(lam=sky_image.array)

    # should go in after poisson noise
    logger.info('adding cosmic rays')
    cosmics.add(image)

    nskipped = nskipped_select + nskipped_low_flux + nskipped_bounds
    logger.info(f'skipped {nskipped}/{nobj}')
    logger.info(f'skipped {nskipped_low_flux}/{nobj} low flux')
    logger.info(f'skipped {nskipped_select}/{nobj} selector')
    logger.info(f'skipped {nskipped_bounds}/{nobj} bounds')

    return image, sky_image, truth


def make_truth(nobj):
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
    st = np.zeros(nobj, dtype=dtype)
    st['skipped'] = True
    st['realized_flux'] = np.nan
    return st
