def make_sky_image(
    sky_model,
    wcs,
    nx,
    ny,
    gradient=True,
    vignetting=True,
    fringing=None,
    sensor=None,
    logger=None,
):
    """
    Convenience function to make a sky image.  This is useful to ensure we get
    the ordering of operations right

    Parameters
    ----------
    sky_model: imsim.sky_model.SkyModel
        The imsim model for the sky
    wcs: galsim.GSFitsWCS
        The WCS for the image
    nx: int
        Image size in x
    ny: int
        Image size in y
    gradient: gradient object, optional
        Must have an apply(image) method
    vignetting: vignetting object, optional
        Must have an apply(image) method
    fringing: fringing object, optional
        Must have an apply(image) method
    sensor: galsim.SiliconSensor, optional
        The sensor to use for modifying pixel areas
    logger: python logger, optional
        logger

    Returns
    -------
    galsim.ImageF
    """
    import galsim

    # constant so we can use value at the center
    x = nx / 2
    y = ny / 2
    image_pos = galsim.PositionD(x=x, y=y)
    sky_pos = wcs.toWorld(image_pos)
    sky_level = sky_model.get_sky_level(sky_pos)

    # this can be variable due to small pixel size variations
    # due to non uniform WCS
    sky_image = galsim.ImageF(nx, ny, wcs=wcs)

    # values will vary due to wcs pixel size variations
    wcs.makeSkyImage(sky_image, sky_level)

    if gradient is not None:
        if logger is not None:
            logger.info('applying sky gradient')
        gradient.apply(sky_image)

    # should come after sky gradient
    if vignetting is not None:
        if logger is not None:
            logger.info('applying sky vignetting')
        vignetting.apply(sky_image)

    if fringing is not None:
        if logger is not None:
            logger.info('applying fringing')
        fringing.apply(sky_image)

    if sensor is not None:
        # TODO this uses orig_center of (0, 0), make sure correct
        if logger is not None:
            logger.info('applying pixel areas')
        area = sensor.calculate_pixel_areas(sky_image)
        sky_image *= area

    return sky_image


class FixedSkyGradient(object):
    """
    A class to apply a sky gradient

    Parameters
    ----------
    sky_model: isim.SkyModel
        The imsim sky model
    """
    def __init__(self, sky_model):
        self.sky_model = sky_model

    def apply(self, image):
        """
        Apply a sky gradient to the input image

        Parameters
        ----------
        image: galsim.Image
            The input galsim image
        """
        import numpy as np
        import imsim

        ny, nx = image.array.shape

        wcs = image.wcs
        world_center = wcs.toWorld(image.true_center)

        sky_gradient = imsim.sky_model.SkyGradient(
            sky_model=self.sky_model,
            wcs=image.wcs,
            world_center=world_center,
            image_xsize=nx,
        )
        xarr, yarr = np.meshgrid(range(nx), range(ny))
        image.array[:] *= sky_gradient(xarr, yarr)
