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
