class CosmicRays(object):
    """
    A class to add cosmic rays to the image

    Parameters
    ----------
    cosmic_ray_rate: float
        Rate per ccd
    exptime: float
        Exposure time in seconds
    image: galsim.Image
        The image
    gs_rng: galsim rng
        A galsim random number generator
    """
    def __init__(
        self,
        cosmic_ray_rate,
        exptime,
        gs_rng,
    ):
        import os
        import imsim
        from imsim.meta_data import data_dir

        self.cosmic_ray_rate = cosmic_ray_rate
        self.exptime = exptime
        self.gs_rng = gs_rng

        if cosmic_ray_rate > 0:
            path = os.path.join(data_dir, 'cosmic_rays_itl_2017.fits.gz')

            self.cosmic_rays = imsim.cosmic_rays.CosmicRays(
                ccd_rate=cosmic_ray_rate,
                catalog_file=path,
            )

    def add(self, image):
        """
        Add cosmic rays to a image

        Parameters
        ----------
        image: galsim.Image
            The image
        """

        if self.cosmic_ray_rate > 0:
            self.cosmic_rays.paint(
                image_array=image.array,
                rng=self.gs_rng,
                exptime=self.exptime,
            )
