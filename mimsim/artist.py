from .defaults import DEFAULT_MAX_FLUX_SIMPLE, DEFAULT_MAXN


class Artist(object):
    """
    A class to draw with photons or ffts.  It takes care of some of the logic
    and optional features involved

    Parameters
    ----------
    bandpass: imsim.bandpass.RubinBandpass
        The bandpass for this image
    sensor: galsim.SiliconSensor
        The sensor
    photon_ops_maker: photon ops maker, e.g. PhotonOpsMakre
        An object that can be called with obj_coord to create photon ops
    diffraction_fft: imsim.stamp.DiffractionFFT
        The diffraction fft object
    gs_rng: galsim rng
        A galsim random number generator
    max_flux_simple: float, optional
        Flux less than this will get a simple SED, no photon ops
        or sensor.  Default 100
    maxN: int, optional
        Max number of photons to generate at once, default 1_000_000
    """
    def __init__(
        self,
        bandpass,
        sensor,
        photon_ops_maker,
        diffraction_fft,
        gs_rng,
        max_flux_simple=DEFAULT_MAX_FLUX_SIMPLE,
        maxN=DEFAULT_MAXN,
    ):
        import logging
        self.bandpass = bandpass
        self.sensor = sensor
        self.photon_ops_maker = photon_ops_maker
        self.diffraction_fft = diffraction_fft
        self.gs_rng = gs_rng
        self.max_flux_simple = max_flux_simple
        self.maxN = maxN
        self.logger = logging.getLogger('Artist')

    def phot_draw(
        self, obj, obj_coord, image_pos, flux, stamp_size, local_wcs, psf,
    ):
        """
        Draw with photons

        obj: galsim object
            The galsim object to draw
        obj_coord: galsim.CelestialCoord
            The ra, dec
        image_pos: galsim.PositionD
            Position on image
        flux: float
            The flux for the object
        stamp_size: int
            The size of the stamp
        local_wcs: galsim local wcs
            From wcs.local
        psf: galsim psf
            The psf

        Returns
        -------
        galsim.Image
        """

        if flux < self.max_flux_simple:
            photon_ops = []
            send_sensor = None
            # this is less important for speed than turning off the photon
            # ops/sensor
            obj = get_faint_version_of_obj(obj=obj, bandpass=self.bandpass)
        else:
            photon_ops = self.photon_ops_maker(obj_coord)
            send_sensor = self.sensor

        return obj.drawImage(
            bandpass=self.bandpass,
            nx=stamp_size,
            ny=stamp_size,
            center=image_pos,
            wcs=local_wcs,
            method='phot',
            rng=self.gs_rng,
            maxN=self.maxN,
            photon_ops=[psf] + photon_ops,
            sensor=send_sensor,
        )

    def fft_draw(self, obj, image_pos, stamp_size, local_wcs, psf):
        """
        Draw with fft

        Parameters
        ----------
        obj: galsim object
            The galsim object to draw
        image_pos: galsim.PositionD
            Position on image
        stamp_size: int
            The size of the stamp
        local_wcs: galsim local wcs
            From wcs.local
        psf: galsim psf
            The fft psf

        Returns
        -------
        galsim.Image
        """
        import galsim

        prof = galsim.Convolve([obj] + [psf])
        try:
            stamp = prof.drawImage(
                bandpass=self.bandpass,
                nx=stamp_size, ny=stamp_size,
                center=image_pos,
                wcs=local_wcs,
                method='fft',
            )
        except galsim.errors.GalSimFFTSizeError as e:
            # I think this shouldn't happen with the updates I made to how the
            # image size is calculated, even for extremely bright things.  So
            # it should be ok to just report what happened, give some extra
            # information to diagonose the problem and raise the error.
            logger = self.logger
            logger.error('Caught error trying to draw using FFT:')
            logger.error('%s', e)
            logger.error('You may need to add a gsparams field with '
                         'maximum_fft_size to')
            logger.error('either the psf or gal field to allow larger FFTs.')
            logger.info('prof = %r', prof)
            logger.info('fft_image = %s', stamp)
            raise
        # Some pixels can end up negative from FFT numerics.  Just set them to
        # 0.
        stamp.array[stamp.array < 0] = 0.

        self.diffraction_fft.apply(stamp, self.bandpass.effective_wavelength)

        stamp.addNoise(galsim.PoissonNoise(rng=self.gs_rng))

        return stamp


def get_faint_version_of_obj(obj, bandpass):
    """
    For faint objects we are currently evaluating it at the effective
    wavelength of the bandpass and assigning a trivial SED

    Parameters
    ----------
    obj: galsim object
        The original chromatic object
    bandpass: imsim.bandpass.RubinBandpass
        The bandpass for the image

    Returns
    -------
    galsim object
    """
    from .seds import get_trivial_sed

    output = obj.evaluateAtWavelength(bandpass.effective_wavelength)
    return output * get_trivial_sed()
