class DCRMaker(object):
    """
    A functor class to create a galsim.PhotonDCR object

    Parameters
    ----------
    bandpass: imsim.bandpass.RubinBandpass
        The bandpass for this image
    hour_angle: galsim.Angle
        The hour angle

    Returns
    --------
    galsim.PhotonDCR
    """

    def __init__(
        self,
        bandpass,
        hour_angle,
    ):
        self.bandpass = bandpass
        self.hour_angle = hour_angle

    def __call__(self, sky_pos):
        import galsim
        from .utils import get_latitude

        return galsim.PhotonDCR(
            base_wavelength=self.bandpass.effective_wavelength,
            latitude=get_latitude(),
            HA=self.hour_angle,
            obj_coord=sky_pos,
        )
