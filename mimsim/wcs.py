import imsim
from imsim.utils import ignore_erfa_warnings


def make_batoid_wcs(
    obsdata,
    dm_detector,
    order=3,
    temperature=280,
    pressure=72.7,
    H2O_pressure=1.0,

):
    """
    Parameters
    ----------
    obsdata: dict
        Information about this observation
    dm_detector: lsst.afw.cameraGeom.Detector
        Data management detector object.  Use make_dm_detector(detnum)
    order: int, optional
        Order of the SIP polynomial [default: 3]
    temperature : float, optional
        Ambient temperature in Kelvin [default: 280 K]
    pressure : float, optional
        Ambient pressure in kPa [default: based on LSST heigh of 2715 meters]
    H2O_pressure : float, optional
        Water vapor pressure in kPa [default: 1.0 kPa]

    Returns
    --------
    wcs, icrf_to_field
        Both are galsim.FittedSIPWCS
    """
    import astropy.time
    import imsim
    from .telescope import make_detector_telescope

    obstime = astropy.time.Time(
        obsdata['mjd'], format='mjd', scale='tai',
    )
    wavelength = obsdata['bandpass'].effective_wavelength

    detector_telescope = make_detector_telescope(
        band=obsdata['band'],
        rot_tel_pos=obsdata['rotTelPos'],
        dm_detector=dm_detector,
    )

    factory = imsim.batoid_wcs.BatoidWCSFactory(
        boresight=obsdata['boresight'],
        obstime=obstime,
        telescope=detector_telescope,
        wavelength=wavelength,
        camera='LsstCam',
        temperature=temperature,
        pressure=pressure,
        H2O_pressure=H2O_pressure,
    )
    return (
        factory.getWCS(det=dm_detector, order=order),
        factory.get_icrf_to_field(det=dm_detector, order=order),
    )


def get_pixel_scale(wcs, bbox):
    """
    get the pixel scale at the center of the image

    Parameters
    ----------
    wcs: galsim.FitsWCS
        The galsim wcs
    bbox: lsst.geom.Box2I
        The image bounding box
    """
    import numpy as np
    import galsim

    x = bbox.getWidth() / 2
    y = bbox.getHeight() / 2
    return np.sqrt(wcs.pixelArea(galsim.PositionD(x, y)))


class BatoidDCRWCSFactory(imsim.batoid_wcs.BatoidWCSFactory):
    """
    Factory for constructing WCS's, including a DCR shift

    Parameters
    ----------
    boresight : galsim.CelestialCoord
        The ICRF coordinate of light that reaches the boresight.  Note that
        this is distinct from the spherical coordinates of the boresight with
        respect to the ICRF axes.
    obstime : astropy.time.Time
        Mean time of observation.
    telescope : batoid.Optic
        Telescope instance. Should include any camera rotation.
    wavelength : float
        Nanometers
    camera : lsst.afw.cameraGeom.Camera
    temperature : float
        Ambient temperature in Kelvin
    pressure : float
        Ambient pressure in kPa
    H2O_pressure : float
        Water vapor pressure in kPa
    dcr: a mimsim.dcr.DCRMaker
        a mimsim.dcr.DCRMaker that can generate dcr at a specified location
    """
    @ignore_erfa_warnings
    def __init__(
        self,
        boresight,
        obstime,
        telescope,
        wavelength,
        camera,
        temperature,
        pressure,
        H2O_pressure,
        dcr,
    ):
        super().__init__(
            boresight=boresight,
            obstime=obstime,
            telescope=telescope,
            wavelength=wavelength,
            camera=camera,
            temperature=temperature,
            pressure=pressure,
            H2O_pressure=H2O_pressure,
        )
        self.dcr = dcr

    def getWCS(self, det, order=3):
        """
        Parameters
        ----------
        det : lsst.afw.cameraGeom.Detector
            Detector of interest.
        order : int
            SIP order for fit.

        Returns
        -------
        wcs : galsim.fitswcs.GSFitsWCS
            WCS transformation between ICRF <-> pixels.
        """
        import galsim

        thxs, thys = self.get_field_samples(det)
        z_offset = imsim.batoid_wcs.det_z_offset(det)

        # trace both directions (field -> ICRF and field -> pixel)
        # then fit TanSIP to ICRF -> pixel.
        fpxs, fpys = self._field_to_focal(thxs, thys, z_offset=z_offset)
        xs, ys = imsim.batoid_wcs.focal_to_pixel(fpxs, fpys, det)
        rob, dob = self._field_to_observed(thxs, thys)
        rc, dc = self._observed_to_ICRF(rob, dob)

        tmp_wcs = galsim.FittedSIPWCS(xs, ys, rc, dc, order=order)

        # now apply DCR shifts. The order of operations is not
        # right here, but I think it is equivalent to the order
        # of photon ops used for imsim runs:
        #
        # TimeSampler: modifies time
        # PupilAnnulusSampler: modifies pupil_u, pupil_v but does not
        #    modify x, y.  This is used for the optics portion, so I
        #    think it could come after dcr.
        # PhotonDCR: modifies x, y but ignores pupil_u, pupil_v
        # RubinDiffractionOptics: uses pupil_u, pupil_v to calculate
        #    a shift in x, y
        # FocusDepth: uses x, y to compute shift in x, y
        #
        # These shifts are all just added, so our approach should get the DCR
        # part right, while missing the effects of other ops

        photon_array = galsim.PhotonArray(N=1, wavelength=self.wavelength)

        for i in range(xs.size):

            image_pos = galsim.PositionD(xs[i], ys[i])

            photon_array.x[0] = image_pos.x
            photon_array.y[0] = image_pos.y

            sky_pos = galsim.CelestialCoord(rc[i], dc[i])
            this_dcr = self.dcr(sky_pos)

            local_wcs = tmp_wcs.local(image_pos)
            this_dcr.applyTo(photon_array, local_wcs=local_wcs)

            xs[i] = photon_array.x[0]
            ys[i] = photon_array.y[0]

        # refit with shifted positions
        return galsim.FittedSIPWCS(xs, ys, rc, dc, order=order)
