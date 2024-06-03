

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

    detector_telescope = make_detector_telescope(
        band=obsdata['band'],
        rot_tel_pos=obsdata['rotTelPos'],
        dm_detector=dm_detector,
    )

    factory = imsim.batoid_wcs.BatoidWCSFactory(
        boresight=obsdata['boresight'],
        obstime=obstime,
        telescope=detector_telescope,
        wavelength=obsdata['bandpass'].effective_wavelength,
        camera='LsstCam',
        temperature=temperature,
        pressure=pressure,
        H2O_pressure=H2O_pressure,
    )
    wcs = factory.getWCS(det=dm_detector, order=order)
    icrf_to_field = factory.get_icrf_to_field(det=dm_detector, order=order)

    return wcs, icrf_to_field


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


class WCSFitterByIndex:
    """
    A WCS fitter with indices into a catalog to provide
    training and validation sets

    Parameters
    ----------
    indices: array
        Array of indices
    reserve: array
        Array of indices used for validation
    order: int, optional
        Default 3
    """
    def __init__(self, indices, reserve, order=3):
        import numpy as np
        self.indices = indices
        self.reserve = reserve

        self.current = 0
        self.x = np.zeros(self.indices.size)
        self.y = np.zeros(self.indices.size)
        self.rarad = np.zeros(self.indices.size)
        self.decrad = np.zeros(self.indices.size)

        self.current_res = 0
        self.x_res = np.zeros(self.reserve.size)
        self.y_res = np.zeros(self.reserve.size)
        self.rarad_res = np.zeros(self.reserve.size)
        self.decrad_res = np.zeros(self.reserve.size)

        self.order = order

    def __contains__(self, index):
        return index in self.indices or index in self.reserve

    def add_entry(self, index, x, y, coord):
        """
        Add an entry to the catalog

        Parameters
        ----------
        index: int
            An index either for the training or validation set
        x: float
            x position of entry
        y: float
            y position of entry
        coord: galsim.CelestialCoord
            ra, dec of object
        """
        if index in self.indices:
            self.x[self.current] = x
            self.y[self.current] = y
            self.rarad[self.current] = coord.ra.rad
            self.decrad[self.current] = coord.dec.rad
            self.current += 1
        elif index in self.reserve:
            self.x_res[self.current_res] = x
            self.y_res[self.current_res] = y
            self.rarad_res[self.current_res] = coord.ra.rad
            self.decrad_res[self.current_res] = coord.dec.rad
            self.current_res += 1

    def fit(self):
        """
        Fit the WCS.  The wcs is returned and is also set
        as an attribute of the object as .wcs

        Returns
        -------
        galsim.GSFitsWCS

        """
        import galsim

        # some may have been skipped
        if self.current < self.x.size:  # pragma: no cover
            self.x = self.x[:self.current]
            self.y = self.y[:self.current]
            self.radrad = self.rarad[:self.current]
            self.decrad = self.decrad[:self.current]

        if self.current_res < self.x.size:  # pragma: no cover
            self.x_res = self.x_res[:self.current]
            self.y_res = self.y_res[:self.current]
            self.radrad_res = self.rarad_res[:self.current]
            self.decrad_res = self.decrad_res[:self.current]

        self.wcs = galsim.FittedSIPWCS(
            self.x, self.y, self.rarad, self.decrad,
            order=self.order,
        )
        return self.wcs

    def get_stats(self):
        """
        Get stats for the fit using the validation set

        Returns
        -------
        stats: dict
            xoff_mean: The sigma clipped mean of the x offset,
                xpredicted - xmeasured
            xoff_err: The sigma clipped erroro on the mean of the x offset,
                xpredicted - xmeasured
            yoff_mean: The sigma clipped mean of the y offset,
                ypredicted - ymeasured
            yoff_err: The sigma clipped erroro on the mean of the y offset,
                ypredicted - ymeasured
        """
        import galsim
        from .utils import sigma_clip

        if not hasattr(self, 'xcheck'):
            self.xcheck, self.ycheck = self.wcs.radecToxy(
                ra=self.rarad_res, dec=self.decrad_res,
                units=galsim.radians,
            )

        xcoff = self.xcheck - self.x_res
        ycoff = self.ycheck - self.y_res

        xcoff_mean, _, xcoff_err = sigma_clip(xcoff)
        ycoff_mean, _, ycoff_err = sigma_clip(ycoff)
        return {
            'ntrain': self.x.size,
            'nreserve': self.x_res.size,
            'xoff_mean': xcoff_mean,
            'xoff_err': xcoff_err,
            'yoff_mean': ycoff_mean,
            'yoff_err': ycoff_err,
        }


# class BatoidDCRWCSFactory(imsim.batoid_wcs.BatoidWCSFactory):
#     """
#     Factory for constructing WCS's, including a DCR shift
#
#     Parameters
#     ----------
#     boresight : galsim.CelestialCoord
#         The ICRF coordinate of light that reaches the boresight.  Note that
#         this is distinct from the spherical coordinates of the boresight with
#         respect to the ICRF axes.
#     obstime : astropy.time.Time
#         Mean time of observation.
#     telescope : batoid.Optic
#         Telescope instance. Should include any camera rotation.
#     wavelength : float
#         Nanometers
#     camera : lsst.afw.cameraGeom.Camera
#     temperature : float
#         Ambient temperature in Kelvin
#     pressure : float
#         Ambient pressure in kPa
#     H2O_pressure : float
#         Water vapor pressure in kPa
#     dcr: a mimsim.dcr.DCRMaker
#         a mimsim.dcr.DCRMaker that can generate dcr at a specified location
#     """
#     @ignore_erfa_warnings
#     def __init__(
#         self,
#         boresight,
#         obstime,
#         telescope,
#         wavelength,
#         camera,
#         temperature,
#         pressure,
#         H2O_pressure,
#         dcr,
#     ):
#         super().__init__(
#             boresight=boresight,
#             obstime=obstime,
#             telescope=telescope,
#             wavelength=wavelength,
#             camera=camera,
#             temperature=temperature,
#             pressure=pressure,
#             H2O_pressure=H2O_pressure,
#         )
#         self.dcr = dcr
#
#     def get_icrf_to_field(self, det, wavelength, order=3):
#         import galsim
#
#         data = self.getWCS(det, wavelength, order=order, more=True)
#         return galsim.FittedSIPWCS(
#             data['thxs'], data['thys'], data['rc'], data['dc'], order=order,
#         )
#
#     def getWCS(self, det, wavelength, order=3, more=False):
#         """
#         Parameters
#         ----------
#         det : lsst.afw.cameraGeom.Detector
#             Detector of interest.
#         wavelength: float
#             Wavelength of light at which to evaluate WCS
#         order : int
#             SIP order for fit.
#
#         Returns
#         -------
#         wcs : galsim.fitswcs.GSFitsWCS
#             WCS transformation between ICRF <-> pixels.
#         """
#         import galsim
#
#         thxs, thys = self.get_field_samples(det)
#         z_offset = imsim.batoid_wcs.det_z_offset(det)
#
#         # trace both directions (field -> ICRF and field -> pixel)
#         # then fit TanSIP to ICRF -> pixel.
#         fpxs, fpys = self._field_to_focal(thxs, thys, z_offset=z_offset)
#         xs, ys = imsim.batoid_wcs.focal_to_pixel(fpxs, fpys, det)
#         rob, dob = self._field_to_observed(thxs, thys)
#         rc, dc = self._observed_to_ICRF(rob, dob)
#
#         tmp_wcs = galsim.FittedSIPWCS(xs, ys, rc, dc, order=order)
#
#         # now apply DCR shifts. The order of operations is not
#         # right here, but I think it is equivalent to the order
#         # of photon ops used for imsim runs:
#         #
#         # TimeSampler: modifies time
#         # PupilAnnulusSampler: modifies pupil_u, pupil_v but does not
#         #    modify x, y.  This is used for the optics portion, so I
#         #    think it could come after dcr.
#         # PhotonDCR: modifies x, y but ignores pupil_u, pupil_v
#         # RubinDiffractionOptics: uses pupil_u, pupil_v to calculate
#         #    a shift in x, y
#         # FocusDepth: uses x, y to compute shift in x, y
#         #
#         # These shifts are all just added, so our approach should get the DCR
#         # part right, while missing the effects of other ops
#
#         photon_array = galsim.PhotonArray(N=1, wavelength=wavelength)
#
#         for i in range(xs.size):
#
#             image_pos = galsim.PositionD(xs[i], ys[i])
#
#             photon_array.x[0] = image_pos.x
#             photon_array.y[0] = image_pos.y
#
#             sky_pos = galsim.CelestialCoord(
#                 rc[i] * galsim.radians, dc[i] * galsim.radians,
#             )
#             this_dcr = self.dcr(sky_pos)
#
#             local_wcs = tmp_wcs.local(image_pos)
#             this_dcr.applyTo(photon_array, local_wcs=local_wcs)
#
#             print(
#                 photon_array.x[0] - xs[i],
#                 photon_array.y[0] - ys[i],
#             )
#             xs[i] = photon_array.x[0]
#             ys[i] = photon_array.y[0]
#
#         # refit with shifted positions
#         wcs = galsim.FittedSIPWCS(xs, ys, rc, dc, order=order)
#
#         if more:
#             return {
#                 'wcs': wcs,
#                 'thxs': thxs,
#                 'thys': thys,
#                 'xs': xs,
#                 'ys': ys,
#                 'rc': rc,
#                 'dc': dc,
#             }
#         else:
#             return wcs
