

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


def fit_wcs(x, y, ra, dec, units, itrain, ireserve, order=3):
    """
    Fit a SIPWCS to the input data

    Parameters
    ----------
    x: array
        Array of x positions for training
    y: array
        Array of y positions for training
    ra: array
        Array of ra for training
    dec: array
        Array of dec for training
    units: galsim.AngleUnit
        The units of the ra/dec input
    itrain: array
        Array of indices to use for training
    ireserve: array
        Array of indices to reserve for validation
    order: int, optional
        Order of fit, default 3

    Returns
    -------
    a tuple of galsim.GSFitsWCS, stats
        stats is a dict with entries
            ntrain: number used for training
            nreserve: number used for validation
            xoff_mean: mean of x offset between predicted and validation
            xoff_err: uncertainty in mean of x offset between predicted
                and validation
            yoff_mean: mean of y offset between predicted and validation
            yoff_err: uncertainty in mean of y offset between predicted
                and validation
    """
    import numpy as np
    import galsim
    from .utils import sigma_clip

    if units == galsim.degrees:
        ra = np.radians(ra)
        dec = np.radians(dec)

    xtrain = x[itrain]
    ytrain = y[itrain]
    ratrain = ra[itrain]
    dectrain = dec[itrain]

    wt, = np.where(np.isfinite(xtrain))
    wcs = galsim.FittedSIPWCS(
        xtrain[wt], ytrain[wt], ratrain[wt], dectrain[wt], order=order,
    )

    xreserve = x[ireserve]
    yreserve = y[ireserve]
    rareserve = ra[ireserve]
    decreserve = dec[ireserve]

    wr, = np.where(np.isfinite(xreserve))
    xcheck, ycheck = wcs.radecToxy(
        ra=rareserve[wr], dec=decreserve[wr],
        units=galsim.radians,
    )

    xcoff = xcheck - xreserve[wr]
    ycoff = ycheck - yreserve[wr]

    xcoff_mean, _, xcoff_err = sigma_clip(xcoff)
    ycoff_mean, _, ycoff_err = sigma_clip(ycoff)
    stats = {
        'ntrain': wt.size,
        'nreserve': wr.size,
        'xoff_mean': xcoff_mean,
        'xoff_err': xcoff_err,
        'yoff_mean': ycoff_mean,
        'yoff_err': ycoff_err,
    }
    return wcs, stats


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
