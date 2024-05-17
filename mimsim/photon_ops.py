"""
The order used in imsim configs

TimeSampler: modifies time
PupilAnnulusSampler: modifies pupil_u, pupil_v but does not modify x, y
   This is used for the optics portion, so I think it could come after dcr.
PhotonDCR: modifies x, y but ignores pupil_u, pupil_v
RubinDiffractionOptics: uses pupil_u, pupil_v to calculate a shift in x, y
FocusDepth: uses x, y to compute shift in x, y
"""


class PhotonOpsMaker(object):
    """
    A functor class to create the photon ops list

    Position dependent ops dcr and optics can be sent
    as keywords.  They must be callable as op(sky_pos)

    Parameters
    ----------
    exptime: float
        Exposure time for image. Used for the TimeSampler
    band: str
        The band, e.g. u, g, r, i, z. Used to get the focus depth
    dcr: a DCR maker, optional
        Should be callable with dcr(sky_pos)
    optics: an optics maker, optional
        Should be callable with optics(sky_pos)

    Returns
    -------
    A list of photon ops
    """
    def __init__(
        self,
        exptime,
        band,
        dcr=None,
        optics=None,
    ):
        import galsim
        from .optics import make_focus_depth, make_refraction

        # these all in opsim data
        self.exptime = exptime
        self.band = band
        self.dcr = dcr
        self.optics = optics

        self.time_sampler = galsim.TimeSampler(
            t0=0.0, exptime=self.exptime,
        )
        self.pupil_sampler = galsim.PupilAnnulusSampler(
            R_outer=4.18,
            # M1 inner diameter is 2.558, but we need a bit of slack for
            # off-axis rays
            R_inner=2.55,
        )

        self.focus_depth = make_focus_depth(self.band)

        # was for Y (Josh) need for other bands
        self.refraction = make_refraction()

    def __call__(self, sky_pos):
        """
        Create the photon ops list

        Parameters
        ----------
        sky_pos: galsim.CelestialCoord
            The sky position at which to create photon ops

        Returns
        -------
        A list of photon ops
        """

        ops = [self.time_sampler, self.pupil_sampler]

        if self.dcr is not None:
            ops += [self.dcr(sky_pos)]

        if self.optics is not None:
            ops += [self.optics(sky_pos)]

        ops += [self.focus_depth, self.refraction]
        return ops
