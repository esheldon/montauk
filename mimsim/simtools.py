"""
tools for quick sims.  These are used in the unit tests
"""

EXAMPLE_HEADER = """rightascension 199.76585516403472
declination 1.9773321936155552
mjd 61578.955910303266
altitude 54.54217760455894
azimuth 27.43507998205512
filter %(filter)s
rotskypos 170.12548239658125
camconfig 1
dist2moon 116.5319374452602
moonalt -25.41110589284
moondec 23.9677727
moonphase 75.76275509748939
moonra 316.82542512143766
nsnap 1
obshistid 488069
rottelpos 13.595319056652864
seed 437745637
seeing 0.7732036296985635
sunalt -13.879135662331112
vistime 30.0
minsource 100
seqnum 0
"""

# from a real instcat, but with mags modified to be brighter
EXAMPLE_DATA = """object 1277244820484 148.86127872008032 -26.12814714857972 21.0 starSED/phoSimMLT/lte027-3.0-0.0a+0.0.BT-Settl.spec.gz 0.0 0.0 0.0 0.0 0 0 point none CCM 0.02817156 3.1
object 1606103731204 145.2544608618256 -24.450305013258443 15.0 starSED/phoSimMLT/lte037-5.5-1.0a+0.4.BT-Settl.spec.gz 0.0 0.0 0.0 0.0 0 0 point none CCM 0.02479101 3.1
object 1277253206020 144.99971675234275 -24.991135983257493 22.0 starSED/phoSimMLT/lte038-5.5-0.5a+0.2.BT-Settl.spec.gz 0.0 0.0 0.0 0.0 0 0 point none CCM 0.02479101 3.1
object 7378523209835 147.95869829097995 -25.155485969898194 21.5 galaxySED/Exp.50E07.0005Z.spec.gz 2.694369 -0.021296607 -0.0052690147 -0.021114068 0 0 sersic2d 0.0676198527 0.0462235659 151.743616 1 CCM 0.5 4.000000000000002 CCM 0.0494707944 3.1
object 4544517269611 148.1382735947749 -25.16234573185714 22.0 galaxySED/Burst.12E09.0005Z.spec.gz 0.5205865 0.0047163125 -0.006862338 0.0041362215 0 0 sersic2d 0.101926744 0.0586509332 82.5908693 1 CCM 0.1 4.000000000000002 CCM 0.0647401028 3.1"""  # noqa


def write_example_instcat_header(fobj, band='i'):
    """
    write an example instcat header to the input file object

    Parameters
    -----------
    fobj: file object
        File object to which we will write the header
    """

    filter = 'ugrizy'.index(band.lower())
    fobj.write(EXAMPLE_HEADER % {'filter': filter})


def write_example_instcat_data(fobj, rng, wcs):
    """
    write example data to the input file object

    Parameters
    -----------
    fobj: file object
        File object to which we will write the header
    """

    radec_gen = CCDRadecGenerator(rng=rng, wcs=wcs)

    lines = EXAMPLE_DATA.split('\n')

    ra, dec = radec_gen(n=len(lines))

    for i, line in enumerate(lines):
        ls = line.split()
        entry = _get_instcat_object_line_as_dict(ls)
        entry['ra'] = ra[i]
        entry['dec'] = dec[i]
        _write_instcat_line(fobj, entry)


def load_example_obsdata(band='i'):
    """
    Load example observation data as one would find in an instcat
    header

    Parameters
    ----------
    band: str
        The band for the observation

    Returns
    -------
    imsim.opsim_data.OpsimDataLoader
    """
    import os
    import tempfile
    from .opsim_data import load_obsdata_from_instcat

    with tempfile.TemporaryDirectory() as dir:
        fname = os.path.join(dir, 'instcat.txt')
        with open(fname, 'w') as fobj:
            write_example_instcat_header(fobj, band=band)

        data = load_obsdata_from_instcat(fname, exptime=30)

    return data


def load_example_instcat(rng, band='i', detnum=88):
    """
    Load example observation data as one would find in an instcat
    header

    Parameters
    ----------
    rng: np.random.default_rng
        The random number generator
    band: str, optional
        The band for the observation.  Default 'i'
    detnum: int, optional
        Id of detector, 1-189.  Default 88, an E2V sensor

    Returns
    -------
    imsim.InstCatalog
    """
    import os
    import tempfile
    from .camera import make_dm_detector
    import imsim
    from .wcs import make_batoid_wcs
    from .opsim_data import load_obsdata_from_instcat

    dm_detector = make_dm_detector(detnum)

    obsdata = load_example_obsdata(band=band)
    wcs, _ = make_batoid_wcs(
        obsdata=obsdata, dm_detector=dm_detector,
    )

    with tempfile.TemporaryDirectory() as dir:
        fname = os.path.join(dir, 'instcat.txt')
        with open(fname, 'w') as fobj:
            write_example_instcat_header(fobj, band=band)
            write_example_instcat_data(fobj=fobj, rng=rng, wcs=wcs)

        obsdata = load_obsdata_from_instcat(fname, exptime=30)

        cat = imsim.instcat.InstCatalog(file_name=fname, wcs=wcs)

    return cat


class CCDRadecGenerator():
    """
    generate positions within the CCD region

    Parameters
    ----------
    rng: np.random.default_rng
        The random number generator
    wcs: galsim.FitsWCS
        The wcs for the image
    """
    def __init__(self, rng, wcs):
        self.rng = rng
        self.wcs = wcs

    def __call__(self, n=None):
        import galsim

        if n is None:
            is_scalar = True
            n = 1
        else:
            is_scalar = False

        # DM detectors are not square
        x = self.rng.uniform(low=1, high=4096, size=n)
        y = self.rng.uniform(low=1, high=4004, size=n)

        ra, dec = self.wcs.xyToradec(
            x=x,
            y=y,
            units=galsim.degrees,
        )

        if is_scalar:
            ra = ra[0]
            dec = dec[0]

        return ra, dec


def _get_instcat_object_line_as_dict(ls):
    """
    Read object entries from an instcat

    Parameters
    ----------
    ls: sequence
        the split line from the instcat
    Returns
    -------
    entry: dict holding data
    """

    return {
        'objid': int(ls[1]),
        'ra': float(ls[2]),
        'dec': float(ls[3]),
        'magnorm': float(ls[4]),
        'sed1': ls[5],
        'sed2': float(ls[6]),
        'gamma1': float(ls[7]),
        'gamma2': float(ls[8]),
        'kappa': float(ls[9]),
        'rest': ' '.join(ls[10:]),
    }


def _write_instcat_line(fobj, entry):
    line = ['object'] + [str(v) for k, v in entry.items()]
    line = ' '.join(line)
    fobj.write(line)
    fobj.write('\n')
