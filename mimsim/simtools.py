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


def write_example_instcat_header(fobj, band='i'):
    """
    write an example instcat header to the input file object

    Parameters
    -----------
    fobj: file object
        File object to which we will write the header
    """

    filter = 'ugrizY'.index(band)
    fobj.write(EXAMPLE_HEADER % {'filter': filter})


def load_example_obsdata(band='i'):
    """
    Load example observation data as one would find in an instcat
    header

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
