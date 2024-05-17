from functools import lru_cache


@lru_cache
def get_trivial_sed():
    """
    Get a "trivial" SED, flat between 100 and 2000

    Returns
    -------
    galsim.SED
    """
    import galsim
    return galsim.SED(
        galsim.LookupTable([100, 2000], [1, 1], interpolant='linear'),
        wave_type='nm', flux_type='fphotons'
    )
