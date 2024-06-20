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


def get_narrow_sed(wavemin, wavemax, wave0, wavesig, npts):
    """
    Make a narrow SED with gaussian profile

    Parameters
    ----------
    wavemin: float
        Lower bound on wavelength
    wavemax: float
        Upper bound on wavelength
    wave0: float
        Center of gaussian profile
    wavesig: float
        Sigma of gaussian profile
    npts: float
        Number of points in sed

    Returns
    -------
    galsim.SED with lookup table
    """
    import galsim
    import numpy as np

    wave = np.linspace(wavemin, wavemax, npts)
    height = np.exp(-0.5 * (wave - wave0)**2 / wavesig**2)

    return galsim.SED(
        galsim.LookupTable(
            wave,
            height,
            interpolant='linear',
        ),
        wave_type='nm',
        flux_type='fphotons'
    )
