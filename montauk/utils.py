def extract_boresight(obsdata):
    """
    Get the boresight from the opsim data as a galsim.CelestialCoord
    """
    import galsim
    return galsim.CelestialCoord(
        obsdata['rightascension'] * galsim.degrees,
        obsdata['declination'] * galsim.degrees,
    )


def get_latitude():
    """
    Get the latitude of the site as a galsim.Angle
    """
    import galsim
    from lsst.obs.lsst.translators.lsst import SIMONYI_LOCATION as RUBIN_LOC

    return RUBIN_LOC.lat.deg * galsim.degrees


def sigma_clip(arrin, niter=4, nsig=4):
    """
    Calculate the mean, sigma, error of an array with sigma clipping.

    parameters
    ----------
    arr: array or sequence
        A numpy array or sequence
    niter: int, optional
        number of iterations, defaults to 4
    nsig: float, optional
        number of sigma, defaults to 4

    returns
    -------
    mean, stdev, err
    """
    import numpy as np

    arr = np.array(arrin, ndmin=1, copy=False)

    if len(arr.shape) > 1:  # pragma: no cover
        raise ValueError(
            'only 1-dimensional arrays suppored, got {arr.shape}'
        )

    indices = np.arange(arr.size)
    nold = arr.size

    mn, sig, err = _get_sigma_clip_stats(arr)

    for i in range(1, niter + 1):

        w, = np.where((np.abs(arr[indices] - mn)) < nsig * sig)

        if w.size == 0:
            # everything clipped, nothing to do but report latest
            # statistics
            break  # pragma: no cover

        if w.size == nold:
            break

        indices = indices[w]
        nold = w.size

        mn, sig, err = _get_sigma_clip_stats(arr[indices])

    return mn, sig, err


def _get_sigma_clip_stats(arr):
    import numpy as np

    mn = arr.mean()
    sig = arr.std()
    err = sig / np.sqrt(arr.size)
    return mn, sig, err
