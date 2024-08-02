def make_fixed_psf(type, options):
    """
    make a fixed psf as a GSObject

    Parameters
    ----------
    type: str
        e.g. 'gauss' 'moffat'
    options: dict
        Options for the construction of the PSF

    Returns
    -------
    GSObject
    """
    import galsim

    if type == 'gauss':
        psf = galsim.Gaussian(**options)
    elif type == 'moffat':
        psf = galsim.Moffat(**options)
    else:
        raise ValueError(f'unsupported psf type: {type}')

    return psf
