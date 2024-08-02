def make_fixed_psf(type, options):
    """
    TODO add more types
    """
    import galsim
    if type == 'gauss':
        psf = galsim.Gaussian(**options)
    elif type == 'moffat':
        psf = galsim.Moffat(**options)
    else:
        raise ValueError(f'unsupported psf type: {psf_type}')
    return psf
