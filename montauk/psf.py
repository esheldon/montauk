def eval_psf(psf, image_pos):
    """
    evaluate the psf at the specified location.  If the psf is a GSObject
    then the position is ignored

    Paramaters
    ----------
    psf: psf object
        A GSObject or an object with a getPSF method

    Returns
    -------
    psf: GSObject
    """
    if hasattr(psf, 'getPSF'):
        return psf.getPSF(image_pos)
    else:
        return psf
