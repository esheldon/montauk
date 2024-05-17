import mimsim
import galsim


def test_make_detector_telescope():

    band = 'r'
    rot_tel_pos = 25 * galsim.degrees
    detnum = 91

    dm_detector = mimsim.camera.make_dm_detector(detnum)

    telescope = mimsim.telescope.make_detector_telescope(
        band=band,
        rot_tel_pos=rot_tel_pos,
        dm_detector=dm_detector,
    )
    # what other tests can we do here?
    assert telescope.name == 'LSST'


def test_pupil_sampler():
    obsdata = mimsim.simtools.load_example_obsdata()

    pa = galsim.PhotonArray(
        N=100, wavelength=obsdata['bandpass'].effective_wavelength,
    )

    time_sampler = galsim.TimeSampler(exptime=30)
    pupil_sampler = mimsim.telescope.make_pupil_sampler()

    time_sampler.applyTo(pa)
    pupil_sampler.applyTo(pa)


if __name__ == '__main__':
    test_make_detector_telescope()
