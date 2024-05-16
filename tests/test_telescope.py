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


if __name__ == '__main__':
    test_make_detector_telescope()
