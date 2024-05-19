import galsim
import imsim
import mimsim
import pytest


STAR_EXPECTED = {
    10: 46,
    99: 46,
    101: 46,
    1000000.0: 68,
    100000000.0: 1332,
}


@pytest.mark.parametrize('obj_type', ['star', 'gal'])
def test_stamp_size(obj_type):
    if obj_type == 'star':
        obj = galsim.DeltaFunction()
    else:
        obj = galsim.Exponential(half_light_radius=0.5)

    obj = obj * mimsim.seds.get_trivial_sed()

    obsdata = mimsim.simtools.load_example_obsdata()

    sky_model = imsim.SkyModel(
        exptime=obsdata['exptime'],
        mjd=obsdata['mjd'],
        bandpass=obsdata['bandpass'],
    )

    sky_pos = obsdata['boresight']
    sky_level_density = sky_model.get_sky_level(sky_pos)
    sky_level = 0.2**2 * sky_level_density

    for flux in [10, 99, 101, 1.e6, 1.e8]:
        stamp_size = mimsim.stamps.get_stamp_size(
            obj=obj, flux=flux, noise_var=sky_level, obsdata=obsdata,
        )
        print('obj_type:', obj_type, 'flux:', flux, 'stamp_size:', stamp_size)
        if obj_type == 'gal':
            assert stamp_size == mimsim.defaults.MIN_STAMP_SIZE
        else:
            assert stamp_size == STAR_EXPECTED[flux]


@pytest.mark.parametrize('flux', [100, 1.e7])
def test_initial_draw_method(flux):
    from mimsim.defaults import FFT_FLUX_THRESH, FFT_SB_THRESH
    draw_method = mimsim.stamps.get_initial_draw_method(flux)

    if flux < FFT_FLUX_THRESH or flux < FFT_SB_THRESH:
        assert draw_method == 'phot'
    else:
        assert draw_method == 'fft'


if __name__ == '__main__':
    test_stamp_size('gal')
    test_stamp_size('star')
