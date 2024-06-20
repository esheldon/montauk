import galsim
import montauk


def test_trivial_sed():
    tsed = montauk.seds.get_trivial_sed()
    expected = galsim.SED(
        galsim.LookupTable([100, 2000], [1, 1], interpolant='linear'),
        wave_type='nm', flux_type='fphotons'
    )

    assert tsed == expected
