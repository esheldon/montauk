import galsim
import montauk
import pytest


def test_make_sensor_smoke():
    seed = 1991
    gs_rng = galsim.BaseDeviate(seed)

    detnum = 15
    dm_detector = montauk.camera.make_dm_detector(detnum)

    _ = montauk.sensor.make_sensor(dm_detector=dm_detector, gs_rng=gs_rng)


@pytest.mark.parametrize('use_tree_rings', [False, True])
def test_make_sensor_tree_rings(use_tree_rings):
    seed = 2112
    gs_rng = galsim.BaseDeviate(seed)

    detnums = [85, 120]
    dm_detectors = [
        montauk.camera.make_dm_detector(detnum)
        for detnum in detnums
    ]
    if use_tree_rings:
        tree_rings = montauk.tree_rings.make_tree_rings(detnums)
        assert len(tree_rings.info) == len(dm_detectors)
        for det in dm_detectors:
            assert det.getName() in tree_rings.info
    else:
        tree_rings = None

    for i in range(len(detnums)):
        dm_detector = dm_detectors[i]

        sensor = montauk.sensor.make_sensor(
            dm_detector=dm_detector, gs_rng=gs_rng, tree_rings=tree_rings,
        )
        if use_tree_rings:
            assert sensor.treering_func is not None
            assert sensor.treering_center is not None


if __name__ == '__main__':
    # test_make_sensor_smoke()
    for ut in [True, False]:
        test_make_sensor_tree_rings(ut)
