import mimsim


def test_make_camera_smoke():
    camera = mimsim.camera.make_camera()

    detnum = 35
    assert len(camera) == 205

    det = camera[detnum]
    assert det.getName() == 'R10_S22'

    det2 = mimsim.camera.make_dm_detector(detnum)
    assert det.getName() == det2.getName()


if __name__ == '__main__':
    test_make_camera_smoke()
