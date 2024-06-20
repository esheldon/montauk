import montauk


def test_make_camera_smoke():
    camera = montauk.camera.make_camera()

    detnum = 35
    assert len(camera) == 205

    det = camera[detnum]
    assert det.getName() == 'R10_S22'

    det2 = montauk.camera.make_dm_detector(detnum)
    assert det.getName() == det2.getName()


if __name__ == '__main__':
    test_make_camera_smoke()
