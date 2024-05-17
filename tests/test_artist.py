import imsim
import galsim
import mimsim


def test_artist():

    seed = 1221
    gs_rng = galsim.BaseDeviate(seed)

    detnum = 91
    dm_detector = mimsim.camera.make_dm_detector(detnum)
    obsdata = mimsim.simtools.load_example_obsdata()

    wcs, icrf_to_field = mimsim.wcs.make_batoid_wcs(
        obsdata=obsdata, dm_detector=dm_detector,
    )

    # image_pos = galsim.PositionD(250.1, 1091.3)
    # sky_pos = wcs.toWorld(image_pos)

    dcr = mimsim.dcr.DCRMaker(
        bandpass=obsdata['bandpass'],
        hour_angle=obsdata['HA'],
    )

    tree_rings = mimsim.tree_rings.make_tree_rings([dm_detector])
    sensor = mimsim.sensor.make_sensor(
        dm_detector=dm_detector,
        gs_rng=gs_rng,
        tree_rings=tree_rings,
    )

    optics = mimsim.optics.OpticsMaker(
        altitude=obsdata['altitude'],
        azimuth=obsdata['azimuth'],
        boresight=obsdata['boresight'],
        rot_tel_pos=obsdata['rotTelPos'],
        band=obsdata['band'],
        dm_detector=dm_detector,
        wcs=wcs,
        icrf_to_field=icrf_to_field,
    )

    photon_ops_maker = mimsim.photon_ops.PhotonOpsMaker(
        exptime=obsdata['exptime'],
        band=obsdata['band'],
        dcr=dcr,
        optics=optics,
    )

    diffraction_fft = imsim.stamp.DiffractionFFT(
        exptime=obsdata['exptime'],
        altitude=obsdata['altitude'],
        azimuth=obsdata['azimuth'],
        rotTelPos=obsdata['rotTelPos'],
    )

    mimsim.artist.Artist(
        bandpass=obsdata['bandpass'],
        sensor=sensor,
        photon_ops_maker=photon_ops_maker,
        diffraction_fft=diffraction_fft,
        gs_rng=gs_rng,
    )


if __name__ == '__main__':
    test_artist()
