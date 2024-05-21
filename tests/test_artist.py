import imsim
import galsim
import mimsim
import pytest


@pytest.mark.parametrize('flux', [1.e6, 10])
def test_artist_smoke(flux):

    seed = 1221
    gs_rng = galsim.BaseDeviate(seed)

    detnum = 91
    dm_detector = mimsim.camera.make_dm_detector(detnum)
    obsdata = mimsim.simtools.load_example_obsdata()

    wcs, icrf_to_field = mimsim.wcs.make_batoid_wcs(
        obsdata=obsdata, dm_detector=dm_detector,
    )

    image_pos = galsim.PositionD(250.1, 1091.3)
    sky_pos = wcs.toWorld(image_pos)

    dcr = mimsim.dcr.DCRMaker(
        bandpass=obsdata['bandpass'],
        hour_angle=obsdata['HA'],
    )

    tree_rings = mimsim.tree_rings.make_tree_rings([detnum])
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
        exptime=obsdata['vistime'],
        band=obsdata['band'],
        dcr=dcr,
        optics=optics,
    )

    diffraction_fft = imsim.stamp.DiffractionFFT(
        exptime=obsdata['vistime'],
        altitude=obsdata['altitude'],
        azimuth=obsdata['azimuth'],
        rotTelPos=obsdata['rotTelPos'],
    )

    artist = mimsim.artist.Artist(
        bandpass=obsdata['bandpass'],
        sensor=sensor,
        photon_ops_maker=photon_ops_maker,
        diffraction_fft=diffraction_fft,
        gs_rng=gs_rng,
    )

    stamp_size = 32

    obj = galsim.Gaussian(fwhm=0.2) * mimsim.seds.get_trivial_sed()

    local_wcs = wcs.local(image_pos)

    psf = galsim.Gaussian(fwhm=0.8)
    artist.phot_draw(
        obj=obj,
        obj_coord=sky_pos,
        image_pos=image_pos,
        flux=flux,
        stamp_size=stamp_size,
        local_wcs=local_wcs,
        psf=psf,
    )

    artist.fft_draw(
        obj=obj,
        image_pos=image_pos,
        stamp_size=stamp_size,
        local_wcs=local_wcs,
        psf=psf,
    )


if __name__ == '__main__':
    test_artist_smoke()
