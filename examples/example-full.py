import galsim
import imsim
import mimsim
import numpy as np


def main():
    mimsim.logging.setup_logging('info')

    seed = 919
    rng = np.random.default_rng(seed)
    gs_rng = galsim.BaseDeviate(rng.integers(0, 2**60))

    band = 'i'
    detnum = 35

    # default is 800, use 100 for speed
    psf_config = {'screen_size': 100}
    cosmic_ray_rate = mimsim.defaults.DEFAULT_COSMIC_RAY_RATE

    dm_detector = mimsim.camera.make_dm_detector(detnum)

    # load the metadata for some example data
    obsdata = mimsim.simtools.load_example_obsdata(band=band)

    # load the objects from the example data
    cat = mimsim.simtools.load_example_instcat(
        rng=rng, band=band, detnum=detnum,
    )

    wcs, icrf_to_field = mimsim.wcs.make_batoid_wcs(
        obsdata=obsdata, dm_detector=dm_detector,
    )

    sky_model = imsim.SkyModel(
        exptime=obsdata['vistime'],
        mjd=obsdata['mjd'],
        bandpass=obsdata['bandpass'],
    )

    gradient = mimsim.sky.FixedSkyGradient(sky_model)
    vignetter = mimsim.vignetting.Vignetter(dm_detector)

    tree_rings = mimsim.tree_rings.make_tree_rings([detnum])
    sensor = mimsim.sensor.make_sensor(
        dm_detector=dm_detector,
        tree_rings=tree_rings,
        gs_rng=gs_rng,
    )

    dcr = mimsim.dcr.DCRMaker(
        bandpass=obsdata['bandpass'],
        hour_angle=obsdata['HA'],
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

    psf = mimsim.psfws.make_psfws_psf(
        obsdata=obsdata,
        gs_rng=gs_rng,
        psf_config=psf_config,
    )

    cosmics = mimsim.cosmic_rays.CosmicRays(
        cosmic_ray_rate=cosmic_ray_rate,
        exptime=obsdata['vistime'],
        gs_rng=gs_rng,
    )

    mimsim.runner.run_sim(
        rng=rng,
        cat=cat,
        obsdata=obsdata,
        artist=artist,
        psf=psf,
        wcs=wcs,
        sky_model=sky_model,
        sensor=sensor,
        dm_detector=dm_detector,
        cosmics=cosmics,
        sky_gradient=gradient,
        vignetting=vignetter,
        apply_pixel_areas=False,  # for speed
    )


if __name__ == '__main__':
    main()
