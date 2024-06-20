# montauk
A library for working with galsim and imsim

Examples
--------
```python
# Example using the provided "run_sim" convenience function
# you don't have to use run_sim, you can write your own
# also, you can replace any of the objects below, such
# as the psf, sky, vignetting etc. as long as they 
# provide the right interface

import galsim
import imsim
import montauk
import numpy as np

montauk.logging.setup_logging('info')

seed = 919
rng = np.random.default_rng(seed)
gs_rng = galsim.BaseDeviate(rng.integers(0, 2**60))

# we will skip fringing, only applies to Y band
band = 'i'
detnum = 35

# default is 800, use 100 for speed
psf_config = {'screen_size': 100}
cosmic_ray_rate = montauk.defaults.DEFAULT_COSMIC_RAY_RATE

dm_detector = montauk.camera.make_dm_detector(detnum)

# load the metadata for some example data
obsdata = montauk.simtools.load_example_obsdata(band=band)

# load the objects from the example data
cat = montauk.simtools.load_example_instcat(
    rng=rng, band=band, detnum=detnum,
)

# this WCS does not include DCR
wcs, icrf_to_field = montauk.wcs.make_batoid_wcs(
    obsdata=obsdata, dm_detector=dm_detector,
)

sky_model = imsim.SkyModel(
    exptime=obsdata['vistime'],
    mjd=obsdata['mjd'],
    bandpass=obsdata['bandpass'],
)

gradient = montauk.sky.FixedSkyGradient(sky_model)
vignetter = montauk.vignetting.Vignetter(dm_detector)

sensor = montauk.sensor.make_sensor(
    dm_detector=dm_detector, gs_rng=gs_rng,
)

dcr = montauk.dcr.DCRMaker(
    bandpass=obsdata['bandpass'],
    hour_angle=obsdata['HA'],
)

optics = montauk.optics.OpticsMaker(
    altitude=obsdata['altitude'],
    azimuth=obsdata['azimuth'],
    boresight=obsdata['boresight'],
    rot_tel_pos=obsdata['rotTelPos'],
    band=obsdata['band'],
    dm_detector=dm_detector,
    wcs=wcs,
    icrf_to_field=icrf_to_field,
)

photon_ops_maker = montauk.photon_ops.PhotonOpsMaker(
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

artist = montauk.artist.Artist(
    bandpass=obsdata['bandpass'],
    sensor=sensor,
    photon_ops_maker=photon_ops_maker,
    diffraction_fft=diffraction_fft,
    gs_rng=gs_rng,
)

psf = montauk.psfws.make_psfws_psf(
    obsdata=obsdata,
    gs_rng=gs_rng,
    psf_config=psf_config,
)

cosmics = montauk.cosmic_rays.CosmicRays(
    cosmic_ray_rate=cosmic_ray_rate,
    exptime=obsdata['vistime'],
    gs_rng=gs_rng,
)

image, sky_image, truth = montauk.runner.run_sim(
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
```
