import galsim
import imsim
import mimsim
from mimsim.utils import sigma_clip
import numpy as np
from tqdm import trange
import pytest


@pytest.mark.parametrize('sed_type', ['narrow', 'trivial'])
def test_wcs_with_dcr(sed_type):
    mimsim.logging.setup_logging('info')

    seed = 9119

    altitude = 25
    nobj = 200
    stamp_size = 33
    n_photons = 10_000
    rng = np.random.default_rng(seed)
    gs_rng = galsim.BaseDeviate(rng.integers(0, 2**30))

    band = 'g'
    detnum = 35

    obsdata = mimsim.simtools.load_example_obsdata(
        band=band,
        altitude=altitude,
    )

    bp = obsdata['bandpass']
    psf = galsim.Gaussian(fwhm=0.6)

    # SED not centered at the effective wavelength
    if sed_type == 'narrow':
        sed = mimsim.seds.get_narrow_sed(
            wavemin=bp.blue_limit,
            wavemax=bp.red_limit,
            wave0=bp.blue_limit + (bp.red_limit - bp.blue_limit) / 4,
            wavesig=5,
            npts=1000,
        )
    else:
        sed = mimsim.seds.get_trivial_sed()

    dm_detector = mimsim.camera.make_dm_detector(detnum)

    # dcr will be relative to the effective wavelength of the bandpass
    dcr = mimsim.dcr.DCRMaker(
        bandpass=obsdata['bandpass'],
        hour_angle=obsdata['HA'],
    )

    tree_rings = mimsim.tree_rings.make_tree_rings([detnum])
    sensor = mimsim.sensor.make_sensor(
        dm_detector=dm_detector,
        tree_rings=tree_rings,
        gs_rng=gs_rng,
    )

    wcs, icrf_to_field = mimsim.wcs.make_batoid_wcs(
        obsdata=obsdata, dm_detector=dm_detector,
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

    obj = galsim.Gaussian(fwhm=0.2, flux=n_photons) * sed

    ras = np.zeros(nobj)
    decs = np.zeros(nobj)
    xmeans = np.zeros(nobj)
    ymeans = np.zeros(nobj)

    for i in trange(nobj):

        # generate random x, y and convert to ra, dec using "wrong wcs"
        x = rng.uniform(low=1, high=4096)
        y = rng.uniform(low=1, high=4096)

        sky_pos = wcs.toWorld(galsim.PositionD(x, y))
        image_pos = wcs.toImage(sky_pos)
        local_wcs = wcs.local(image_pos)

        pos = artist.get_pos(
            obj=obj, obj_coord=sky_pos, image_pos=image_pos,
            stamp_size=stamp_size, local_wcs=local_wcs, psf=psf,
        )

        ras[i] = sky_pos.ra.rad
        decs[i] = sky_pos.dec.rad

        xmeans[i] = pos['x']
        ymeans[i] = pos['y']

    nuse = xmeans.size // 2

    new_wcs = galsim.FittedSIPWCS(
        xmeans[:nuse], ymeans[:nuse], ras[:nuse], decs[:nuse], order=3,
    )
    xcheck, ycheck = new_wcs.radecToxy(
        ra=ras[nuse:], dec=decs[nuse:],
        units=galsim.radians,
    )
    xcoff = xmeans[nuse:] - xcheck
    ycoff = ymeans[nuse:] - ycheck

    xcoff_mean, _, xcoff_err = sigma_clip(xcoff)
    ycoff_mean, _, ycoff_err = sigma_clip(ycoff)
    check_range('xcoff', xcoff_mean, xcoff_err)
    check_range('ycoff', ycoff_mean, ycoff_err)


def check_range(name, val, err, nsigma=3.5):
    low = val - nsigma * err
    high = val + nsigma * err
    mess = f'{low:.4f} < {name} < {high:.4f}'

    if low > 0 or high < 0:
        raise ValueError(mess)
    else:
        print(mess)


def print_range(name, val, err, nsigma=3.5):
    low = val - nsigma * err
    high = val + nsigma * err
    print(f'{low:.4f} < {name} < {high:.4f}')


if __name__ == '__main__':
    test_wcs_with_dcr('narrow')
    test_wcs_with_dcr('trivial')
