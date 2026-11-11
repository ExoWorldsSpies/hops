"""A small synthetic data set for the regression tests.

12 light frames (900 x 1200 pixels) of a star field with a transiting planet, drift and seeing
changes, plus bias, dark and flat frames. Some flat pixels are dead (also in the first two rows),
so the master flat cleaning is exercised. Everything is generated at test time from a fixed seed;
no image data is stored in the repository.
"""

import os
import datetime
import numpy as np

from astropy.io import fits as pf
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u

TARGET_RA = '11:00:40.18'
TARGET_DEC = '+64:57:50.3'
SITE_LAT = '+48:05:17.15'
SITE_LONG = '+15:45:24.40'
SITE_LAT_DEG = 48 + 5 / 60 + 17.15 / 3600
SITE_LONG_DEG = 15 + 45 / 60 + 24.40 / 3600

SHAPE = (900, 1200)
N_FRAMES = 12
CADENCE = 180.0
EXPTIME = 120.0
START = datetime.datetime(2026, 3, 10, 21, 0, 0)

TARGET_XY = (600.3, 449.6)
COMPARISON_OFFSETS = [(-180.0, 95.0), (210.0, -120.0), (-250.0, -160.0), (150.0, 170.0)]

# planet (KPS-1b like), claret limb-darkening coefficients of a sun-like star in R
PLANET = {'rp_over_rs': 0.1145, 'period': 1.706291, 'sma_over_rs': 6.37, 'eccentricity': 0.0,
          'inclination': 83.2, 'periastron': 0.0}
LDC = [0.55, -0.15, 0.62, -0.25]


def _frame_times():
    starts = [START + datetime.timedelta(seconds=i * CADENCE) for i in range(N_FRAMES)]
    mid = Time([s + datetime.timedelta(seconds=EXPTIME / 2) for s in starts], scale='utc')
    target = SkyCoord(TARGET_RA, TARGET_DEC, unit=(u.hourangle, u.deg))
    site = EarthLocation.from_geodetic(SITE_LONG_DEG * u.deg, SITE_LAT_DEG * u.deg, 300 * u.m)
    bjd = mid.tdb + mid.light_travel_time(target, location=site)
    return starts, bjd.jd


def _stamp(image, x, y, flux, sigma):
    r = int(5 * sigma) + 2
    x0, y0, x1, y1 = int(x) - r, int(y) - r, int(x) + r + 1, int(y) + r + 1
    xs, ys = slice(max(x0, 0), min(x1, image.shape[1])), slice(max(y0, 0), min(y1, image.shape[0]))
    xx, yy = np.meshgrid(np.arange(xs.start, xs.stop) + 0.5, np.arange(ys.start, ys.stop) + 0.5)
    image[ys, xs] += flux / (2 * np.pi * sigma ** 2) * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * sigma ** 2))


def _write(path, data, header):
    hdu = pf.PrimaryHDU(np.clip(np.round(data), 0, 65535).astype(np.uint16))
    for key, value in header.items():
        hdu.header[key] = value
    hdu.writeto(path, overwrite=True)


def generate(directory, seed=7):
    """Write light, bias, dark and flat frames into `directory` (deterministic for a seed)."""

    import hops.pylightcurve41 as plc

    os.makedirs(directory, exist_ok=True)
    rng = np.random.default_rng(seed)
    ny, nx = SHAPE

    yy, xx = np.mgrid[0:ny, 0:nx]
    rr = ((xx - nx / 2) ** 2 + (yy - ny / 2) ** 2) / ((nx / 2) ** 2 + (ny / 2) ** 2)
    flat_field = (1 - 0.12 * rr) * (1 + 0.01 * rng.standard_normal(SHAPE))
    bias_level = 300 + 3 * np.sin(xx / 37.0)
    dark_rate = 0.4 + 0.05 * rng.random(SHAPE)
    read_noise = 8.0

    n_stars = 140
    stars_x = rng.uniform(30, nx - 30, n_stars)
    stars_y = rng.uniform(30, ny - 30, n_stars)
    stars_f = 10 ** rng.uniform(3.3, 5.6, n_stars)
    keep = np.hypot(stars_x - TARGET_XY[0], stars_y - TARGET_XY[1]) > 40
    for dx, dy in COMPARISON_OFFSETS:
        keep &= np.hypot(stars_x - TARGET_XY[0] - dx, stars_y - TARGET_XY[1] - dy) > 40
    stars_x, stars_y, stars_f = stars_x[keep], stars_y[keep], stars_f[keep]

    starts, bjd = _frame_times()
    mid_time = float(bjd[0] + 0.55 * (bjd[-1] - bjd[0]))
    relative_flux = plc.transit(LDC, PLANET['rp_over_rs'], PLANET['period'], PLANET['sma_over_rs'],
                                PLANET['eccentricity'], PLANET['inclination'], PLANET['periastron'],
                                mid_time, bjd, method='claret')

    common = {'OBJCTRA': TARGET_RA.replace(':', ' '), 'OBJCTDEC': TARGET_DEC.replace(':', ' '),
              'SITELAT': SITE_LAT, 'SITELONG': SITE_LONG, 'FILTER': 'R'}

    def noisy(electrons):
        return rng.poisson(np.clip(electrons, 0, None)).astype(float) + read_noise * rng.standard_normal(SHAPE)

    for i, start in enumerate(starts):
        dx = 0.12 * i + rng.normal(0, 0.4)
        dy = -0.07 * i + rng.normal(0, 0.4)
        sigma = 1.8 + 0.25 * np.sin(i / 6.0) + rng.normal(0, 0.03)
        image = np.full(SHAPE, (120 + 15 * i / N_FRAMES) * EXPTIME / 60.0)
        for x, y, f in zip(stars_x, stars_y, stars_f):
            _stamp(image, x + dx, y + dy, f, sigma)
        _stamp(image, TARGET_XY[0] + dx, TARGET_XY[1] + dy, 1.6e5 * relative_flux[i], sigma)
        for (ox, oy), f in zip(COMPARISON_OFFSETS, [1.4e5, 1.75e5, 1.5e5, 1.9e5]):
            _stamp(image, TARGET_XY[0] + ox + dx, TARGET_XY[1] + oy + dy, f, sigma)
        raw = bias_level + noisy(image * flat_field + dark_rate * EXPTIME)
        header = dict(common, EXPTIME=EXPTIME, **{'DATE-OBS': start.isoformat(timespec='milliseconds')})
        _write(os.path.join(directory, 'synthetic_light_{0:04d}.fits'.format(i + 1)), raw, header)

    # dead pixels in every flat (two of them in rows 0 and 1) give bad pixels in the master flat
    dead_y, dead_x = rng.integers(0, ny, 60), rng.integers(1, nx - 1, 60)
    dead_y[:3] = [0, 1, 1]

    for i in range(5):
        _write(os.path.join(directory, 'synthetic_bias_{0:04d}.fits'.format(i + 1)),
               bias_level + read_noise * rng.standard_normal(SHAPE), dict(common, EXPTIME=0.0))
        _write(os.path.join(directory, 'synthetic_dark_{0:04d}.fits'.format(i + 1)),
               bias_level + noisy(dark_rate * EXPTIME), dict(common, EXPTIME=EXPTIME))
        flat = bias_level + noisy(20000 * flat_field + dark_rate * 5.0)
        flat[dead_y, dead_x] = 0
        _write(os.path.join(directory, 'synthetic_flat_{0:04d}.fits'.format(i + 1)), flat, dict(common, EXPTIME=5.0))


def configure(log, directory):
    """The settings the Data & Target window would save for the synthetic data set."""

    os.chdir(directory)
    log.set_param('directory', directory)
    log.set_param('directory_short', os.path.split(directory)[1])
    log.load_main_log()
    for key, value in {
        'observation_files': 'light', 'bias_files': 'bias', 'dark_files': 'dark', 'darkf_files': '---',
        'flat_files': 'flat', 'exposure_time_key': 'EXPTIME', 'observation_date_key': 'DATE-OBS',
        'observation_time_key': 'DATE-OBS', 'time_stamp': 'exposure start',
        'target_ra_dec': '{0} {1}'.format(TARGET_RA, TARGET_DEC), 'target_name': 'Synthetic-1',
        'location': '{0} {1}'.format(SITE_LAT, SITE_LONG), 'filter': 'R', 'bin_fits': 1,
        'crop_x1': 0, 'crop_x2': 0, 'crop_y1': 0, 'crop_y2': 0, 'crop_edge_pixels': 0,
        'faint_target_mode': False, 'moving_target_mode': False, 'rotating_field_mode': False,
        'colour_camera_mode': False, 'data_target_complete': True,
    }.items():
        log.set_param(key, value)
    log.set_param('data_target_version', log.version)
    log.save_local_log()
