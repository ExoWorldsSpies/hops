
import sys
import numpy as np
import twirl
import hops.pylightcurve41 as plc

from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import fit_wcs_from_points, pixel_to_skycoord
from photutils.aperture import CircularAperture, aperture_photometry

from .centroids_and_stars import _find_centroids, _star_from_centroid, _separation, _get_gaia_stars


def image_mean_std(fits_data,
                   samples=10000, mad_filter=5.0):

    test_data = 1.0 * fits_data.flatten()
    test_data = test_data[::int(len(test_data)/samples + 1)]
    quantiles = np.quantile(1.0 * test_data, [0.5 - 0.341, 0.5, 0.5 + 0.341])

    if quantiles[0] == quantiles[1]:
        if quantiles[1] == quantiles[2]:
            test2 = test_data[np.where(test_data > quantiles[2])]
            quantiles2 = np.quantile(test2, [0.5 - 0.341, 0.5, 0.5 + 0.341])
            return quantiles[1], quantiles2[2] - quantiles[1]
        else:
            return quantiles[1], quantiles[2] - quantiles[1]

    else:
        try:
            if np.sum((test_data/64.0-np.int_(test_data/64.0)) == 0) > len(test_data)/2:
                distribution = plc.one_d_distribution(fits_data/64.0, samples=samples, gaussian_fit=True,
                                                      mad_filter=mad_filter)
                return 64.0 * distribution[2], 64.0 * distribution[3]
            elif np.sum((test_data/16.0-np.int_(test_data/16.0)) == 0) > len(test_data)/2:
                distribution = plc.one_d_distribution(fits_data/16.0, samples=samples, gaussian_fit=True,
                                                      mad_filter=mad_filter)
                return 16.0 * distribution[2], 16.0 * distribution[3]
            elif np.sum((test_data/4.0-np.int_(test_data/4.0)) == 0) > len(test_data)/2:
                distribution = plc.one_d_distribution(fits_data/4.0, samples=samples, gaussian_fit=True,
                                                      mad_filter=mad_filter)
                return 4.0 * distribution[2], 4.0 * distribution[3]
            else:
                distribution = plc.one_d_distribution(fits_data, samples=samples, gaussian_fit=True,
                                                      mad_filter=mad_filter)
                return distribution[2], distribution[3]
        except Exception as e:
            import traceback
            print(f'{traceback.format_exc()}')
            return quantiles[1], quantiles[2] - quantiles[1]


def image_burn_limit(fits_header, key=None):

    if key and key in fits_header:
        return(fits_header[key])

    else:
        return min(65535, int(2 ** abs(fits_header['BITPIX']) - 1))


def image_psf(fits_data, fits_header,
              mean=None, std=None, burn_limit=None,
              sample=10, psf_guess=2):

    if mean is None or std is None:
        mean, std = image_mean_std(fits_data)

    if burn_limit is None:
        burn_limit = image_burn_limit(fits_header)

    # t0 = time.time()

    centroids = []
    snr = 54

    while len(centroids) < sample * 2 and snr >= 4.0:
        centroids = _find_centroids(fits_data, -10, 10**10, -10, 10**10, mean, std, 0.9 * burn_limit, psf_guess, snr)
        snr -= 10

    # print('1', 1000*(time.time() - t0))
    # t0 = time.time()

    psf = []
    psf_err = []

    for centroid in centroids:
        star = _star_from_centroid(fits_data, centroid[1], centroid[2], mean, std, burn_limit, psf_guess)

        if star:

            psf.append(star[0][4])
            psf_err.append(np.sqrt(star[1][4][4]))

        if len(psf) == sample:
            break

    # print('2', 1000*(time.time() - t0))
    # t0 = time.time()

    psf = plc.waverage(np.array(psf), np.array(psf_err))[0]

    # print('3', 1000*(time.time() - t0))

    return psf


def image_find_stars(fits_data, fits_header, x_low=0, x_upper=None, y_low=0, y_upper=None, x_centre=None, y_centre=None,
                     mean=None, std=None, burn_limit=None, psf=None,
                     snr=4.0, psf_variation_allowed=0.5,
                     aperture=3, sky_inner_aperture=1.7, sky_outer_aperture=2.4,
                     order_by_flux=True, star_limit=None,
                     progress_window=None, verbose=False):

    if verbose:
        print('\nAnalysing frame...')

    if mean is None or std is None:
        mean, std = image_mean_std(fits_data)

    if burn_limit is None:
        burn_limit = image_burn_limit(fits_header)

    if psf is None:
        psf = image_psf(fits_data, fits_header, mean, std, burn_limit)

    if x_upper is None:
        x_upper = fits_data.shape[1]

    if y_upper is None:
        y_upper = fits_data.shape[0]

    if x_centre is None:
        x_centre = fits_data.shape[1] / 2

    if y_centre is None:
        y_centre = fits_data.shape[0] / 2

    fits_data_size_y, fits_data_size_x = fits_data.shape

    centroids = _find_centroids(fits_data, x_low, x_upper, y_low, y_upper, mean, std, burn_limit, psf, snr)

    if progress_window:
        if progress_window.exit:
            return None

    stars = []
    if verbose:
        print('Verifying stars...')

    for num, centroid in enumerate(centroids):

        star = _star_from_centroid(fits_data, centroid[1], centroid[2], mean, std, burn_limit, psf, snr)

        if star:

            if abs(star[0][4] - psf) < (3 * np.sqrt(star[1][4][4]) + psf * psf_variation_allowed):

                if len(stars) == 0 or np.min([np.sqrt((ff[0]-star[0][2])**2 + (ff[1]-star[0][3])**2) for ff in stars]) > psf:

                    x_mean, y_mean = star[0][2], star[0][3]
                    ap = aperture * star[0][4]

                    if x_mean - ap < 0 or x_mean + ap > fits_data_size_x or y_mean - ap < 0 or y_mean + ap > fits_data_size_y:
                        pass
                    else:
                        total_flux = aperture_photometry(fits_data, CircularAperture(np.array([x_mean-0.5, y_mean-0.5]),
                                                                                     ap))['aperture_sum'][0]
                        total_area = np.pi * ap * ap

                        sky_area_1 = ap * sky_inner_aperture
                        sky_area_2 = ap * sky_outer_aperture
                        sky_area_1, sky_area_2 = min(sky_area_2, sky_area_1), max(sky_area_2, sky_area_1)

                        y_min = int(max(int(y_mean - sky_area_2), 0))
                        y_max = int(min(int(y_mean + sky_area_2), len(fits_data) - 1))
                        x_min = int(max(int(x_mean - sky_area_2), 0))
                        x_max = int(min(int(x_mean + sky_area_2), len(fits_data[0]) - 1))
                        datax, datay = np.meshgrid(np.arange(x_min, x_max + 1) + 0.5,
                                                   np.arange(y_min, y_max + 1) + 0.5)

                        sky_area = fits_data[y_min:y_max + 1, x_min:x_max + 1]
                        sky_area = np.ones_like(sky_area) * sky_area

                        sky_area = sky_area[np.where((((datax - x_mean)**2 + (datay - y_mean)**2) > sky_area_1**2) * (sky_area < star[0][1] + 3 * std))]

                        sky, sky_error = plc.mean_std_from_median_mad(sky_area)

                        sky_flux = sky * total_area
                        sky_flux_unc = np.sqrt(total_area * sky_error * sky_error)

                        target_flux = total_flux - sky_flux

                        stars.append([
                            star[0][2], star[0][3], star[0][0], star[0][1], star[0][4], star[0][5],
                            total_flux, target_flux, sky_flux, sky_flux_unc
                        ])

            if verbose:
                sys.stdout.write('\r\033[K')
                sys.stdout.write('{0}/{1} '.format(num + 1, len(centroids)))
                sys.stdout.flush()

        if star_limit and len(stars) >= star_limit:
            break

    if verbose:
        print('')

    if len(stars) > 0:

        if order_by_flux is True:
            stars = sorted(stars, key=lambda x: -x[7])
        else:

            stars = sorted(stars, key=lambda x: np.sqrt((x[0] - x_centre) ** 2 + (x[1] - y_centre) ** 2))

        return stars

    else:
        return None


def image_plate_solve(fits_data, fits_header, ra, dec,
                      mean=None, std=None, burn_limit=None, psf=None, stars=None, n=20, pixel=None,
                      progress_window=None, verbose=False, gaia_engine=_get_gaia_stars):

    if verbose:
        print('\nAnalysing frame...')

    if mean is None or std is None:
        mean, std = image_mean_std(fits_data)

    if burn_limit is None:
        burn_limit = image_burn_limit(fits_header)

    if psf is None:
        psf = image_psf(fits_data, fits_header, mean, std, burn_limit)

    if type(stars) == type(None):
        stars = image_find_stars(fits_data, fits_header,
                                 mean=mean, std=std, burn_limit=burn_limit, psf=psf,
                                 progress_window=progress_window, verbose=verbose)

    stars = np.array(stars)

    if pixel is None:
        pixel = 2.0 / psf
    default_wcs = wcs.WCS(naxis=2)
    default_wcs.wcs.ctype = ["RA---ARC", "DEC--ARC"]
    default_wcs.wcs.crpix=np.array(fits_data.shape)[::-1]/2
    default_wcs.wcs.crval=np.array([ra, dec])
    default_wcs.wcs.pc[0][0] = -pixel/60/60
    default_wcs.wcs.pc[1][1] = -pixel/60/60

    ra1, dec1 = default_wcs.all_pix2world([[0, 0]], 0)[0]
    ra2, dec2 = default_wcs.all_pix2world([[0, len(fits_data)]], 0)[0]
    ra3, dec3 = default_wcs.all_pix2world([[len(fits_data[0]), len(fits_data)]], 0)[0]
    ra4, dec4 = default_wcs.all_pix2world([[len(fits_data[0]), 0]], 0)[0]

    fov_radius = max([_separation(ra, dec, ra1, dec1), _separation(ra, dec, ra2, dec2),
                      _separation(ra, dec, ra3, dec3), _separation(ra, dec, ra4, dec4)])

    if verbose:
        print('Initial FOV radius guess: ', fov_radius)

    x0, y0 = len(fits_data[0])/2, len(fits_data)/2

    gaia_query = gaia_engine(ra, dec, 2.0 * fov_radius, 9 * len(stars))
    gaia_stars = np.array([[star['ra'], star['dec']] for star in gaia_query])

    xx = np.array(default_wcs.wcs_world2pix(gaia_stars[:,0], gaia_stars[:,1], 1))

    gaia_stars_in = gaia_stars[np.where(np.sqrt((xx[0] - x0)**2 + (xx[1] - y0)**2) < min(x0, y0))]
    gaia_stars_in = gaia_stars_in[:2*n]

    projected_gaia_stars_in = np.array(default_wcs.wcs_world2pix(gaia_stars_in[:,0], gaia_stars_in[:,1], 1)).T

    detected_stars = np.array([[star[0], star[1]] for star in stars])
    detected_stars_in = detected_stars[np.where(np.sqrt((detected_stars[:,0] - x0)**2 + (detected_stars[:,1] - y0)**2) < min(x0, y0))]
    detected_stars_in = detected_stars_in[:2*n]

    tolerance = 3 * psf
    X = twirl.utils.find_transform(projected_gaia_stars_in, detected_stars_in, n=n, tolerance=tolerance)

    projected_gaia_stars = np.array(default_wcs.wcs_world2pix(gaia_stars[:,0], gaia_stars[:,1], 1)).T
    transformed_projected_gaia_stars = twirl.utils.affine_transform(X)(projected_gaia_stars)

    central_transformed_projected_gaia_star = gaia_stars[sorted(
        np.arange(len(transformed_projected_gaia_stars)),
        key=lambda x: np.sqrt(
            (transformed_projected_gaia_stars[x][0] - len(fits_data[0])/2)**2 +
            (transformed_projected_gaia_stars[x][1] - len(fits_data)/2)**2)
    )[0]]

    s1, s2 = twirl.utils.cross_match(
        transformed_projected_gaia_stars,
        detected_stars,
        return_ixds=True, tolerance=tolerance).T

    plate_solution = fit_wcs_from_points(
        detected_stars[s2].T,
        SkyCoord(gaia_stars[s1], unit="deg"),
        proj_point = SkyCoord(*(central_transformed_projected_gaia_star + 0.01), unit="deg"),
        sip_degree=3,
    )

    ra1, dec1 = plate_solution.all_pix2world([[0, 0]], 0)[0]
    ra2, dec2 = plate_solution.all_pix2world([[0, len(fits_data)]], 0)[0]
    ra3, dec3 = plate_solution.all_pix2world([[len(fits_data[0]), len(fits_data)]], 0)[0]
    ra4, dec4 = plate_solution.all_pix2world([[len(fits_data[0]), 0]], 0)[0]

    fov_radius = max([_separation(ra, dec, ra1, dec1), _separation(ra, dec, ra2, dec2),
                      _separation(ra, dec, ra3, dec3), _separation(ra, dec, ra4, dec4)])

    if verbose:
        print('Final FOV radius: ', fov_radius)

    gaia_query = gaia_engine(ra, dec, fov_radius, 8 * len(stars))
    gaia_stars = np.array([[star['ra'], star['dec']] for star in gaia_query])

    transformed_projected_gaia_stars = np.array(plate_solution.wcs_world2pix(gaia_stars[:,0], gaia_stars[:,1], 1)).T

    s1, s2 = twirl.utils.cross_match(
        transformed_projected_gaia_stars,
        detected_stars,
        return_ixds=True, tolerance=tolerance).T

    plate_solution = fit_wcs_from_points(
        detected_stars[s2].T,
        SkyCoord(gaia_stars[s1], unit="deg"),
        proj_point = SkyCoord(*(central_transformed_projected_gaia_star + 0.01), unit="deg"),
        sip_degree=3,
    )

    source_id_key = 'source_id'
    if source_id_key not in gaia_query.keys():
        source_id_key = 'SOURCE_ID'

    identified_stars = {
        'plate_solution': plate_solution,
        'identified_stars': {
            gaia_query[source_id_key][s1][nn]:{
                'ra': gaia_query['ra'][s1][nn],
                'dec': gaia_query['dec'][s1][nn],
                'phot_g_mean_mag': gaia_query['phot_g_mean_mag'][s1][nn],
                'phot_bp_mean_mag': gaia_query['phot_bp_mean_mag'][s1][nn],
                'phot_rp_mean_mag': gaia_query['phot_rp_mean_mag'][s1][nn],
                'x': stars[s2][nn][0],
                'y': stars[s2][nn][1],
                'peak': stars[s2][nn][2],
                'floor': stars[s2][nn][3],
                'x-sigma': stars[s2][nn][4],
                'y-sigma': stars[s2][nn][5],
                'total_flux': stars[s2][nn][6],
                'target_flux': stars[s2][nn][7],
                'sky_flux': stars[s2][nn][8],
                'sky_flux_unc': stars[s2][nn][9],
            }
            for nn in range(len(s1))
        }
    }

    return identified_stars


def bin_frame(data_frame, binning):

    binning = int(binning)

    if binning <= 1:
        return data_frame

    new_frame = []
    for xx in range(len(data_frame))[::binning]:
        new_frame.append(list(np.sum(data_frame[xx:xx + binning], 0)))

    new_frame = np.swapaxes(np.array(new_frame), 0, 1)

    new_frame2 = []
    for xx in range(len(new_frame))[::binning]:
        new_frame2.append(list(np.sum(new_frame[xx:xx + binning], 0)))

    data_frame = np.swapaxes(np.array(new_frame2), 0, 1)

    return data_frame


