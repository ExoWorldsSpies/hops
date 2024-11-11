
import numpy as np
import warnings
import hops.pylightcurve41 as plc


from scipy.optimize import minimize
from astroquery.gaia import Gaia


def _find_centroids(data_array, x_low, x_upper, y_low, y_upper, mean, std, burn_limit, psf, snr=4):

    psf = max(1, int(round(psf)))
    limit = mean + snr * std

    x_upper = int(min(x_upper, len(data_array[0])))
    y_upper = int(min(y_upper, len(data_array)))
    x_low = int(max(0, x_low))
    y_low = int(max(0, y_low))

    data_array = np.full_like(data_array[y_low:y_upper + 1, x_low:x_upper + 1],
                              data_array[y_low:y_upper + 1, x_low:x_upper + 1])

    bright = np.where(data_array[psf:-psf, psf:-psf] > limit)
    bright = (bright[0] + psf, bright[1] + psf)

    test = []
    for i in range(-psf, psf + 1):
        for j in range(-psf, psf + 1):
            test.append(data_array[bright[0] + i, bright[1]+j])

    test = np.array(test)
    data_array_test = data_array[bright]

    min_test = np.sum(test > mean + snr * std, 0) >= 0.5 * (2 * psf + 1)**2
    max_test = np.max(test, 0)

    centroids = np.where((max_test < burn_limit) * (max_test == data_array_test) * min_test)[0]
    centroids = np.swapaxes([data_array_test[centroids], bright[1][centroids] + x_low, bright[0][centroids] + y_low], 0, 1)
    centroids = np.int_(centroids)
    centroids = np.array(sorted(centroids, key=lambda x: -x[0]))

    del test
    del data_array_test

    # import matplotlib.pyplot as plt
    # import matplotlib.patches as mpatches
    # fig = plt.figure()
    # ax = fig.add_subplot(1,1,1)
    # plt.imshow(data_array, vmax=mean+3*std, origin='lower', extent=(x_low, x_upper, y_low, y_upper))
    # for centroid in stars:
    #     patch = mpatches.Circle((centroid[1], centroid[2]), 5, ec='r', fill=False)
    #     ax.add_patch(patch)
    # plt.savefig('test.pdf')
    # plt.show()

    return centroids


def two_d_gaussian(x_array, y_array, model_norm, model_floor, model_x_mean, model_y_mean, model_x_sigma, model_y_sigma, model_theta):

    xt_array = x_array - model_x_mean
    yt_array = y_array - model_y_mean
    coss = np.cos(model_theta)
    sinn = np.sin(model_theta)

    return model_floor + model_norm * np.exp(-0.5 * (((-xt_array * sinn + yt_array * coss) / model_y_sigma) ** 2 + ((xt_array * coss + yt_array * sinn) / model_x_sigma) ** 2))


def _star_from_centroid(data_array, centroid_x, centroid_y, mean, std, burn_limit, psf, snr=4, search_window=10):

    star = None
    try:
        # t0 = time.time()

        bright = data_array[centroid_y][centroid_x]

        search_window = int(round(search_window * psf))
        y_min = int(max(int(centroid_y) - search_window, 0))
        y_max = int(min(int(centroid_y) + search_window + 1, len(data_array)))
        x_min = int(max(int(centroid_x) - search_window, 0))
        x_max = int(min(int(centroid_x) + search_window + 1, len(data_array[0])))

        datax, datay = np.meshgrid(np.arange(x_min, x_max) + 0.5,
                                   np.arange(y_min, y_max) + 0.5)

        datax = datax.flatten()
        datay = datay.flatten()
        dataz = data_array[y_min: y_max, x_min: x_max].flatten()
        datae = np.sqrt(np.abs(dataz) + 1)

        # print('0', 1000*(time.time() - t0))
        # t0 = time.time()

        initials = np.array([bright - mean, mean, centroid_x + 0.5, centroid_y + 0.5, psf, psf, 0])
        bounds_1 = np.array([0,             mean - 10 * std, x_min,            y_min,            0,              0,              -np.pi/2])
        bounds_2 = np.array([np.max(dataz), mean + 10 * std, x_max, y_max, 100 * psf, 100 * psf, np.pi / 2])

        # print('1', 1000*(time.time() - t0))
        # t0 = time.time()

        def function_to_fit(xy_array, model_norm, model_floor, model_x_mean, model_y_mean, model_x_sigma, model_y_sigma, model_theta):
            return two_d_gaussian(datax, datay, model_norm, model_floor, model_x_mean, model_y_mean, model_x_sigma, model_y_sigma, model_theta)

        popt, pcov = plc.curve_fit(function_to_fit, [0], dataz, p0=initials, maxfev=int(1600/(snr**2)),
                                   sigma=datae,
                                   bounds=(np.array(bounds_1), np.array(bounds_2))
                                   )

        # print(popt, pcov)

        # print('2', 1000*(time.time() - t0))
        # t0 = time.time()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if popt[0] > snr * std:
                if popt[0] < burn_limit:
                    if np.nan not in [np.sqrt(abs(pcov[ff][ff])) for ff in range(len(pcov) - 1)]:
                        if np.inf not in [np.sqrt(abs(pcov[ff][ff])) for ff in range(len(pcov) - 1)]:
                            if 0 not in [np.sqrt(abs(pcov[ff][ff])) for ff in range(len(pcov) - 1)]:
                                star = (popt, pcov)
            #             else:
            #                 print('Error not estimated.')
            #         else:
            #             print('Error not estimated.')
            #     else:
            #         print('Star saturated.')
            # else:
            #     print('Low peak')

        if popt[5] > popt[4]:
            popt[4], popt[5] = popt[5], popt[4]
            pcov[4][4], pcov[5][5] = pcov[5][5], pcov[4][4]
            popt[6] -= np.pi/2

        # print('3', 1000*(time.time() - t0))

    except Exception as e:
        # print(e)
        pass

    return star


def _separation(ra1, dec1, ra2, dec2):
    return (180.0/np.pi) * np.arccos(np.minimum(1, np.sin(dec1* np.pi / 180.0) * np.sin(dec2* np.pi / 180.0) +
                                                np.cos(dec1* np.pi / 180.0) * np.cos(dec2* np.pi / 180.0) * np.cos(ra1* np.pi / 180.0 - ra2* np.pi / 180.0)))


def _get_gaia_stars(ra_0, dec_0, radius, limit=100):

    Gaia.ROW_LIMIT = limit
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"

    query = "SELECT TOP {0} {1} FROM {2} WHERE 1 = CONTAINS(POINT('ICRS', {3}, {4}),CIRCLE('ICRS', ra, dec, {5})) ORDER BY phot_g_mean_mag ASC".format(
        limit,
        'source_id, ra, dec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag',
        Gaia.MAIN_GAIA_TABLE,
        ra_0,
        dec_0,
        radius
    )
    job = Gaia.launch_job_async(query)

    return job.get_results()
