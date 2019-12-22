from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np

from .analysis_functions_and_fit import fit_two_d_gaussian
from .analysis_distributions import one_d_distribution
from .tools_maths import waverage


def _star_from_centroid(data_array, centroid_x, centroid_y, mean, std, burn_limit, star_std, std_limit):

    star = None
    try:
        search_window = int(round(5 * star_std))
        y_min = int(max(int(centroid_y) - search_window, 0))
        y_max = int(min(int(centroid_y) + search_window, len(data_array) - 1))
        x_min = int(max(int(centroid_x) - search_window, 0))
        x_max = int(min(int(centroid_x) + search_window, len(data_array[0]) - 1))

        datax, datay = np.meshgrid(np.arange(x_min, x_max + 1) + 0.5,
                                   np.arange(y_min, y_max + 1) + 0.5)

        dataz = data_array[y_min: y_max + 1, x_min: x_max + 1]
        popt, pcov = fit_two_d_gaussian(datax, datay, dataz, positive=True, maxfev=1000)

        if (popt[0] > mean + std_limit * std and popt[0] + popt[1] < burn_limit
                and pcov[0][0] > 0 and popt[0] > std_limit * np.sqrt(pcov[0][0])):

                star = (popt, pcov)

    except:
        pass

    return star


def find_single_star(data_array, predicted_x, predicted_y, mean=None, std=None, burn_limit=50000, star_std=2,
                     std_limit=3.0):

    if mean is None or std is None:
        try:
            fit_mean, fit_std = one_d_distribution(data_array, gaussian_fit=True)[2:4]
        except:
            fit_mean = np.mean(data_array)
            fit_std = np.std(data_array)

        if not mean:
            mean = fit_mean

        if not std:
            std = fit_std

    centroids = find_centroids(data_array, predicted_x - 5 * star_std, predicted_x + 5 * star_std,
                               predicted_y - 5 * star_std, predicted_y + 5 * star_std, mean, std, burn_limit, star_std,
                               std_limit)

    centroids = sorted(centroids, key=lambda x: np.sqrt((x[0] - predicted_x) ** 2 + (x[1] - predicted_y) ** 2))

    star = None
    for centroid in centroids:
        star = _star_from_centroid(data_array, centroid[0], centroid[1], mean, std, burn_limit, star_std, std_limit)
        if star:
            star = [star[0][2], star[0][3], star[0][0], star[0][1], star[0][4], star[0][5], centroid[0], centroid[1]]
            break

    return star


def find_all_stars(data_array, x_low=0, x_upper=None, y_low=0, y_upper=None, x_centre=None, y_centre=None,
                   mean=None, std=None, burn_limit=50000, star_std=2, std_limit=3.0, order_by_flux=False):

    if mean is None or std is None:
        try:
            fit_mean, fit_std = one_d_distribution(data_array, gaussian_fit=True)[2:4]
        except:
            fit_mean = np.mean(data_array)
            fit_std = np.std(data_array)

        if not mean:
            mean = fit_mean

        if not std:
            std = fit_std

    if x_upper is None:
        x_upper = data_array.shape[1]

    if y_upper is None:
        y_upper = data_array.shape[0]

    if x_centre is None:
        x_centre = data_array.shape[1] / 2

    if y_centre is None:
        y_centre = data_array.shape[0] / 2

    centroids = find_centroids(data_array, x_low, x_upper, y_low, y_upper, mean, std, burn_limit, star_std, std_limit)

    stars = []
    psf_x = []
    psf_x_err = []
    psf_y = []
    psf_y_err = []
    for centroid in centroids:
        star = _star_from_centroid(data_array, centroid[0], centroid[1], mean, std, burn_limit, star_std, std_limit)

        if star:
            stars.append([star[0][2], star[0][3], star[0][0], star[0][1], star[0][4], star[0][5],
                          np.sqrt((star[0][2] - x_centre) ** 2 + (star[0][3] - y_centre) ** 2),
                          2 * np.pi * star[0][0] * star[0][4] * star[0][5]])
            psf_x.append(star[0][4])
            psf_x_err.append(np.sqrt(star[1][4][4]))
            psf_y.append(star[0][5])
            psf_y_err.append(np.sqrt(star[1][5][5]))

    if len(stars) > 0:

        psf = (waverage(psf_x, psf_x_err)[0], waverage(psf_y, psf_y_err)[0])

        if order_by_flux:
            stars = sorted(stars, key=lambda x: -x[-1])
        else:
            stars = sorted(stars, key=lambda x: x[-2])

        return stars, psf

    else:
        return None, None


def find_centroids(data_array, x_low, x_upper, y_low, y_upper, mean, std, burn_limit, star_std, std_limit):

    x_upper = int(min(x_upper, len(data_array[0])))
    y_upper = int(min(y_upper, len(data_array)))
    x_low = int(max(0, x_low))
    y_low = int(max(0, y_low))

    data_array = np.full_like(data_array[y_low:y_upper + 1, x_low:x_upper + 1],
                              data_array[y_low:y_upper + 1, x_low:x_upper + 1])

    test = []

    for i in range(-star_std, star_std + 1):
        for j in range(-star_std, star_std + 1):
            rolled = np.roll(np.roll(data_array, i, 0), j, 1)
            test.append(rolled)

    median_test = np.median(test, 0)
    max_test = np.max(test, 0)
    del test
    stars = np.where((data_array < burn_limit) & (data_array > mean + std_limit * std) & (max_test == data_array)
                     & (median_test > mean + std))
    del data_array

    stars = [stars[1] + x_low, stars[0] + y_low]
    stars = np.swapaxes(stars, 0, 1)

    return stars
