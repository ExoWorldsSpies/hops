__all__ = ['find_all_stars', 'find_single_star', 'find_centroids', 'fast_psf_find', 'pixel_to_aperture_overlap']

import sys
import numpy as np
import warnings

import hops.pylightcurve41 as plc


def _star_from_centroid(data_array, centroid_x, centroid_y, mean, std, star_std, snr=3,
                        force_circles=False):

    star = None
    star_std = int(round(star_std))
    try:
        search_window = int(round(6 * star_std))
        y_min = int(max(int(centroid_y) - search_window, 0))
        y_max = int(min(int(centroid_y) + search_window, len(data_array) - 1))
        x_min = int(max(int(centroid_x) - search_window, 0))
        x_max = int(min(int(centroid_x) + search_window, len(data_array[0]) - 1))

        datax, datay = np.meshgrid(np.arange(x_min, x_max + 1) + 0.5,
                                   np.arange(y_min, y_max + 1) + 0.5)

        dataz = data_array[y_min: y_max + 1, x_min: x_max + 1]

        popt, pcov = plc.fit_two_d_gaussian(datax, datay, dataz, point_xy=(centroid_x, centroid_y),
                                        sigma=star_std, positive=True, floor=mean, maxfev=1000,
                                        symmetric=force_circles)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if popt[0] > snr * std:
                if np.nan not in [np.sqrt(abs(pcov[ff][ff])) for ff in range(len(pcov) - 1)]:
                    if np.inf not in [np.sqrt(abs(pcov[ff][ff])) for ff in range(len(pcov) - 1)]:
                        if 0 not in [np.sqrt(abs(pcov[ff][ff])) for ff in range(len(pcov) - 1)]:
                            star = (popt, pcov)

    except:
        pass

    return star


def find_single_star(data_array, predicted_x, predicted_y, window=5,
                     mean=None, std=None, burn_limit=65000, star_std=2,
                     std_limit=3.0):

    star = None
    star_std = int(round(star_std))

    if 0 < predicted_x < len(data_array[0]) and 0 < predicted_y < len(data_array):
        if mean is None or std is None:
            try:
                fit_mean, fit_std = plc.one_d_distribution(data_array, gaussian_fit=True)[2:4]
            except:
                fit_mean = np.mean(data_array)
                fit_std = np.std(data_array)

            if not mean:
                mean = fit_mean

            if not std:
                std = fit_std

        centroids = find_centroids(data_array, predicted_x - window * star_std, predicted_x + window * star_std,
                                   predicted_y - 5 * star_std, predicted_y + 5 * star_std, mean, std, burn_limit, star_std,
                                   std_limit)

        centroids = sorted(centroids, key=lambda x: np.sqrt((x[1] - predicted_x) ** 2 + (x[2] - predicted_y) ** 2))

        for centroid in centroids:
            star = _star_from_centroid(data_array, centroid[1], centroid[2], mean, std, star_std, std_limit)
            if star:
                star = [star[0][2], star[0][3], star[0][0], star[0][1], star[0][4], star[0][5], centroid[0], centroid[1]]
                break

    return star


def find_all_stars(data_array, x_low=0, x_upper=None, y_low=0, y_upper=None, x_centre=None, y_centre=None,
                   mean=None, std=None, burn_limit=65000, star_std=None, std_limit=5.0,
                   force_circles=False, psf_variation_allowed=0.5,
                   order_by_flux=False, order_by_distance_and_flux=False,
                   progressbar=None, progress_window=None, verbose=False):

    star_std = int(round(star_std))

    if mean is None or std is None:
        try:
            fit_mean, fit_std = plc.one_d_distribution(data_array, samples=10000, gaussian_fit=True)[2:4]
        except:
            fit_mean = np.mean(data_array)
            fit_std = np.std(data_array)

        if not mean:
            mean = fit_mean

        if not std:
            std = fit_std

    if not star_std:
        star_std = fast_psf_find(data_array, mean, std, burn_limit)

    if x_upper is None:
        x_upper = data_array.shape[1]

    if y_upper is None:
        y_upper = data_array.shape[0]

    if x_centre is None:
        x_centre = data_array.shape[1] / 2

    if y_centre is None:
        y_centre = data_array.shape[0] / 2

    if verbose:
        print('\nAnalysing frame...')

    centroids = find_centroids(data_array, x_low, x_upper, y_low, y_upper, mean, std, burn_limit, star_std, std_limit)

    if progress_window:
        if progress_window.exit:
            return None, None

    stars = []
    psf_x = []
    psf_x_err = []
    psf_y = []
    psf_y_err = []
    if verbose:
        print('Verifying stars...')

    for num, centroid in enumerate(centroids):

        star = _star_from_centroid(data_array, centroid[1], centroid[2], mean, std, star_std, force_circles=force_circles)

        if star:

            if force_circles:
                star = [star[0], star[1]]
                star[0] = list(star[0])
                star[0].append(star[0][4])

            stars.append([star[0][2], star[0][3], star[0][0], star[0][1], star[0][4], star[0][5],
                          np.sqrt((star[0][2] - x_centre) ** 2 + (star[0][3] - y_centre) ** 2),
                          2 * np.pi * star[0][0] * star[0][4] * star[0][5]])

            if force_circles:
                psf_x.append(star[0][4])
                psf_x_err.append(np.sqrt(abs(star[1][4][4])))
                psf_y.append(star[0][4])
                psf_y_err.append(np.sqrt(abs(star[1][4][4])))

            else:

                if star[0][4] > star[0][5]:
                    psf_x.append(star[0][4])
                    psf_x_err.append(np.sqrt(abs(star[1][4][4])))
                    psf_y.append(star[0][5])
                    psf_y_err.append(np.sqrt(abs(star[1][5][5])))
                else:
                    psf_y.append(star[0][4])
                    psf_y_err.append(np.sqrt(abs(star[1][4][4])))
                    psf_x.append(star[0][5])
                    psf_x_err.append(np.sqrt(abs(star[1][5][5])))

        if verbose:
            sys.stdout.write('\r\033[K')
            sys.stdout.write('{0}/{1}'.format(num + 1, len(centroids)))
            sys.stdout.flush()

    if verbose:
        print('')

    if len(stars) > 0:

        psf = (plc.waverage(np.array(psf_x), np.array(psf_x_err))[0], plc.waverage(np.array(psf_y), np.array(psf_y_err))[0])

        not_trails = np.where((psf_x < (1.0 + psf_variation_allowed) * psf[0]) * (psf_x > (1 - psf_variation_allowed) * psf[0]))

        stars = np.array(stars)[not_trails]

        data_array_len = 1.0 * len(data_array)

        if order_by_flux:
            if order_by_flux == True:
                stars = sorted(stars, key=lambda x: -x[-1])
            else:
                stars = sorted(stars, key=lambda x: abs(x[-1] - order_by_flux))
        elif order_by_distance_and_flux:
            if order_by_distance_and_flux == True:
                stars = sorted(stars, key=lambda x: -x[-1] / (x[-2] ** 3))
            else:
                stars = sorted(stars, key=lambda x: (x[-2]/data_array_len) + 100 * abs((x[-1] - order_by_distance_and_flux)/order_by_distance_and_flux))
        else:
            stars = sorted(stars, key=lambda x: x[-2])

        return stars, psf

    else:
        return None, None


def find_centroids(data_array, x_low, x_upper, y_low, y_upper, mean, std, burn_limit, psf, snr):

    psf = int(round(psf))
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

    median_test = np.median(test, 0)
    max_test = np.max(test, 0)
    del test

    data_array_test = data_array[bright]

    stars = np.where((max_test < burn_limit) & (max_test == data_array_test) & (median_test > mean + 2 * std))[0]
    stars = np.swapaxes([data_array_test[stars], bright[1][stars] + x_low, bright[0][stars] + y_low], 0, 1)

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

    return stars


def fast_psf_find(data_array, mean, std, burn_limit, sample=10):

    psf = 2

    stars = []
    snr = 53

    while len(stars) < sample * 2 and snr >= 3:

        limit = mean + snr * std

        bright = np.where(data_array[psf:-psf, psf:-psf] > limit)
        bright = (bright[0] + psf, bright[1] + psf)

        test = []
        for i in range(-psf, psf + 1):
            for j in range(-psf, psf + 1):
                test.append(data_array[bright[0] + i, bright[1]+j])

        median_test = np.median(test, 0)
        max_test = np.max(test, 0)
        data_array_test = data_array[bright]

        del test
        stars = np.where((max_test < burn_limit) & (max_test == data_array_test) & (median_test > mean + 2 * std))[0]
        stars = np.swapaxes([data_array_test[stars], bright[1][stars], bright[0][stars]], 0, 1)

        stars = sorted(stars, key=lambda x: -x[0])

        snr -= 10

    psf_x = []
    psf_x_err = []
    psf_y = []
    psf_y_err = []

    for centroid in stars:
        star = _star_from_centroid(data_array, centroid[1], centroid[2], mean, std, psf)

        if star:
            if star[0][4] > star[0][5]:
                psf_x.append(star[0][4])
                psf_x_err.append(np.sqrt(star[1][4][4]))
                psf_y.append(star[0][5])
                psf_y_err.append(np.sqrt(star[1][5][5]))
            else:
                psf_y.append(star[0][4])
                psf_y_err.append(np.sqrt(star[1][4][4]))
                psf_x.append(star[0][5])
                psf_x_err.append(np.sqrt(star[1][5][5]))

        if len(psf_x) == sample:
            break

    psf = (plc.waverage(np.array(psf_x), np.array(psf_x_err))[0], plc.waverage(np.array(psf_y), np.array(psf_y_err))[0])

    return max(psf)


# circle to square

def solve(c1_x, c1_y, c2_x, c2_y):
    if c1_x > c2_x:
        c1_x, c1_y, c2_x, c2_y = c2_x, c2_y, c1_x, c1_y

    a = (c2_y - c1_y)/(c2_x - c1_x)
    b = (c1_y - c1_x * a)
    D = (2 * a * b) ** 2 - 4 * (1 + a ** 2) * (b ** 2 - 1)
    cuts = []
    if D == 0:
        sol_x = - (a * b) / (1 + a ** 2)
        if (sol_x > c1_x) * (sol_x < c2_x):
            cuts.append(sol_x)
    elif D >0 :
        sol_x = (- (2 * a * b) + np.sqrt(D)) / (2 * (1 + a ** 2))
        if (sol_x > c1_x) * (sol_x < c2_x):
            cuts.append(sol_x)
        sol_xx = (- (2 * a * b) - np.sqrt(D)) / (2 * (1 + a ** 2))
        if (sol_xx > c1_x) * (sol_xx < c2_x):
            cuts.append(sol_xx)

    cuts = [(x, a * x + b) for x in cuts]

    if len(cuts) == 2:
        if cuts[0][1] < cuts[1][1]:
            cuts = [cuts[1], cuts[0]]

    return cuts


def find_corners_and_cuts(a, d, theta):

    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)

    aa = np.sqrt(0.5) * a

    c1_x = d - aa * cos_theta
    c1_y = aa * sin_theta

    c2_x = d + aa * sin_theta
    c2_y = aa * cos_theta

    c3_x = d + aa * cos_theta
    c3_y = - aa * sin_theta

    c4_x = d - aa * sin_theta
    c4_y = - aa * cos_theta

    line_1 = solve(c1_x, c1_y, c2_x, c2_y)
    line_2 = solve(c2_x, c2_y, c3_x, c3_y)
    line_3 = solve(c3_x, c3_y, c4_x, c4_y)
    line_4 = solve(c4_x, c4_y, c1_x, c1_y)

    return ([(c1_x, c1_y), (c2_x, c2_y), (c3_x, c3_y), (c4_x, c4_y)], [line_1, line_2, line_3, line_4])


def get_distance(p1, p2):
    return np.sqrt((p1[0] - p2[0]) **2 + (p1[1] - p2[1]) ** 2)

def get_circular_sector_area(p1, p2):
    return np.arcsin(0.5 * get_distance(p1, p2))

def get_triangle_area(p1, p2, p3):
    return np.abs(0.5 * (p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1])))

def get_circle_to_chord_area(p1, p2):
    return get_circular_sector_area(p1, p2) - get_triangle_area((0, 0), p1, p2)


def get_area(a, d, theta):
    if d <= 1 - np.sqrt(2)*a/2:
        return 1
    elif d >= 1 + np.sqrt(2)*a/2:
        return 0

    elif theta == np.pi / 4:
        if (d + a/2.0) ** 2 + (a/2.0)**2 <= 1:
            return 1
        elif d + a/2.0 < 1:
            cut_1 = (np.sqrt(1 - (a/2.0)**2), a/2.0)
            cut_2 = (d + a/2.0, np.sqrt(1 - (d + a/2.0)**2))
            area = (
                    get_triangle_area((d - a/2.0, a/2.0), (d - a/2.0, 0.0), cut_1) +
                    get_triangle_area(cut_1, (d - a/2.0, 0.0), cut_2) +
                    get_triangle_area(cut_2, (d - a/2.0, 0.0), (d + a/2.0, 0)) +
                    get_circle_to_chord_area(cut_1, cut_2)
            )
            return 2 * area / a ** 2
        elif (d - 0.5 * a) ** 2 + (0.5 * a)**2 < 1:
            cut_1 = (np.sqrt(1 - (a/2.0)**2), a/2.0)
            area = (
                    get_triangle_area((d - a/2.0, a/2.0), (d - a/2.0, 0.0), cut_1) +
                    get_triangle_area(cut_1, (d - a/2.0, 0.0), (1, 0)) +
                    get_circle_to_chord_area(cut_1, (1, 0))
            )
            return 2 * area / a ** 2
        elif d - a/2.0 < 1:
            cut_1 = (d - a/2.0, np.sqrt(1 - (d - a/2.0)**2))
            area = (
                get_circle_to_chord_area(cut_1, (1, 0))
            )
            return 2 * area / a ** 2
        else:
            return 0

    else:

        corners, cuts = find_corners_and_cuts(a, d, theta)
        corner_1, corner_2, corner_3, corner_4 = corners
        cuts_1, cuts_2, cuts_3, cuts_4 = cuts
        cuts_numbers = [len(cut) for cut in cuts]

        if np.sum(cuts_numbers) == 0:
            if np.sqrt(corner_1[0] ** 2 + corner_1[1] ** 2) < 1:
                return 1
            else:
                return 0

        elif (corner_2[0] ** 2 + corner_2[1] ** 2) == 1:
            area = (
                    get_triangle_area(corner_1, corner_4, corner_2) +
                    get_triangle_area(corner_4, corner_2, cuts_3[0]) +
                    get_circle_to_chord_area(corner_2, cuts_3[0])
            )
            return area / a ** 2

        elif (corner_4[0] ** 2 + corner_4[1] ** 2) == 1:
            area = (
                    get_triangle_area(corner_1, cuts_1[0], corner_4) +
                    get_circle_to_chord_area(cuts_1[0], corner_4)
            )
            return area / a ** 2

        elif cuts_numbers == [1, 0, 0, 1]:
            area = (
                    get_triangle_area(corner_1, cuts_1[0], cuts_4[0]) +
                    get_circle_to_chord_area(cuts_1[0], cuts_4[0])
            )
            return area / a ** 2

        elif cuts_numbers == [1, 0, 1, 0]:
            area = (
                    get_triangle_area(corner_1, corner_4, cuts_1[0]) +
                    get_triangle_area(corner_4, cuts_1[0], cuts_3[0]) +
                    get_circle_to_chord_area(cuts_1[0], cuts_3[0])
            )
            return area / a ** 2

        elif cuts_numbers == [0, 1, 1, 0]:
            area = (
                    get_triangle_area(corner_1, corner_2, corner_4) +
                    get_triangle_area(corner_2, corner_4, cuts_2[0]) +
                    get_triangle_area(corner_4, cuts_2[0], cuts_3[0]) +
                    get_circle_to_chord_area(cuts_2[0], cuts_3[0])
            )
            return area / a ** 2

        elif cuts_numbers == [1, 2, 1, 0]:
            area = (
                    get_triangle_area(corner_1, corner_4, cuts_1[0]) +
                    get_triangle_area(corner_4, cuts_1[0], cuts_3[0]) +
                    get_triangle_area(cuts_1[0], cuts_3[0], cuts_2[0]) +
                    get_triangle_area(cuts_3[0], cuts_2[0], cuts_2[1]) +
                    get_circle_to_chord_area(cuts_1[0], cuts_2[0]) +
                    get_circle_to_chord_area(cuts_2[1], cuts_3[0])
            )
            return area / a ** 2

        elif cuts_numbers == [0, 0, 0, 2]:
            area = (
                get_circle_to_chord_area(cuts_4[0], cuts_4[1])
            )
            return area / a ** 2

        else:
            print('pending')


def pixel_to_aperture_overlap(x_pix, y_pix, x_ap, y_ap, ap):

    x1, y1 = x_pix - 0.5, y_pix + 0.5
    x2, y2 = x_pix + 0.5, y_pix + 0.5
    x3, y3 = x_pix + 0.5, y_pix - 0.5
    x4, y4 = x_pix - 0.5, y_pix - 0.5

    rc = np.sqrt( (x_pix - x_ap) ** 2 + (y_pix - y_ap) ** 2)
    r1 = np.sqrt( (x1 - x_ap) ** 2 + (y1 - y_ap) ** 2)
    r2 = np.sqrt( (x2 - x_ap) ** 2 + (y2 - y_ap) ** 2)
    if r2 < r1:
        r1 = r2
    r3 = np.sqrt( (x3 - x_ap) ** 2 + (y3 - y_ap) ** 2)
    if r3 < r1:
        r1 = r3
    r4 = np.sqrt( (x4 - x_ap) ** 2 + (y4 - y_ap) ** 2)
    if r4 < r1:
        r1 = r4

    theta = np.arccos((r1 ** 2 - rc ** 2 - 0.5) / (-np.sqrt(2) * rc))
    a = 1 / ap
    d = np.sqrt((x_pix - x_ap) ** 2 + (y_pix - y_ap) ** 2) / ap

    return get_area(a, d, theta)
