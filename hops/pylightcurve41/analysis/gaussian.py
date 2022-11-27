
__all__ = ['gaussian', 'fit_gaussian', 'two_d_gaussian', 'fit_two_d_gaussian']


import numpy as np
import warnings

from .curve_fit import curve_fit
from .numerical_integration import sample_function


# gaussian 1D

def gaussian(x_array, model_norm, model_floor, model_mean, model_std):

    return model_floor + (model_norm *
                          np.exp(- 0.5 * (model_mean - x_array) * (model_mean - x_array) / (model_std * model_std)))


def fit_gaussian(datax, datay, positive=False, sampled=False, sampled_precision=3, maxfev=10000, expected_mean=None):

    # TODO option to restrict the space searched

    datax = np.array(datax, dtype=np.float)
    datay = np.array(datay, dtype=np.float)

    # TODO option to point towards the solution

    initial_floor = np.median(datay)
    initial_norm = np.max(datay) - initial_floor
    if expected_mean:
        initial_mean = expected_mean
    else:
        initial_mean = datax[np.argmax(datay)]
    initial_sigma = np.abs(datax[np.argmin(np.abs(datay - np.max(datay) / 2))] - initial_mean)

    initial_values = [initial_norm, initial_floor, initial_mean, initial_sigma]

    if positive:
        def gaussian_to_fit(x_array, model_norm, model_floor, model_mean, model_std):
            return gaussian(x_array, np.abs(model_norm), np.abs(model_floor), model_mean, model_std)
    else:
        gaussian_to_fit = gaussian

    if sampled:
        splits = 0.5 * (datax[1:] + datax[:-1])
        datax1 = np.array([datax[0] - (splits[0] - datax[0])] + list(splits))
        datax2 = np.array(list(splits) + [datax[-1] + (datax[-1] - splits[-1])])

        datax = (datax1, datax2)
        gaussian_to_fit = sample_function(gaussian_to_fit, precision=sampled_precision)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        popt, pcov = curve_fit(gaussian_to_fit, datax, datay, p0=initial_values, maxfev=maxfev)

    popt[3] = np.abs(popt[3])
    if positive:
        popt[0] = np.abs(popt[0])
        popt[1] = np.abs(popt[1])

    return popt, pcov


# gaussian 2D

def two_d_gaussian(x_array, y_array, model_norm, model_floor, model_x_mean, model_y_mean, model_x_sigma, model_y_sigma,
                   model_theta, extend=2):

    model_floor = np.abs(model_floor)
    model_norm = np.abs(model_norm)

    xt_array = x_array - model_x_mean
    yt_array = y_array - model_y_mean

    xt_array, yt_array = (xt_array * np.cos(model_theta) + yt_array * np.sin(model_theta),
                          - xt_array * np.sin(model_theta) + yt_array * np.cos(model_theta))

    return (model_floor + model_norm *
            np.exp(-0.5*((yt_array/model_y_sigma)**2 + (xt_array/model_x_sigma)**extend))
            )


def fit_two_d_gaussian(datax, datay, dataz, errors=None, positive=False, symmetric=False, point_xy=None, sigma=None, floor=None,
                       maxfev=10000, extend=2):

    # TODO option to restrict the space searched

    datax = np.array(datax, dtype=np.float)
    datay = np.array(datay, dtype=np.float)
    dataz = np.array(dataz, dtype=np.float)

    if not floor:
        initial_floor = np.median(dataz)
    else:
        initial_floor = floor
    initial_norm = np.max(dataz) - initial_floor
    if sigma is not None:
        initial_x_sigma = sigma
        initial_y_sigma = sigma
    else:
        initial_x_sigma = 1.0
        initial_y_sigma = 1.0

    if not point_xy:
        test_data_x = np.sum(dataz, 0)
        test_floor = np.median(test_data_x)
        test_norm = np.max(test_data_x) - test_floor

        initial_x_mean_arg = int(len(datax[0]) / 2)
        model = gaussian(datax[0], test_norm, test_floor, datax[0][initial_x_mean_arg], initial_x_sigma)
        dx = np.argmax(np.convolve(test_data_x, model)) - np.argmax(np.convolve(model, model))
        initial_x_mean = datax[0][initial_x_mean_arg + dx]

        test_data_y = np.sum(dataz, 1)
        test_floor = np.median(test_data_y)
        test_norm = np.max(test_data_y) - test_floor

        initial_y_mean_arg = int(len(datay[:, 0]) / 2)
        model = gaussian(datay[:, 0], test_norm, test_floor, datay[:, 0][initial_y_mean_arg], initial_y_sigma)
        dy = np.argmax(np.convolve(np.sum(dataz, 1), model)) - np.argmax(np.convolve(model, model))
        initial_y_mean = datay[:, 0][initial_y_mean_arg + dy]

    else:
        initial_x_mean, initial_y_mean = np.int_(point_xy)

    if sigma is None:
        initial_x_sigma = np.abs(datax[0][np.argmin(np.abs(np.sum(dataz, 0) - np.max(np.sum(dataz, 1)) / 2))] -
                                 initial_x_mean)

        initial_y_sigma = np.abs(datay[:, 0][np.argmin(np.abs(np.sum(dataz, 1) - np.max(np.sum(dataz, 0)) / 2))] -
                                 initial_y_mean)

    if positive and symmetric:

        def gaussian_to_fit(xy_array, model_norm, model_floor, model_x_mean, model_y_mean, model_sigma):
            return two_d_gaussian(xy_array[0], xy_array[1], np.abs(model_norm), np.abs(model_floor), model_x_mean,
                                  model_y_mean, model_sigma, model_sigma, 0, extend=extend)

        initial_values = [initial_norm, initial_floor, initial_x_mean, initial_y_mean, initial_x_sigma]

    elif positive:

        def gaussian_to_fit(xy_array, model_norm, model_floor, model_x_mean, model_y_mean, model_x_sigma, model_y_sigma,
                            model_theta):
            return two_d_gaussian(xy_array[0], xy_array[1], np.abs(model_norm), np.abs(model_floor), model_x_mean,
                                  model_y_mean, model_x_sigma, model_y_sigma, model_theta, extend=extend)

        initial_values = [initial_norm, initial_floor, initial_x_mean, initial_y_mean,
                          initial_x_sigma, initial_y_sigma, 0]

    elif symmetric:

        def gaussian_to_fit(xy_array, model_norm, model_floor, model_x_mean, model_y_mean, model_sigma):
            return two_d_gaussian(xy_array[0], xy_array[1], model_norm, model_floor, model_x_mean,
                                  model_y_mean, model_sigma, model_sigma, 0, extend=extend)

        initial_values = [initial_norm, initial_floor, initial_x_mean, initial_y_mean, initial_x_sigma]

    else:

        def gaussian_to_fit(xy_array, model_norm, model_floor, model_x_mean, model_y_mean, model_x_sigma, model_y_sigma,
                            model_theta):
            return two_d_gaussian(xy_array[0], xy_array[1], model_norm, model_floor, model_x_mean, model_y_mean,
                                  model_x_sigma, model_y_sigma, model_theta, extend=extend)

        initial_values = [initial_norm, initial_floor, initial_x_mean, initial_y_mean,
                          initial_x_sigma, initial_y_sigma, 0]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if errors is not None:
            popt, pcov = curve_fit(gaussian_to_fit, (datax.flatten(), datay.flatten()), dataz.flatten(), p0=initial_values, sigma=errors.flatten(), maxfev=maxfev)
        else:
            popt, pcov = curve_fit(gaussian_to_fit, (datax.flatten(), datay.flatten()), dataz.flatten(), p0=initial_values, maxfev=maxfev)

    if positive:
        popt[0] = np.abs(popt[0])
        popt[1] = np.abs(popt[1])

    popt[4] = np.abs(popt[4])
    if not symmetric:
        popt[5] = np.abs(popt[5])

        if popt[5] > popt[4]:
            popt[4], popt[5] = popt[5], popt[4]
            pcov[4][4], pcov[5][5] = pcov[5][5], pcov[4][4]
            popt[6] -= np.pi/2

    return popt, pcov

