__all__ = ['waverage', 'chi_squared', 'reduced_chi_squared', 'mad', 'mean_std_from_median_mad',
           'cartesian_to_polar', 'correlation', 'drift_rotate']

import numpy as np


def waverage(data, uncertainties, axis=None):

    weights = 1.0 / (uncertainties ** 2)

    weighted_average = np.sum(data * weights, axis=axis) / np.sum(weights, axis=axis)
    weighted_average_uncertainty = 1.0 / np.sqrt((np.sum(weights, axis=axis)))

    return weighted_average, weighted_average_uncertainty


def chi_squared(data, model, uncertainties):
    return np.sum(((data - model) ** 2) / (uncertainties ** 2))


def reduced_chi_squared(data, model, uncertainties, free_parameters):
    return chi_squared(data, model, uncertainties) / len(data - free_parameters)


def mad(datax):

    datax = 1.0 * datax.flatten()

    return np.sqrt(np.median((datax - np.median(datax)) ** 2))


def mean_std_from_median_mad(datax, samples=None):

    datax = 1.0 * datax.flatten()

    if samples:
        datax = datax[::int(len(datax)/samples + 1)]

    median = np.median(datax)
    dx = datax - median
    mad = np.sqrt(np.median(dx * dx))

    return median, mad * 0.6745


def cartesian_to_polar(x_position, y_position, x_ref_position, y_ref_position):

    radius = np.sqrt((x_position - x_ref_position) ** 2 + (y_position - y_ref_position) ** 2)

    if (x_position - x_ref_position) > 0:
        if (y_position - y_ref_position) >= 0:
            angle = np.arctan((y_position - y_ref_position) / (x_position - x_ref_position))
        else:
            angle = 2.0 * np.pi + np.arctan((y_position - y_ref_position) / (x_position - x_ref_position))
    elif (x_position - x_ref_position) < 0:
        angle = np.arctan((y_position - y_ref_position) / (x_position - x_ref_position)) + np.pi
    else:
        if (y_position - y_ref_position) >= 0:
            angle = np.pi / 2
        else:
            angle = 3 * np.pi / 2

    return radius, angle


def drift_rotate(x_position, y_position, x_ref_position, y_ref_position, dx, dy, dtheta):

    rr, tt = cartesian_to_polar(x_position, y_position, x_ref_position, y_ref_position)

    tt = tt + dtheta

    xx = x_ref_position + rr * np.cos(tt) + dx
    yy = y_ref_position + rr * np.sin(tt) + dy

    return xx, yy


def correlation(data_x, data_y):

    data_x = np.array(data_x)
    data_y = np.array(data_y)
    correlation_xy = (np.sum((data_x - np.mean(data_x)) * (data_y - np.mean(data_y))) /
                      ((len(data_x) - 1) * np.std(data_x) * np.std(data_y)))

    return correlation_xy