from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import warnings


def cartesian_to_polar(x_position, y_position, x_ref_position, y_ref_position):

    x_position, y_position = float(x_position), float(y_position)
    x_ref_position, y_ref_position = float(x_ref_position), float(y_ref_position)

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


def correlation(data_x, data_y):

    data_x = np.array(data_x)
    data_y = np.array(data_y)
    correlation_xy = (np.sum((data_x - np.mean(data_x)) * (data_y - np.mean(data_y))) /
                      ((len(data_x) - 1) * np.std(data_x) * np.std(data_y)))

    return correlation_xy


def waverage(data, uncertainties):

    data = np.array(data)
    uncertainties = np.array(uncertainties)

    weights = 1.0 / (uncertainties ** 2)

    weighted_average = np.sum(data * weights, 0) / np.sum(weights, 0)
    weighted_average_uncertainty = 1.0 / np.sqrt((np.sum(weights, 0)))

    return weighted_average, weighted_average_uncertainty


# decimal points and rounding

def values_to_print(value, error_minus, error_plus):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        try:
            digit1 = abs(int(np.log10(error_minus - int(error_minus)))) + 1
        except OverflowError:
            if error_minus < 1.0:
                digit1 = 3
            else:
                digit1 = 1
        try:
            digit2 = abs(int(np.log10(error_plus - int(error_plus)))) + 1
        except OverflowError:
            if error_plus < 1.0:
                digit2 = 3
            else:
                digit2 = 1
        try:
            done1 = 1 // int(error_minus * (10 ** digit1))
        except ZeroDivisionError:
            done1 = 0
        try:
            done2 = 1 // int(error_plus * (10 ** digit2))
        except ZeroDivisionError:
            done2 = 0

    width = max(digit1 + done1, digit2 + done2)

    print_value = '{0:.{width}f}'.format(round(value, width), width=width)
    print_m_error = '{0:.{width}f}'.format(round(error_minus, width), width=width)
    print_p_error = '{0:.{width}f}'.format(round(error_plus, width), width=width)

    return print_value, print_m_error, print_p_error


def find_decimals(number):

    xx = 30

    while round(number, xx) == round(number, xx - 1):
        xx -= 1

    return xx


def mad(datax):

    datax = np.array(datax, dtype=np.float)
    datax = datax.flatten()

    return np.sqrt(np.median((datax - np.median(datax)) ** 2))
