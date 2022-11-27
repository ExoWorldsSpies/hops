
__all__ = ['planet_orbit', 'planet_star_projected_distance', 'planet_phase',
           'transit', 'transit_integrated', 'transit_depth', 'transit_duration',
           'eclipse', 'eclipse_integrated', 'eclipse_depth', 'eclipse_duration', 'eclipse_mid_time',
           'fp_over_fs','transit_t12', 'exotethys']


import os
import glob
import warnings
import pickle
import numpy as np

from ..errors import *
from ..__databases__ import plc_data
from ..analysis.numerical_integration import gauss_numerical_integration
from ..analysis.curve_fit import curve_fit
from ..processes.files import open_dict


# orbit


def planet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array, ww=0):

    inclination = inclination * np.pi / 180.0
    periastron = periastron * np.pi / 180.0
    ww = ww * np.pi / 180.0

    if eccentricity == 0 and ww == 0:
        vv = 2 * np.pi * (time_array - mid_time) / period
        bb = sma_over_rs * np.cos(vv)
        return [bb * np.sin(inclination), sma_over_rs * np.sin(vv), - bb * np.cos(inclination)]

    if periastron < np.pi / 2:
        aa = 1.0 * np.pi / 2 - periastron
    else:
        aa = 5.0 * np.pi / 2 - periastron
    bb = 2 * np.arctan(np.sqrt((1 - eccentricity) / (1 + eccentricity)) * np.tan(aa / 2))
    if bb < 0:
        bb += 2 * np.pi
    mid_time = float(mid_time) - (period / 2.0 / np.pi) * (bb - eccentricity * np.sin(bb))
    m = (time_array - mid_time - np.int_((time_array - mid_time) / period) * period) * 2.0 * np.pi / period
    u0 = m
    stop = False
    u1 = 0
    for ii in range(10000):  # setting a limit of 1k iterations - arbitrary limit
        u1 = u0 - (u0 - eccentricity * np.sin(u0) - m) / (1 - eccentricity * np.cos(u0))
        stop = (np.abs(u1 - u0) < 10 ** (-7)).all()
        if stop:
            break
        else:
            u0 = u1
    if not stop:
        raise RuntimeError('Failed to find a solution in 10000 loops')

    vv = 2 * np.arctan(np.sqrt((1 + eccentricity) / (1 - eccentricity)) * np.tan(u1 / 2))
    #
    rr = sma_over_rs * (1 - (eccentricity ** 2)) / (np.ones_like(vv) + eccentricity * np.cos(vv))
    aa = np.cos(vv + periastron)
    bb = np.sin(vv + periastron)
    x = rr * bb * np.sin(inclination)
    y = rr * (-aa * np.cos(ww) + bb * np.sin(ww) * np.cos(inclination))
    z = rr * (-aa * np.sin(ww) - bb * np.cos(ww) * np.cos(inclination))

    return [x, y, z]


def planet_star_projected_distance(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array):

    position_vector = planet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array)

    return np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2])


def planet_phase(period, mid_time, time_array):
    return (time_array - mid_time)/period


# flux drop


def integral_r_claret(limb_darkening_coefficients, r):
    a1, a2, a3, a4 = limb_darkening_coefficients
    mu44 = 1.0 - r * r
    mu24 = np.sqrt(mu44)
    mu14 = np.sqrt(mu24)
    return - (2.0 * (1.0 - a1 - a2 - a3 - a4) / 4) * mu44 \
           - (2.0 * a1 / 5) * mu44 * mu14 \
           - (2.0 * a2 / 6) * mu44 * mu24 \
           - (2.0 * a3 / 7) * mu44 * mu24 * mu14 \
           - (2.0 * a4 / 8) * mu44 * mu44


def num_claret(r, limb_darkening_coefficients, rprs, z):
    a1, a2, a3, a4 = limb_darkening_coefficients
    rsq = r * r
    mu44 = 1.0 - rsq
    mu24 = np.sqrt(mu44)
    mu14 = np.sqrt(mu24)
    return ((1.0 - a1 - a2 - a3 - a4) + a1 * mu14 + a2 * mu24 + a3 * mu24 * mu14 + a4 * mu44) \
        * r * np.arccos(np.minimum((-rprs ** 2 + z * z + rsq) / (2.0 * z * r), 1.0))


def integral_r_f_claret(limb_darkening_coefficients, rprs, z, r1, r2, precision=3):
    return gauss_numerical_integration(num_claret, r1, r2, precision, limb_darkening_coefficients, rprs, z)


# integral definitions for zero method


def integral_r_zero(limb_darkening_coefficients, r):
    musq = 1 - r * r
    return (-1.0 / 6) * musq * 3.0


def num_zero(r, limb_darkening_coefficients, rprs, z):
    rsq = r * r
    return r * np.arccos(np.minimum((-rprs ** 2 + z * z + rsq) / (2.0 * z * r), 1.0))


def integral_r_f_zero(limb_darkening_coefficients, rprs, z, r1, r2, precision=3):
    return gauss_numerical_integration(num_zero, r1, r2, precision, limb_darkening_coefficients, rprs, z)


# integral definitions for linear method


def integral_r_linear(limb_darkening_coefficients, r):
    a1 = limb_darkening_coefficients[0]
    musq = 1 - r * r
    return (-1.0 / 6) * musq * (3.0 + a1 * (-3.0 + 2.0 * np.sqrt(musq)))


def num_linear(r, limb_darkening_coefficients, rprs, z):
    a1 = limb_darkening_coefficients[0]
    rsq = r * r
    return (1.0 - a1 * (1.0 - np.sqrt(1.0 - rsq))) \
        * r * np.arccos(np.minimum((-rprs ** 2 + z * z + rsq) / (2.0 * z * r), 1.0))


def integral_r_f_linear(limb_darkening_coefficients, rprs, z, r1, r2, precision=3):
    return gauss_numerical_integration(num_linear, r1, r2, precision, limb_darkening_coefficients, rprs, z)


# integral definitions for quadratic method


def integral_r_quad(limb_darkening_coefficients, r):
    a1, a2 = limb_darkening_coefficients[:2]
    musq = 1 - r * r
    mu = np.sqrt(musq)
    return (1.0 / 12) * (-4.0 * (a1 + 2.0 * a2) * mu * musq + 6.0 * (-1 + a1 + a2) * musq + 3.0 * a2 * musq * musq)


def num_quad(r, limb_darkening_coefficients, rprs, z):
    a1, a2 = limb_darkening_coefficients[:2]
    rsq = r * r
    cc = 1.0 - np.sqrt(1.0 - rsq)
    return (1.0 - a1 * cc - a2 * cc * cc) \
        * r * np.arccos(np.minimum((-rprs ** 2 + z * z + rsq) / (2.0 * z * r), 1.0))


def integral_r_f_quad(limb_darkening_coefficients, rprs, z, r1, r2, precision=3):
    return gauss_numerical_integration(num_quad, r1, r2, precision, limb_darkening_coefficients, rprs, z)


# integral definitions for square root method


def integral_r_sqrt(limb_darkening_coefficients, r):
    a1, a2 = limb_darkening_coefficients[:2]
    musq = 1 - r * r
    mu = np.sqrt(musq)
    return ((-2.0 / 5) * a2 * np.sqrt(mu) - (1.0 / 3) * a1 * mu + (1.0 / 2) * (-1 + a1 + a2)) * musq


def num_sqrt(r, limb_darkening_coefficients, rprs, z):
    a1, a2 = limb_darkening_coefficients[:2]
    rsq = r * r
    mu = np.sqrt(1.0 - rsq)
    return (1.0 - a1 * (1 - mu) - a2 * (1.0 - np.sqrt(mu))) \
        * r * np.arccos(np.minimum((-rprs ** 2 + z * z + rsq) / (2.0 * z * r), 1.0))


def integral_r_f_sqrt(limb_darkening_coefficients, rprs, z, r1, r2, precision=3):
    return gauss_numerical_integration(num_sqrt, r1, r2, precision, limb_darkening_coefficients, rprs, z)


# dictionaries containing the different methods,
# if you define a new method, include the functions in the dictionary as well

integral_r = {
    'claret': integral_r_claret,
    'linear': integral_r_linear,
    'quad': integral_r_quad,
    'sqrt': integral_r_sqrt,
    'zero': integral_r_zero
}

integral_r_f = {
    'claret': integral_r_f_claret,
    'linear': integral_r_f_linear,
    'quad': integral_r_f_quad,
    'sqrt': integral_r_f_sqrt,
    'zero': integral_r_f_zero,
}


def integral_centred(method, limb_darkening_coefficients, rprs, ww1, ww2):
    return (integral_r[method](limb_darkening_coefficients, rprs)
            - integral_r[method](limb_darkening_coefficients, 0.0)) * np.abs(ww2 - ww1)


def integral_plus_core(method, limb_darkening_coefficients, rprs, z, ww1, ww2, precision=3):
    if len(z) == 0:
        return z
    rr1 = z * np.cos(ww1) + np.sqrt(np.maximum(rprs ** 2 - (z * np.sin(ww1)) ** 2, 0))
    rr1 = np.clip(rr1, 0, 1)
    rr2 = z * np.cos(ww2) + np.sqrt(np.maximum(rprs ** 2 - (z * np.sin(ww2)) ** 2, 0))
    rr2 = np.clip(rr2, 0, 1)
    w1 = np.minimum(ww1, ww2)
    r1 = np.minimum(rr1, rr2)
    w2 = np.maximum(ww1, ww2)
    r2 = np.maximum(rr1, rr2)
    parta = integral_r[method](limb_darkening_coefficients, 0.0) * (w1 - w2)
    partb = integral_r[method](limb_darkening_coefficients, r1) * w2
    partc = integral_r[method](limb_darkening_coefficients, r2) * (-w1)
    partd = integral_r_f[method](limb_darkening_coefficients, rprs, z, r1, r2, precision=precision)
    return parta + partb + partc + partd


def integral_minus_core(method, limb_darkening_coefficients, rprs, z, ww1, ww2, precision=3):
    if len(z) == 0:
        return z
    rr1 = z * np.cos(ww1) - np.sqrt(np.maximum(rprs ** 2 - (z * np.sin(ww1)) ** 2, 0))
    rr1 = np.clip(rr1, 0, 1)
    rr2 = z * np.cos(ww2) - np.sqrt(np.maximum(rprs ** 2 - (z * np.sin(ww2)) ** 2, 0))
    rr2 = np.clip(rr2, 0, 1)
    w1 = np.minimum(ww1, ww2)
    r1 = np.minimum(rr1, rr2)
    w2 = np.maximum(ww1, ww2)
    r2 = np.maximum(rr1, rr2)
    parta = integral_r[method](limb_darkening_coefficients, 0.0) * (w1 - w2)
    partb = integral_r[method](limb_darkening_coefficients, r1) * (-w1)
    partc = integral_r[method](limb_darkening_coefficients, r2) * w2
    partd = integral_r_f[method](limb_darkening_coefficients, rprs, z, r1, r2, precision=precision)
    return parta + partb + partc - partd


def transit_flux_drop(limb_darkening_coefficients, rp_over_rs, z_over_rs, method='claret', precision=3):

    z_over_rs = np.where(z_over_rs < 0, 1.0 + 100.0 * rp_over_rs, z_over_rs)
    z_over_rs = np.maximum(z_over_rs, 10**(-10))

    # cases
    zsq = z_over_rs * z_over_rs
    sum_z_rprs = z_over_rs + rp_over_rs
    dif_z_rprs = rp_over_rs - z_over_rs
    sqr_dif_z_rprs = zsq - rp_over_rs ** 2
    case0 = np.where((z_over_rs == 0) & (rp_over_rs <= 1))
    case1 = np.where((z_over_rs < rp_over_rs) & (sum_z_rprs <= 1))
    casea = np.where((z_over_rs < rp_over_rs) & (sum_z_rprs > 1) & (dif_z_rprs < 1))
    caseb = np.where((z_over_rs < rp_over_rs) & (sum_z_rprs > 1) & (dif_z_rprs > 1))
    case2 = np.where((z_over_rs == rp_over_rs) & (sum_z_rprs <= 1))
    casec = np.where((z_over_rs == rp_over_rs) & (sum_z_rprs > 1))
    case3 = np.where((z_over_rs > rp_over_rs) & (sum_z_rprs < 1))
    case4 = np.where((z_over_rs > rp_over_rs) & (sum_z_rprs == 1))
    case5 = np.where((z_over_rs > rp_over_rs) & (sum_z_rprs > 1) & (sqr_dif_z_rprs < 1))
    case6 = np.where((z_over_rs > rp_over_rs) & (sum_z_rprs > 1) & (sqr_dif_z_rprs == 1))
    case7 = np.where((z_over_rs > rp_over_rs) & (sum_z_rprs > 1) & (sqr_dif_z_rprs > 1) & (-1 < dif_z_rprs))
    plus_case = np.concatenate((case1[0], case2[0], case3[0], case4[0], case5[0], casea[0], casec[0]))
    minus_case = np.concatenate((case3[0], case4[0], case5[0], case6[0], case7[0]))
    star_case = np.concatenate((case5[0], case6[0], case7[0], casea[0], casec[0]))

    # cross points
    ph = np.arccos(np.clip((1.0 - rp_over_rs ** 2 + zsq) / (2.0 * z_over_rs), -1, 1))
    theta_1 = np.zeros(len(z_over_rs))
    ph_case = np.concatenate((case5[0], casea[0], casec[0]))
    theta_1[ph_case] = ph[ph_case]
    theta_2 = np.arcsin(np.minimum(rp_over_rs / z_over_rs, 1))
    theta_2[case1] = np.pi
    theta_2[case2] = np.pi / 2.0
    theta_2[casea] = np.pi
    theta_2[casec] = np.pi / 2.0
    theta_2[case7] = ph[case7]

    # flux_upper
    plusflux = np.zeros(len(z_over_rs))
    plusflux[plus_case] = integral_plus_core(method, limb_darkening_coefficients, rp_over_rs, z_over_rs[plus_case],
                                             theta_1[plus_case], theta_2[plus_case], precision=precision)
    if len(case0[0]) > 0:
        plusflux[case0] = integral_centred(method, limb_darkening_coefficients, rp_over_rs, 0.0, np.pi)
    if len(caseb[0]) > 0:
        plusflux[caseb] = integral_centred(method, limb_darkening_coefficients, 1, 0.0, np.pi)

    # flux_lower
    minsflux = np.zeros(len(z_over_rs))
    minsflux[minus_case] = integral_minus_core(method, limb_darkening_coefficients, rp_over_rs,
                                               z_over_rs[minus_case], 0.0, theta_2[minus_case], precision=precision)

    # flux_star
    starflux = np.zeros(len(z_over_rs))
    starflux[star_case] = integral_centred(method, limb_darkening_coefficients, 1, 0.0, ph[star_case])

    # flux_total
    total_flux = integral_centred(method, limb_darkening_coefficients, 1, 0.0, 2.0 * np.pi)

    return 1 - (2.0 / total_flux) * (plusflux + starflux - minsflux)


# transit

def transit(limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron,
            mid_time, time_array, method='claret', precision=3):

    position_vector = planet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return transit_flux_drop(limb_darkening_coefficients, rp_over_rs, projected_distance,
                             method=method, precision=precision)


def transit_integrated(limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                       periastron, mid_time, time_array, exp_time, max_sub_exp_time=10, method='claret', precision=3):

    time_factor = int(exp_time / max_sub_exp_time) + 1
    exp_time /= (60.0 * 60.0 * 24.0)

    time_array_hr = (time_array[:, None] + np.arange(-exp_time / 2 + exp_time / time_factor / 2, exp_time / 2,
                                                     exp_time / time_factor)).flatten()

    position_vector = planet_orbit(period, sma_over_rs, eccentricity,
                                   inclination, periastron, mid_time, time_array_hr)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return np.mean(np.reshape(
        transit_flux_drop(limb_darkening_coefficients, rp_over_rs, projected_distance,
                          method=method, precision=precision),
        (len(time_array), time_factor)), 1)


def transit_duration(rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron):

    ww = periastron * np.pi / 180
    ii = inclination * np.pi / 180
    ee = eccentricity
    aa = sma_over_rs
    ro_pt = (1 - ee ** 2) / (1 + ee * np.sin(ww))
    b_pt = aa * ro_pt * np.cos(ii)
    if b_pt > 1:
        b_pt = 0.5
    s_ps = 1.0 + rp_over_rs
    df = np.arcsin(np.sqrt((s_ps ** 2 - b_pt ** 2) / ((aa ** 2) * (ro_pt ** 2) - b_pt ** 2)))
    aprox = (period * (ro_pt ** 2)) / (np.pi * np.sqrt(1 - ee ** 2)) * df * 60 * 60 * 24

    def function_to_fit(x, t):
        return planet_star_projected_distance(period, sma_over_rs, eccentricity, inclination, periastron,
                                              10000, np.array(10000 + t / 24 / 60 / 60))

    popt1, pcov1 = curve_fit(function_to_fit, [0], [1.0 + rp_over_rs], p0=[-aprox / 2])
    popt2, pcov2 = curve_fit(function_to_fit, [0], [1.0 + rp_over_rs], p0=[aprox / 2])

    return (popt2[0] - popt1[0]) / 24 / 60 / 60


def transit_t12(rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron):

    ww = periastron * np.pi / 180
    ii = inclination * np.pi / 180
    ee = eccentricity
    aa = sma_over_rs
    ro_pt = (1 - ee ** 2) / (1 + ee * np.sin(ww))
    b_pt = aa * ro_pt * np.cos(ii)
    if b_pt > 1:
        b_pt = 0.5
    s_ps = 1.0 + rp_over_rs
    df = np.arcsin(np.sqrt((s_ps ** 2 - b_pt ** 2) / ((aa ** 2) * (ro_pt ** 2) - b_pt ** 2)))
    aprox = (period * (ro_pt ** 2)) / (np.pi * np.sqrt(1 - ee ** 2)) * df * 60 * 60 * 24

    def function_to_fit(x, t):
        return planet_star_projected_distance(period, sma_over_rs, eccentricity, inclination, periastron,
                                              10000, np.array(10000 + t / 24 / 60 / 60))

    popt1, pcov1 = curve_fit(function_to_fit, [0], [1.0 + rp_over_rs], p0=[-aprox / 2])
    popt2, pcov2 = curve_fit(function_to_fit, [0], [1.0 - rp_over_rs], p0=[-aprox / 2])

    return min((popt2[0] - popt1[0]) / 24 / 60 / 60,
               0.5*transit_duration(rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron)
               )


def transit_depth(limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                  periastron, method='claret', precision=6):

    return 1 - transit(limb_darkening_coefficients, rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                       periastron, 10000, np.array([10000]), method=method, precision=precision)[0]


# eclipse

def eclipse(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron, mid_time,
            time_array, precision=3):

    position_vector = planet_orbit(period, sma_over_rs / rp_over_rs, eccentricity, inclination, periastron + 180,
                                   mid_time, time_array)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 / rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return (1.0 + fp_over_fs * transit_flux_drop([0, 0, 0, 0], 1 / rp_over_rs, projected_distance,
                                                 precision=precision, method='zero')) / (1.0 + fp_over_fs)


def eclipse_integrated(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron,
                       mid_time, time_array, exp_time, max_sub_exp_time=10, precision=3):

    time_factor = int(exp_time / max_sub_exp_time) + 1
    exp_time /= (60.0 * 60.0 * 24.0)

    time_array_hr = (time_array[:, None] + np.arange(-exp_time / 2 + exp_time / time_factor / 2, exp_time / 2,
                                                     exp_time / time_factor)).flatten()

    position_vector = planet_orbit(period, sma_over_rs / rp_over_rs, eccentricity, inclination, periastron + 180,
                                   mid_time, time_array_hr)

    projected_distance = np.where(
        position_vector[0] < 0, 1.0 + 5.0 * rp_over_rs,
        np.sqrt(position_vector[1] * position_vector[1] + position_vector[2] * position_vector[2]))

    return np.mean(np.reshape(
        (1.0 + fp_over_fs * transit_flux_drop([0, 0, 0, 0], 1 / rp_over_rs, projected_distance, method='zero',
                                              precision=precision)) / (1.0 + fp_over_fs),
        (len(time_array), time_factor)), 1)


def eclipse_mid_time(period, sma_over_rs, eccentricity, inclination, periastron, mid_time):
    test_array = np.arange(0, period, 0.001)
    xx, yy, zz = planet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time,
                              test_array + mid_time)

    test1 = np.where(xx < 0)
    yy = yy[test1]
    test_array = test_array[test1]

    aprox = test_array[np.argmin(np.abs(yy))]

    def function_to_fit(x, t):
        xx, yy, zz = planet_orbit(period, sma_over_rs, eccentricity, inclination, periastron, mid_time,
                            np.array(mid_time + t))
        return yy

    popt, pcov = curve_fit(function_to_fit, [0], [0], p0=[aprox])

    return mid_time + popt[0]


def eclipse_duration(rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron):

    ww = periastron * np.pi / 180
    ii = inclination * np.pi / 180
    ee = eccentricity
    aa = sma_over_rs
    ro_pt = (1 - ee ** 2) / (1 + ee * np.sin(ww))
    b_pt = aa * ro_pt * np.cos(ii)
    if b_pt > 1:
        b_pt = 0.5
    s_ps = 1.0 + rp_over_rs
    df = np.arcsin(np.sqrt((s_ps ** 2 - b_pt ** 2) / ((aa ** 2) * (ro_pt ** 2) - b_pt ** 2)))
    aprox = (period * (ro_pt ** 2)) / (np.pi * np.sqrt(1 - ee ** 2)) * df * 60 * 60 * 24

    emt = eclipse_mid_time(period, sma_over_rs, eccentricity, inclination, periastron, 10000)

    def function_to_fit(x, t):
        xx = planet_star_projected_distance(period, sma_over_rs, eccentricity, inclination, periastron,
                                              10000, np.array(emt + t / 24 / 60 / 60))
        return xx

    popt1, pcov1 = curve_fit(function_to_fit, [0], [1.0 + rp_over_rs], p0=[-aprox / 2])
    popt2, pcov2 = curve_fit(function_to_fit, [0], [1.0 + rp_over_rs], p0=[aprox / 2])

    return (popt2[0] - popt1[0]) / 24 / 60 / 60


def eclipse_depth(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity, inclination, periastron, precision=6):

    return 1 - eclipse(fp_over_fs, rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                       periastron, 10000, np.array([10000]), precision=precision)[0]


def _get_filter(photometric_filter):

    if photometric_filter not in plc_data.all_filters():
        raise PyLCInputError('{0} is not available. Available filters: {1}'.format(
            photometric_filter, ','.join(plc_data.all_filters())))


def fp_over_fs(rp_over_rs, sma_over_rs, albedo, emissivity, stellar_temperature, filter_name, wlrange=None):

    def _black_body(w, t):
        # w in mu
        w = w / (10 ** 6)
        h = 6.62607004 * (10 ** (-34))
        c = 3 * (10 ** 8)
        w5 = w ** 5
        k = 1.38064852 * (10 ** (-23))
        return (2 * h * c * c / w5) / (np.exp(h * c / w / k / t) - 1)

    planet_temperature = stellar_temperature * np.sqrt(0.5 / sma_over_rs) * (((1 - albedo) / emissivity) ** 0.25)

    if isinstance(filter_name, str):
        if filter_name in plc_data.all_filters():
            photometry_path = plc_data.photometry()
            filter_name = np.loadtxt(os.path.join(photometry_path, filter_name + '.pass'))

    passband_wlrange = [min(filter_name[:, 0]), max(filter_name[:, 0])]
    if not wlrange:
        wlrange = passband_wlrange
    else:
        wlrange = [min(wlrange), max(wlrange)]

        if wlrange[0] < passband_wlrange[0] or wlrange[1] > passband_wlrange[1]:
            raise PyLCInputError('Wavelength is not compatible with the filter.')

    sub_passband = filter_name[np.where((filter_name[:, 0] <= wlrange[1]) * (filter_name[:, 0] >= wlrange[0]))]
    wavelength_array, band = np.swapaxes(sub_passband, 0, 1)

    binsedge = 0.5 * (wavelength_array[:-1] + wavelength_array[1:])
    binsedge1 = np.append(wavelength_array[0] - (binsedge[0] - wavelength_array[0]), binsedge)
    binsedge2 = np.append(binsedge, wavelength_array[-1] + (wavelength_array[-1] - binsedge[-1]))
    binswidth = binsedge2 - binsedge1

    weights = band * binswidth / 10000
    emission = ((rp_over_rs ** 2) *
                _black_body(wavelength_array / 10000, planet_temperature) /
                _black_body(wavelength_array / 10000, stellar_temperature))

    emission = np.sum(emission * weights) / np.sum(weights)

    reflection = albedo * (rp_over_rs ** 2) / (sma_over_rs ** 2)

    return reflection + emission


def exotethys(stellar_logg, stellar_temperature, stellar_metallicity, filter_name, wlrange=None, method='claret', stellar_model='Phoenix_2018'):

    if 'phoenix' in stellar_model.lower() and stellar_metallicity != 0:
        warnings.warn('PHOENIX models are only computed for solar metallicity stars. Setting stellar_metallicity = 0.')
        stellar_metallicity = 0

    from exotethys import sail, ls_database

    path = plc_data.exotethys()

    available_stellar_models = ['Atlas_2000', 'Phoenix_2012_13', 'Phoenix_2018', 'Phoenix_drift_2012',
                                'Stagger_2015', 'Stagger_2018']
    if not isinstance(stellar_model, str) or stellar_model not in available_stellar_models:
        raise PyLCInputError('Stellar_model {0} is not available. Available models: {1}. '
                             'Please consult EXOTETHYS documentation for more details'.format(
            stellar_model, ','.join(available_stellar_models)))

    method_map = {
        'claret': 'claret4',
        'sqrt': 'square_root',
        'quad': 'quadratic',
        'linear': 'linear'
    }
    if method not in method_map:
        raise PyLCInputError('Method {0} is not valid. Available methods: {1}'.format(method, ','.join(list(method_map.keys()))))

    run_id = ''.join([str(np.random.randint(0, 99)) for ff in range(10)])
    bandpass_file = os.path.join(path, 'ww_{0}.pass'.format(run_id))
    parameters_file = os.path.join(path, 'ww_{0}_parameters.txt'.format(run_id))
    output_file = os.path.join(path, 'ww_{0}_ldc.pickle'.format(run_id))

    if isinstance(filter_name, str):
        if filter_name in plc_data.all_filters():
            photometry_path = plc_data.photometry()
            filter_name = np.loadtxt(os.path.join(photometry_path, filter_name + '.pass'))

    passband_wlrange = [min(filter_name[:, 0]), max(filter_name[:, 0])]
    if not wlrange:
        wlrange = passband_wlrange
    else:
        wlrange = [min(wlrange), max(wlrange)]

    if wlrange[0] < passband_wlrange[0] or wlrange[1] > passband_wlrange[1]:
        print('Wavelength range requested: {0}, Passband wavelength range: {1}.'.format(wlrange, passband_wlrange))
        raise PyLCInputError('Wavelength is not compatible with the filter.')

    sub_passband = filter_name[np.where((filter_name[:, 0] <= wlrange[1]) * (filter_name[:, 0] >= wlrange[0]))]

    np.savetxt(bandpass_file, sub_passband)

    r = open(os.path.join(path, 'ww_parameters.txt')).read()
    r = r.replace('{{output}}', path)
    r = r.replace('{{id}}', str(run_id))
    r = r.replace('{{model}}', stellar_model)
    r = r.replace('{{law}}', method_map[method])
    r = r.replace('{{teff}}', str(stellar_temperature))
    r = r.replace('{{logg}}', str(stellar_logg))
    r = r.replace('{{meta}}', str(stellar_metallicity))

    w = open(parameters_file, 'w')
    w.write(r)
    w.close()

    try:
        sail.ldc_calculate(parameters_file)

        results = open_dict(output_file)

        for ff in glob.glob(os.path.join(path, '*{0}*'.format(run_id))):
            os.remove(ff)

        ldcs = results['passbands']['ww_{0}.pass'.format(run_id)][method_map[method]]['coefficients']

        while len(ldcs) < 4:
            ldcs = np.append(ldcs, 0)

        return ldcs

    except pickle.UnpicklingError:
        for file in glob.glob(os.path.join(ls_database()[0], '*', '*')):
            try:
                _ = open_dict(file)
            except pickle.UnpicklingError:
                os.remove(file)
        return exotethys(stellar_logg, stellar_temperature, stellar_metallicity, filter_name, wlrange, method, stellar_model)
