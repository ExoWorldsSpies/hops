
__all__ = ['curve_fit', 'interp1d']

import warnings

from scipy.optimize import curve_fit as scipy_curve_fit
from scipy.interpolate import interp1d as scipy_interp1d


def curve_fit(*args, **kwargs):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message='Covariance of the parameters could not be estimated')
        return scipy_curve_fit(*args, **kwargs)


def interp1d(*args, **kwargs):
    return scipy_interp1d(*args, **kwargs)