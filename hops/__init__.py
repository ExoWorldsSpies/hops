
__version__ = '2.6.7'
__message__ = 'Fixed link to IERS earth rotation data.\nFixed Simbad astroquery.\nAdded filter for nan values.\nIgnoring trails in alignment'

from .__run__ import run_app


def __get_abspath__():
    import os
    return os.path.abspath(os.path.dirname(__file__))
