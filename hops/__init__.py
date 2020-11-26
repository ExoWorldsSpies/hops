
__version__ = '2.6.6'
__message__ = 'Fixed link to IERS earth rotation data.\nFixed Simbad astroquery.\nFixed numpy version for windows\nAdded filter for nan values.'

from .__run__ import run_app


def __get_abspath__():
    import os
    return os.path.abspath(os.path.dirname(__file__))
