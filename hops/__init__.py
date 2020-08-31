
__version__ = '2.6.2'
__message__ = 'Includes uncertainties in the final photometry output'

from .__run__ import run_app


def __get_abspath__():
    import os
    return os.path.abspath(os.path.dirname(__file__))
