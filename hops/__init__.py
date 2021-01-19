
__version__ = '2.6.8'
__message__ = 'Fixed issue with matplotlib and math characters'

from .__run__ import run_app


def __get_abspath__():
    import os
    return os.path.abspath(os.path.dirname(__file__))
