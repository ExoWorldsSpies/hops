
__version__ = '3.3.1'
__message__ = 'Version 3.3 is now online with new features!\n- UltraShortExposure mode, for stellar occultations.'

from .__run__ import run_app


def __get_abspath__():
    import os
    return os.path.abspath(os.path.dirname(__file__))
