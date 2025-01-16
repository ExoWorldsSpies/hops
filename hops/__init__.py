
__version__ = '3.3.1'
__message__ = 'Version 3.3 is now online with new features!\n- UltraShortExposure mode, for stellar occultations.\nv3.3.2 - Updated installation process\nv3.3.3 - Handling of saturated stars and better gary scale for plots.'

from .__run__ import run_app


def __get_abspath__():
    import os
    return os.path.abspath(os.path.dirname(__file__))
