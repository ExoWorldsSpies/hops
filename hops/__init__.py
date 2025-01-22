
__version__ = '3.3.3'
__message__ = 'Version 3.3 is now online with new features!\n- UltraShortExposure mode, for stellar occultations.\nv3.3.2 - Handling of saturated stars and better gary scale for plots\nv3.3.3 - Moving target mode.'

from .__run__ import run_app


def __get_abspath__():
    import os
    return os.path.abspath(os.path.dirname(__file__))
