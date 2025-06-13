
__version__ = '3.3.6'
__message__ = 'Version 3.3 is now online with new features!\n3.3.1 - UltraShortExposure mode, for stellar occultations.\nv3.3.2 - Handling of saturated stars and better gray scale for plots\nv3.3.3 - Moving target mode.\nv3.3.4 - Colour camera mode.\nv.3.3.5 - Filter out frames with saturated pixels in Photometry.\nv.3.3.6 - Bypasss exotethys loading of paths with spaces.'

from .__run__ import run_app


def __get_abspath__():
    import os
    return os.path.abspath(os.path.dirname(__file__))
