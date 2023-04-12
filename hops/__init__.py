
__version__ = '3.1.1'
__message__ = 'Version 3.1 is now online with many new features!\n3.1.1 - Many bugs fixed, improved alignment, option to load location from profile, option to crop image edges.'

from .__run__ import run_app


def __get_abspath__():
    import os
    return os.path.abspath(os.path.dirname(__file__))
