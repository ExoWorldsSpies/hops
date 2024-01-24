
__version__ = '3.2.0'
__message__ = 'Version 3.2 is now online with new features, including connection with GAIA DR3!'

from .__run__ import run_app


def __get_abspath__():
    import os
    return os.path.abspath(os.path.dirname(__file__))
