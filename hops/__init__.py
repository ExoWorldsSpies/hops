
__version__ = '3.0.1'
__message__ = 'Version 3 is now online with many new features!\n3.0.1 - Add filter for overbloomig during alignment.'

from .__run__ import run_app


def __get_abspath__():
    import os
    return os.path.abspath(os.path.dirname(__file__))
