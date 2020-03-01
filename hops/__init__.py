from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

__version__ = '2.5.7'
__message__ = 'Important - ExoClock database connected\nPhotometry update - suggesting comparison stars'

from .__run__ import run_app


def __get_abspath__():
    import os
    return os.path.abspath(os.path.dirname(__file__))
