
import platform
import numpy
import matplotlib
import scipy
import tkinter
import emcee
import astroquery
import sys
import astropy

tcl = tkinter.Tcl()
tkinter__version__ = tcl.call("info", "patchlevel")


def run_app():

    print('\nSystem setup:\n')

    print('OS: {0} {1}\n'.format(platform.system(), platform.release()))
    print('Python: {0}\n'.format(sys.version))
    print('Python paths: {0}\n'.format('\n              '.join(sys.path)))

    print('Tkinter version: {0}'.format(tkinter__version__))
    print('    @ {0}'.format(tkinter.__file__).replace('/__init__.py', ''))
    print('Numpy version: {0}'.format(numpy.__version__))
    print('    @ {0}'.format(numpy.__file__).replace('/__init__.py', ''))
    print('Scipy version: {0}'.format(scipy.__version__))
    print('    @ {0}'.format(scipy.__file__).replace('/__init__.py', ''))
    print('Matplotlib version: {0}'.format(matplotlib.__version__))
    print('    @ {0}'.format(matplotlib.__file__).replace('/__init__.py', ''))
    print('Emcee version: {0}'.format(emcee.__version__))
    print('    @ {0}'.format(emcee.__file__).replace('/__init__.py', ''))
    print('Astroquery version: {0}'.format(astroquery.__version__))
    print('    @ {0}'.format(astroquery.__file__).replace('/__init__.py', ''))
    print('Astropy version: {0}'.format(astropy.__version__))
    print('    @ {0}'.format(astropy.__file__).replace('/__init__.py', ''))

    print('\n\nLoading... Please wait for the main window to appear.')

    from .application import HOPS

    app = HOPS()
    app.run(f_after=app.log.check_for_updates)
