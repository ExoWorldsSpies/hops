from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import sys

if sys.version_info[0] > 2:
    from tkinter import *
    import tkinter.ttk as ttk
    import tkinter.filedialog as tkFileDialog
    from tkinter.messagebox import *
    from urllib.request import urlopen
else:
    import ttk
    from Tkinter import *
    import tkFileDialog
    from tkMessageBox import *
    from urllib import urlopen

import warnings
warnings.filterwarnings(
    'ignore', message='Matplotlib is building the font cache using fc-list. This may take a moment.')
warnings.filterwarnings(
    'ignore', message='The installed version of numexpr 2.4.4 is not supported in pandas and will be not be used')

import matplotlib
matplotlib.use('TkAgg')

import datetime
import os
import sys
import glob
import time
import yaml
import numpy as np
import shutil
import hops.pylightcurve3 as plc
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.patches as mpatch

from astropy.io import fits as pf
from scipy.optimize import curve_fit
from matplotlib.figure import Figure
from matplotlib.offsetbox import AnchoredText
from matplotlib.backend_bases import key_press_handler, MouseEvent
import matplotlib.gridspec as gridspec
try:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    NavigationToolbar2TkAgg = NavigationToolbar2Tk
except ImportError:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import FigureCanvasBase
import matplotlib.image as mpimg

from astroquery.simbad import Simbad
import requests
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation

import tkinter.scrolledtext as scrolledtext

import webbrowser

def openweb():
    webbrowser.open("https://www.exoworldsspies.com/en/software", new=1)

def openweb_simbad(radec_string):

    def mock(radec_string=radec_string):
        radec_string = radec_string.replace('+', '%2B').replace(' ', '+')
        webbrowser.open("http://simbad.u-strasbg.fr/simbad/sim-coo?Coord={0}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=20&Radius.unit=arcmin&submit=submit+query&CoordList=".format(radec_string), new=1)

    return mock

__location__ = os.path.abspath(os.path.dirname(__file__))
__home__ = os.path.expanduser('~')

holomon_logo = glob.glob(os.path.join(__location__, 'holomon.gif'))[0]
holomon_logo_jpg = mpimg.imread(glob.glob(os.path.join(__location__, 'logo.jpg'))[0])

log_file = os.path.join(__location__, 'log.yaml')
log_profile_file = os.path.join(__location__, 'log_profile.yaml')
log_user_file = os.path.join(__location__, 'log_user.yaml')

local_log_profile_file = os.path.join(__home__, 'hops_log_profile.yaml')
local_log_user_file = os.path.join(__home__, 'hops_log_user.yaml')
local_log_file = 'log.yaml'

if not os.path.isfile(local_log_profile_file):
    shutil.copy(log_profile_file, local_log_profile_file)

if not os.path.isfile(local_log_user_file):
    shutil.copy(log_user_file, local_log_user_file)
else:
    try:
        test = yaml.load(open(local_log_profile_file, 'r'), Loader=yaml.SafeLoader)
    except:
        shutil.copy(log_user_file, local_log_user_file)


def initiate_local_log_file():
    # if not os.path.isfile(local_log_file):
    #     shutil.copy(log_file, os.path.join(os.path.abspath('.'), local_log_file))
    try:
        test = yaml.load(open(local_log_file, 'r'), Loader=yaml.SafeLoader)
    except:
        shutil.copy(log_file, os.path.join(os.path.abspath('.'), local_log_file))


def read_local_log_user(keyword):
    x = yaml.load(open(local_log_user_file, 'r'), Loader=yaml.SafeLoader)
    return x[keyword]


def write_local_log_user(keyword, value):
    x = yaml.load(open(local_log_user_file, 'r'), Loader=yaml.SafeLoader)
    x[keyword] = value

    yaml.dump(x, open(local_log_user_file, 'w'), default_flow_style=False)


def read_log_profile(keyword):
    x = yaml.load(open(log_profile_file, 'r'), Loader=yaml.SafeLoader)
    if not x[keyword]:
        return ''
    else:
        test = str(x[keyword])
        if test[0] == ' ':
            test = test[1:]
        return test


def read_local_log_profile(keyword):
    try:
        x = yaml.load(open(local_log_profile_file, 'r'), Loader=yaml.SafeLoader)
        test = str(x[keyword])
        if len(str(test)) > 0:
            if test[0] == ' ':
                test = test[1:]
        return test
    except KeyError:
        return read_log_profile(keyword)


def write_local_log_profile(keyword, value):
    x = yaml.load(open(local_log_profile_file, 'r'), Loader=yaml.SafeLoader)
    x[keyword] = value

    yaml.dump(x, open(local_log_profile_file, 'w'), default_flow_style=False)


def read_log(keyword, keyword2=None):
    x = yaml.load(open(log_file, 'r'), Loader=yaml.SafeLoader)
    if keyword2:
        return x[keyword][keyword2]
    else:
        return x[keyword]


def read_local_log(keyword, keyword2=None):
    try:
        x = yaml.load(open(local_log_file, 'r'), Loader=yaml.SafeLoader)
        if keyword2:
            return x[keyword][keyword2]
        else:
            return x[keyword]
    except KeyError:
        try:
            value = read_log(keyword, keyword2)
            write_local_log(keyword, value, keyword2)
            return value
        except KeyError:
            return False


def write_local_log(keyword, value, keyword2=None):
    x = yaml.load(open(local_log_file, 'r'), Loader=yaml.SafeLoader)
    if keyword2:
        x[keyword][keyword2] = value
    else:
        x[keyword] = value

    yaml.dump(x, open(local_log_file, 'w'), default_flow_style=False)


def find_fits_files(fits_file):

    fits_list = glob.glob('*{0}*.f*t*'.format(fits_file)) + glob.glob('*{0}*.F*T*'.format(fits_file))
    fits_list = list(np.unique(fits_list))
    fits_list.sort()
    return fits_list


def test_fits_keyword(fits_file, keyword):

    if len(fits_file) == 0:
        return [False, 'No keyword found']

    else:
        try:
            fits_file = find_fits_files(fits_file)[0]

            fits = pf.open(fits_file)

            try:
                fits = [fits['SCI']]
            except KeyError:
                sci_id = 0
                for sci_id in range(len(fits)):
                    try:
                        if (fits[sci_id].data).all():
                            break
                    except:
                        pass
                fits = [fits[sci_id]]

            if fits[0].header[str(keyword)]:
                return [True, 'Keyword found', fits[0].header[str(keyword)]]

            else:
                return [False, 'No keyword found']

        except (KeyError, IndexError):
            return [False, 'No keyword found']


def test_file_number(fits_file):

    if len(fits_file) == 0:
        test = 0
    else:
        test = len(find_fits_files(fits_file))

    if test > 0:
        return [True, '{0} files found'.format(test)]
    else:
        return [False, 'No files found']


def test_coordinates(ra_dec_string):

    try:
        ra_dec_string = ra_dec_string.replace(':', ' ').split(' ')
        target = plc.Target(plc.Hours(*ra_dec_string[:3]), plc.Degrees(*ra_dec_string[3:]))
        return[True, 'Coordinates\naccepted']
    except:
        return [False, 'Wrong\ncoordinates']

# def initialise_window(window, window_name=None, exit_command=None):
#
#     if not window_name:
#         window_name = read_main_log('windows', 'software_window')
#
#     if not exit_command:
#         def exit_command():
#             os._exit(-1)
#
#     window.wm_title(window_name)
#     window.protocol('WM_DELETE_WINDOW', exit_command)
#
#     window.withdraw()
#
#
# def setup_window(window, objects):
#
#     main_font = tuple(read_main_log('windows', 'main_font'))
#     title_font = tuple(read_main_log('windows', 'title_font'))
#     button_font = tuple(read_main_log('windows', 'button_font'))
#     entries_bd = read_main_log('windows', 'entries_bd')
#
#     for row in range(len(objects)):
#         if len(objects[row]) == 0:
#             label_empty = Label(window, text='')
#             label_empty.grid(row=row, column=100)
#         else:
#             for obj in objects[row]:
#
#                 if obj[0].winfo_class() == 'Button':
#                     obj[0].configure(font=button_font)
#                 elif obj[0].winfo_class() == 'Entry':
#                     obj[0].configure(bd=entries_bd, font=main_font)
#                 elif obj[0].winfo_class() in ['Label', 'Radiobutton']:
#                     if len(obj) == 5:
#                         if obj[4] == 'title':
#                             obj[0].configure(font=title_font)
#                         else:
#                             obj[0].configure(font=main_font)
#                     else:
#                         obj[0].configure(font=main_font)
#
#                 if len(obj) >= 4:
#                     obj[0].grid(row=row, column=obj[1], columnspan=obj[2], rowspan=obj[3])
#                 elif len(obj) == 3:
#                     obj[0].grid(row=row, column=obj[1], columnspan=obj[2])
#                 else:
#                     obj[0].grid(row=row, column=obj[1])
#
#
# def finalise_window(window, center=True, topmost=False):
#
#     window.update_idletasks()
#
#     if center:
#         x = (window.winfo_screenwidth() - window.winfo_reqwidth()) / 2
#         y = (window.winfo_screenheight() - window.winfo_reqheight()) / 2
#         window.geometry('+%d+%d' % (x, y))
#
#     else:
#         window.geometry('+%d+%d' % (0, 0))
#
#     window.update_idletasks()
#
#     window.lift()
#     window.wm_attributes("-topmost", 1)
#     # if not topmost:
#     window.after_idle(window.attributes, '-topmost', 0)
#
#     window.deiconify()


def test_float_input(input_str, typing):

    if typing == '1':
        try:
            if float(input_str) >= 0:
                return True            
            elif float(input_str) < 0:
                return True
            else:
                return False
        except ValueError:
            return False

    else:
        return True


def test_float_positive_input(input_str, typing):

    if typing == '1':
        try:
            if float(input_str) >= 0:
                return True
            else:
                return False
        except ValueError:
            return False

    else:
        return True


def test_int_input(input_str, typing):

    if typing == '1':
        try:
            if int(input_str):
                return True
            else:
                return False
        except ValueError:
            return False

    else:
        return True


def test_int_positive_input(input_str, typing):

    if typing == '1':
        try:
            if int(input_str) >= 0:
                return True
            else:
                return False
        except ValueError:
            return False

    else:
        return True


def test_int_positive_non_zero_input(input_str, typing):

    if typing == '1':
        try:
            if int(input_str) > 0:
                return True
            else:
                return False
        except ValueError:
            return False

    else:
        return True


def get_reduction_alignment_help():
    return open(os.path.join(__location__, 'reduction_alignment_help.txt')).read()


def get_photometry_help():
    return open(os.path.join(__location__, 'photometry_help.txt')).read()


def get_fitting_help():
    return open(os.path.join(__location__, 'fitting_help.txt')).read()


filter_map = {'Clear': 'V', 'Luminance': 'V',
              'U': 'U', 'B': 'B', 'V': 'V', 'R': 'R', 'I': 'I', 'H': 'H', 'J': 'J', 'K': 'K',
              'u': 'u', 'b': 'b', 'v': 'v', 'y': 'y',
              'u\'': 'u,', 'g\'': 'g,', 'r\'': 'r,', 'i\'': 'i,', 'z\'': 'z,',
              'Astrodon ExoPlanet-BB': 'R',
              'UV': 'U', 'Rc': 'R', 'Ic': 'I', 'Re': 'R', 'Ie': 'I', 'Y': 'y,', 'r': 'r,', 'z': 'z,', 'i': 'i,',
              }
