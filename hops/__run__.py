



import warnings
warnings.filterwarnings(
    'ignore', message='Matplotlib is building the font cache using fc-list. This may take a moment.')
warnings.filterwarnings(
    'ignore', message='The installed version of numexpr 2.4.4 is not supported in pandas and will be not be used')
warnings.filterwarnings(
    'ignore', message='Covariance of the parameters could not be estimated')

import matplotlib
matplotlib.use('TkAgg')

import os
import glob
import numpy as np
from matplotlib.cm import Greys, Greys_r
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
import hops.pylightcurve3 as plc
from urllib.request import urlopen


from tkinter import Tk, TclError
from tkinter import Label, Button, Scale
from tkinter import DISABLED, NORMAL, END, RIGHT, LEFT, BOTH, Y, HORIZONTAL, StringVar
from tkinter.messagebox import showinfo

from hops.hops_tools.windows import MainWindow, AddOnWindow, SideWindow, setup_window, openweb_simbad
from hops.hops_tools.logs import log
from hops.hops_tools.tests import *
from hops.hops_tools.observing_planner import ObservingPlanner

from .reduction_routines import reduction as rdr_reduction
from .alignment_routines import alignment as alr_alignment
from .photometry_routines import photometry as phr_photometry
from .fitting_routines import fitting as ftr_fitting



class ReductionWindow(MainWindow):

    def __init__(self, run):

        self.run = run

        MainWindow.__init__(self, 'HOPS - Reduction & Alignment')

        self.show_content_window = AddOnWindow('Files list', 3, 3, 1)
        self.show_header_window = AddOnWindow('Header keywords list', 3, 3, 7)
        self.observing_planner = ObservingPlanner()

        # set variables, create and place widgets, main window

        self.observing_planner_button = self.Button(text='OBSERVATION\nPLANNER', command=self.observing_planner.show)

        self.update_directory = self.BooleanVar(False)
        self.directory = self.StringVar(log.read_local_log_user('directory'))
        self.directory_short = self.StringVar(log.read_local_log_user('directory_short'))
        self.directory_entry = self.Button(textvariable=self.directory_short, command=self.choose_directory)

        self.observation_files = self.StringVar(log.read_local_log_profile('observation_files'))
        self.observation_files_entry = self.Entry(textvariable=self.observation_files)
        self.observation_files_entry.bind(sequence='<KeyRelease>', func=self.update_window)
        self.observation_files_test = self.StringVar(' ')

        self.bias_files = self.StringVar(log.read_local_log_profile('bias_files'))
        self.bias_files_entry = self.Entry(textvariable=self.bias_files)
        self.bias_files_entry.bind(sequence='<KeyRelease>', func=self.update_window)
        self.bias_files_test = self.StringVar(' ')

        self.dark_files = self.StringVar(log.read_local_log_profile('dark_files'))
        self.dark_files_entry = self.Entry(textvariable=self.dark_files)
        self.dark_files_entry.bind(sequence='<KeyRelease>', func=self.update_window)
        self.dark_files_test = self.StringVar(' ')

        self.flat_files = self.StringVar(log.read_local_log_profile('flat_files'))
        self.flat_files_entry = self.Entry(textvariable=self.flat_files)
        self.flat_files_entry.bind(sequence='<KeyRelease>', func=self.update_window)
        self.flat_files_test = self.StringVar(' ')

        self.bin_fits = self.StringVar(log.read_local_log_profile('bin_fits'))
        self.bin_fits_entry = self.Entry(textvariable=self.bin_fits, validate='key')
        self.bin_fits_entry.bind(sequence='<KeyRelease>', func=self.update_window)
        self.bin_fits_entry['validatecommand'] = (self.bin_fits_entry.register(test_int_positive_non_zero_input),
                                                  '%P', '%d')

        self.show_files_button = self.Button(text='Show files', command=self.show_content_window.show)

        self.exposure_time_key = self.StringVar(log.read_local_log_profile('exposure_time_key'))
        self.exposure_time_key_entry = self.Entry(textvariable=self.exposure_time_key)
        self.exposure_time_key_entry.bind(sequence='<KeyRelease>', func=self.update_window)
        self.exposure_time_key_test = self.StringVar(' ')

        self.observation_date_key = self.StringVar(log.read_local_log_profile('observation_date_key'))
        self.observation_date_key_entry = self.Entry(textvariable=self.observation_date_key)
        self.observation_date_key_entry.bind(sequence='<KeyRelease>', func=self.update_window)
        self.observation_date_key_test = self.StringVar(' ')

        self.observation_time_key = self.StringVar(log.read_local_log_profile('observation_time_key'))
        self.observation_time_key_entry = self.Entry(textvariable=self.observation_time_key)
        self.observation_time_key_entry.bind(sequence='<KeyRelease>', func=self.update_window)
        self.observation_time_key_test = self.StringVar(' ')

        self.auto_target_ra_dec = self.StringVar(log.read_log('photometry', 'auto_target_ra_dec'))
        self.use_auto_target_ra_dec = self.BooleanVar(log.read_log('photometry', 'use_auto_target_ra_dec'))
        self.use_auto_target_ra_dec_entry = self.Checkbutton(text='Use detected values',
                                                             variable=self.use_auto_target_ra_dec,
                                                             command=self.update_window)

        self.target_ra_dec = self.StringVar(log.read_log('photometry', 'target_ra_dec'))
        self.target_ra_dec_entry = self.Entry(textvariable=self.target_ra_dec)
        self.target_ra_dec_entry.bind(sequence='<KeyRelease>', func=self.update_window)
        self.target_ra_dec_test = self.StringVar(' ')

        self.show_header_button = self.Button(text='Show header', command=self.show_header_window.show)

        self.mid_exposure = self.BooleanVar(log.read_log('photometry', 'mid_exposure'))
        self.mid_exposure_entry = self.Checkbutton(text='If the time stamp in your fits files refers to the '
                                                        'mid-exposure\ninstead of the exposure start, please tick '
                                                        'this box.',
                                                   variable=self.mid_exposure)

        self.run_reduction_alignment_button = self.Button(text='RUN REDUCTION & ALIGNMENT',
                                                          command=self.run_reduction_alignment)
        self.proceed_to_photometry_button = self.Button(text='PROCEED TO PHOTOMETRY',
                                                        command=self.proceed_to_photometry)
        self.proceed_to_fitting_button = self.Button(text='PROCEED TO FITTING', command=self.proceed_to_fitting)

        self.setup_window([
            [[self.get_logo(), 0, 1, 8]],
            [],
            [[self.Label(text='Reduction & Alignment'), 1, 3, 1, 'title']],
            [],
            [[self.Label(text='Directory'), 1], [self.directory_entry, 2, 2]],

            [[self.Label(text='Name identifier for observation files'), 1],
             [self.observation_files_entry, 2], [self.Label(textvar=self.observation_files_test), 3]],

            [[self.Label(text='Name identifier for bias files'), 1],
             [self.bias_files_entry, 2], [self.Label(textvar=self.bias_files_test), 3]],

            [[self.Label(text='Name identifier for dark files'), 1],
             [self.dark_files_entry, 2], [self.Label(textvar=self.dark_files_test), 3]],

            [[self.created_by_label, 0, 1, 3],
             [self.Label(text='Name identifier for flat files'), 1],
             [self.flat_files_entry, 2], [self.Label(textvar=self.flat_files_test), 3]],

            [[self.Label(text='Bin fits files (reduced only)'), 1], [self.bin_fits_entry, 2]],

            [[self.show_files_button, 2]],
            [[self.my_profile_button, 0]],

            [[self.Label(text='Detected target RA DEC'), 1],
             [self.Label(textvar=self.auto_target_ra_dec), 2], [self.use_auto_target_ra_dec_entry, 3]],
            [[self.Btn, 0],
             [self.Label(text='Manual target RA DEC\n(hh:mm:ss +/-dd:mm:ss)'), 1],
             [self.target_ra_dec_entry, 2], [self.Label(textvar=self.target_ra_dec_test), 3]],

            [[self.Label(text='Exposure time header keyword'), 1],
             [self.exposure_time_key_entry, 2], [self.Label(textvar=self.exposure_time_key_test), 3]],

            [[self.Label(text='Observation date header keyword'), 1],
             [self.observation_date_key_entry, 2], [self.Label(textvar=self.observation_date_key_test), 3]],

            [[self.Label(text='Observation time header keyword'), 1],
             [self.observation_time_key_entry, 2], [self.Label(textvar=self.observation_time_key_test), 3]],

            [[self.mid_exposure_entry, 1, 3]],
            [[self.show_header_button, 2]],

            [[self.Label(text='Extra tools:'), 0]],
            [[self.observing_planner_button, 0, 1, 2], [self.run_reduction_alignment_button, 1, 3]],
            [[self.proceed_to_photometry_button, 1], [self.proceed_to_fitting_button, 2, 2]],
            [],
        ])

        # set variables, create and place widgets, side windows

        scrollbar = self.show_content_window.Scrollbar()
        scrollbar.pack(side=RIGHT, fill=Y)
        self.files_list = self.show_content_window.Listbox(yscrollcommand=scrollbar.set, font='Courier')
        self.files_list.pack(side=LEFT, fill=BOTH, expand=True)
        scrollbar.config(command=self.files_list.yview)

        scrollbar = self.show_header_window.Scrollbar()
        scrollbar.pack(side=RIGHT, fill=Y)
        self.header_list = self.show_header_window.Listbox(yscrollcommand=scrollbar.set, font='Courier')
        self.header_list.pack(side=LEFT, fill=BOTH, expand=True)
        scrollbar.config(command=self.header_list.yview)

        # run

        self.update_directory.set(True)
        self.update_window(None)
        self.update_directory.set(False)

        self.run.check_for_update()
        self.loop()

    # define functions

    def update_window(self, *entry):

        if not entry:
            pass

        if self.running.get():

            self.directory_entry['state'] = DISABLED
            self.observation_files_entry['state'] = DISABLED
            self.bias_files_entry['state'] = DISABLED
            self.dark_files_entry['state'] = DISABLED
            self.flat_files_entry['state'] = DISABLED
            self.bin_fits_entry['state'] = DISABLED
            self.show_files_button['state'] = DISABLED
            self.exposure_time_key_entry['state'] = DISABLED
            self.observation_date_key_entry['state'] = DISABLED
            self.observation_time_key_entry['state'] = DISABLED
            self.target_ra_dec_entry['state'] = DISABLED
            self.use_auto_target_ra_dec_entry['state'] = DISABLED
            self.mid_exposure_entry['state'] = DISABLED
            self.show_header_button['state'] = DISABLED
            self.run_reduction_alignment_button['state'] = DISABLED
            self.proceed_to_photometry_button['state'] = DISABLED
            self.proceed_to_fitting_button['state'] = DISABLED
            self.observing_planner_button['state'] = DISABLED
            self.my_profile_button['state'] = DISABLED

        elif not os.path.isdir(self.directory.get()):

            self.directory_entry['state'] = NORMAL
            self.observation_files_entry['state'] = DISABLED
            self.bias_files_entry['state'] = DISABLED
            self.dark_files_entry['state'] = DISABLED
            self.flat_files_entry['state'] = DISABLED
            self.bin_fits_entry['state'] = DISABLED
            self.show_files_button['state'] = DISABLED
            self.exposure_time_key_entry['state'] = DISABLED
            self.observation_date_key_entry['state'] = DISABLED
            self.observation_time_key_entry['state'] = DISABLED
            self.target_ra_dec_entry['state'] = DISABLED
            self.use_auto_target_ra_dec_entry['state'] = DISABLED
            self.mid_exposure_entry['state'] = DISABLED
            self.show_header_button['state'] = DISABLED
            self.run_reduction_alignment_button['state'] = DISABLED
            self.proceed_to_photometry_button['state'] = DISABLED
            self.proceed_to_fitting_button['state'] = DISABLED
            self.observing_planner_button['state'] = NORMAL
            self.my_profile_button['state'] = NORMAL

        else:

            self.directory_entry['state'] = NORMAL

            if self.update_directory.get():

                os.chdir(self.directory.get())
                log.initiate_local_log_file()

                log.write_local_log('directory', self.directory.get())
                log.write_local_log('directory_short', self.directory_short.get())

                self.observation_files.set(log.read_local_log('pipeline', 'observation_files'))
                self.bias_files.set(log.read_local_log('reduction', 'bias_files'))
                self.dark_files.set(log.read_local_log('reduction', 'dark_files'))
                self.flat_files.set(log.read_local_log('reduction', 'flat_files'))
                self.bin_fits.set(log.read_local_log('reduction', 'bin_fits'))
                self.exposure_time_key.set(log.read_local_log('pipeline_keywords', 'exposure_time_key'))
                self.observation_date_key.set(log.read_local_log('pipeline_keywords', 'observation_date_key'))
                self.observation_time_key.set(log.read_local_log('pipeline_keywords', 'observation_time_key'))
                self.target_ra_dec.set(log.read_local_log('photometry', 'target_ra_dec'))
                self.auto_target_ra_dec.set(log.read_local_log('photometry', 'auto_target_ra_dec'))
                self.use_auto_target_ra_dec.set(log.read_local_log('photometry', 'use_auto_target_ra_dec'))
                self.mid_exposure.set(log.read_local_log('photometry', 'mid_exposure'))

            self.observation_files_entry['state'] = NORMAL
            self.bias_files_entry['state'] = NORMAL
            self.dark_files_entry['state'] = NORMAL
            self.flat_files_entry['state'] = NORMAL
            self.bin_fits_entry['state'] = NORMAL
            self.show_files_button['state'] = NORMAL
            self.observing_planner_button['state'] = NORMAL
            self.my_profile_button['state'] = NORMAL

            self.files_list.delete(0, END)
            self.files_list.insert(END, '  List of files in your directory:')
            self.files_list.insert(END, '  ')

            xx = find_fits_files('*')

            for ii in xx:
                self.files_list.insert(END, '  {0}'.format(str(ii).split(os.sep)[-1]))

            check_science = test_file_number(self.observation_files_entry.get())
            self.observation_files_test.set(check_science[1])

            check_bias = test_file_number(self.bias_files_entry.get())
            self.bias_files_test.set(check_bias[1])

            check_dark = test_file_number(self.dark_files_entry.get())
            self.dark_files_test.set(check_dark[1])

            check_flat = test_file_number(self.flat_files_entry.get())
            self.flat_files_test.set(check_flat[1])

            if not check_science[0]:

                self.exposure_time_key_entry['state'] = DISABLED
                self.observation_date_key_entry['state'] = DISABLED
                self.observation_time_key_entry['state'] = DISABLED
                self.target_ra_dec_entry['state'] = DISABLED
                self.use_auto_target_ra_dec_entry['state'] = DISABLED
                self.mid_exposure_entry['state'] = DISABLED
                self.show_header_button['state'] = DISABLED
                self.run_reduction_alignment_button['state'] = DISABLED

            else:

                self.target_ra_dec_entry['state'] = NORMAL
                self.use_auto_target_ra_dec_entry['state'] = NORMAL
                self.mid_exposure_entry['state'] = NORMAL
                self.exposure_time_key_entry['state'] = NORMAL
                self.observation_date_key_entry['state'] = NORMAL
                self.observation_time_key_entry['state'] = NORMAL
                self.show_header_button['state'] = NORMAL

                fits = plc.open_fits(find_fits_files(self.observation_files_entry.get())[0])

                try:
                    fits = [fits['SCI']]
                except KeyError:
                    sci_id = 0
                    for sci_id in range(len(fits)):
                        try:
                            if fits[sci_id].data.all():
                                break
                        except:
                            pass
                    fits = [fits[sci_id]]

                science_header = fits[0].header

                self.header_list.delete(0, END)
                self.header_list.insert(END, '  Keywords:      Values:')
                self.header_list.insert(END, '  ')

                for ii in science_header:

                    if ii != '':
                        self.header_list.insert(END, '  {0}{1}{2}'.format(str(ii[:10]), ' ' * (15 - len(str(ii[:10]))),
                                                                     str(science_header[ii])))

                check_ra = [None, None]
                for key in log.read_local_log_profile('target_ra_key').split(','):
                    check_ra = test_fits_keyword(self.observation_files_entry.get(), key)
                    if check_ra[0]:
                        break

                check_dec = [None, None]
                for key in log.read_local_log_profile('target_dec_key').split(','):
                    check_dec = test_fits_keyword(self.observation_files_entry.get(), key)
                    if check_dec[0]:
                        break

                if check_ra[0] and check_dec[0]:
                    try:
                        if isinstance(check_ra[2], str):
                            target = plc.Target(plc.Hours(check_ra[2].replace(',', '.')),
                                                plc.Degrees(check_dec[2].replace(',', '.')))
                            self.auto_target_ra_dec.set(target.coord)
                            self.use_auto_target_ra_dec_entry['state'] = NORMAL
                        elif isinstance(check_ra[2], float):
                            target = plc.Target(plc.Degrees(check_ra[2]), plc.Degrees(check_dec[2]))
                            self.auto_target_ra_dec.set(target.coord)
                            self.use_auto_target_ra_dec_entry['state'] = NORMAL
                        else:
                            self.auto_target_ra_dec.set('None detected')
                            self.use_auto_target_ra_dec.set(0)
                            self.use_auto_target_ra_dec_entry['state'] = DISABLED
                    except:
                        self.auto_target_ra_dec.set('None detected')
                        self.use_auto_target_ra_dec.set(0)
                        self.use_auto_target_ra_dec_entry['state'] = DISABLED
                else:
                    self.auto_target_ra_dec.set('None detected')
                    self.use_auto_target_ra_dec.set(0)
                    self.use_auto_target_ra_dec_entry['state'] = DISABLED

                if self.use_auto_target_ra_dec.get():
                    self.target_ra_dec.set(self.auto_target_ra_dec.get())
                    self.target_ra_dec_entry['state'] = DISABLED
                else:
                    self.target_ra_dec_entry['state'] = NORMAL

                check_ra_dec = test_coordinates(self.target_ra_dec_entry.get())
                self.target_ra_dec_test.set(check_ra_dec[1])

                check_exposure_time = test_fits_keyword(self.observation_files_entry.get(),
                                                        self.exposure_time_key_entry.get())
                self.exposure_time_key_test.set(check_exposure_time[1])

                check_observation_date = test_fits_keyword(self.observation_files_entry.get(),
                                                           self.observation_date_key_entry.get())
                self.observation_date_key_test.set(check_observation_date[1])

                if check_observation_date[0]:
                    if len(check_observation_date[2].split('T')) == 2:
                        self.observation_time_key.set(self.observation_date_key_entry.get())
                        self.observation_time_key_entry['state'] = DISABLED

                check_observation_time = test_fits_keyword(self.observation_files_entry.get(),
                                                           self.observation_time_key_entry.get())
                self.observation_time_key_test.set(check_observation_time[1])

                if (check_ra_dec[0] and check_exposure_time[0] and
                        check_observation_date[0] and check_observation_time[0]):
                    self.run_reduction_alignment_button['state'] = NORMAL

                else:
                    self.run_reduction_alignment_button['state'] = DISABLED

                if log.read_local_log('pipeline', 'reduction_complete') and log.read_local_log('pipeline',
                                                                                               'alignment_complete'):
                    try:
                        all_stars = plc.open_dict('all_stars.pickle')
                        in_fov = all_stars['in_fov']
                        all_stars = all_stars['all_stars']
                        self.proceed_to_photometry_button['state'] = NORMAL
                    except:
                        self.proceed_to_photometry_button['state'] = DISABLED
                else:
                    self.proceed_to_photometry_button['state'] = DISABLED

                if log.read_local_log('pipeline', 'photometry_complete'):
                    self.proceed_to_fitting_button['state'] = NORMAL
                else:
                    self.proceed_to_fitting_button['state'] = DISABLED

    # define actions for the different buttons, including calls to the function that updates the window

    def choose_directory(self):

        new_directory = self.askdirectory()

        if len(new_directory) > 0:
            self.directory.set(new_directory)
            self.directory_short.set('/'.join(new_directory.split('/')[-2:]))
            log.write_local_log_user('directory', self.directory.get())
            log.write_local_log_user('directory_short', self.directory_short.get())
            self.directory_entry.update()

            self.update_directory.set(True)
            self.update_window('a')
            self.update_directory.set(False)

    def proceed_to_photometry(self):
        self.run.run_from_reduction = False
        self.run.run_from_photometry = True
        self.run.run_from_fitting = False
        self.show_content_window.close()
        self.show_header_window.close()
        self.my_profile_window.close()
        self.close()

    def proceed_to_fitting(self):
        self.run.run_from_reduction = False
        self.run.run_from_photometry = False
        self.run.run_from_fitting = True
        self.show_content_window.close()
        self.show_header_window.close()
        self.my_profile_window.close()
        self.close()

    def run_reduction_alignment(self):

        self.running.set(True)
        self.update_window(None)

        log.write_local_log('pipeline', self.directory.get(), 'directory')
        log.write_local_log_user('directory', self.directory.get())
        log.write_local_log('pipeline', self.directory_short.get(), 'directory_short')
        log.write_local_log_user('directory_short', self.directory_short.get())
        log.write_local_log('pipeline', self.observation_files.get(), 'observation_files')
        log.write_local_log('reduction', self.bias_files.get(), 'bias_files')
        log.write_local_log('reduction', self.dark_files.get(), 'dark_files')
        log.write_local_log('reduction', self.flat_files.get(), 'flat_files')
        try:
            bins_test = int(self.bin_fits.get())
        except:
            bins_test = 1
        log.write_local_log('reduction', bins_test, 'bin_fits')
        log.write_local_log('photometry', self.target_ra_dec.get(), 'target_ra_dec')
        log.write_local_log('photometry', self.use_auto_target_ra_dec.get(), 'use_auto_target_ra_dec')
        log.write_local_log('photometry', self.mid_exposure.get(), 'mid_exposure')
        log.write_local_log('pipeline_keywords', self.exposure_time_key.get(), 'exposure_time_key')
        log.write_local_log('pipeline_keywords', self.observation_date_key.get(), 'observation_date_key')
        log.write_local_log('pipeline_keywords', self.observation_time_key.get(), 'observation_time_key')

        rdr_reduction()

        if log.read_local_log('pipeline', 'reduction_complete'):
            alr_alignment()

        if log.read_local_log('pipeline', 'reduction_complete') and log.read_local_log('pipeline',
                                                                                       'alignment_complete'):
            self.run.run_from_reduction = False
            self.run.run_from_photometry = True
            self.run.run_from_fitting = False
            self.show_content_window.close()
            self.show_header_window.close()
            self.my_profile_window.close()
            self.close()
        else:
            self.running.set(False)
            self.update_window(None)


class PhotometryWindow(MainWindow):

    def __init__(self, run):

        self.run = run
        MainWindow.__init__(self, 'HOPS - Photometry', position=11)
        self.show_fov_window = SideWindow('FOV', position=10)
        frame1 = self.show_fov_window.Frame()

        # set variables, create and place widgets, main window

        self.reduction_directory = log.read_local_log('pipeline', 'reduction_directory')
        self.light_curve_aperture_file = log.read_local_log('pipeline', 'light_curve_aperture_file')
        self.photometry_directory = log.read_local_log('pipeline', 'photometry_directory')
        self.fov_figure = log.read_local_log('pipeline', 'fov_figure')
        self.mean_key = log.read_local_log('pipeline_keywords', 'mean_key')
        self.std_key = log.read_local_log('pipeline_keywords', 'std_key')
        self.align_x0_key = log.read_local_log('pipeline_keywords', 'align_x0_key')
        self.align_y0_key = log.read_local_log('pipeline_keywords', 'align_y0_key')
        self.frame_low_std = log.read_local_log('windows', 'frame_low_std')
        self.frame_upper_std = log.read_local_log('windows', 'frame_upper_std')
        self.bin_fits = int(log.read_local_log('reduction', 'bin_fits'))
        self.burn_limit = int(log.read_local_log('alignment', 'burn_limit')) * self.bin_fits * self.bin_fits
        self.star_std = log.read_local_log('alignment', 'star_std')
        self.star_psf = log.read_local_log('alignment', 'star_psf')
        self.max_comparisons = log.read_local_log('photometry', 'max_comparisons')
        self.max_targets = self.max_comparisons + 1
        self.target_ra_dec = log.read_local_log('photometry', 'target_ra_dec')
        self.visible_fov_x_min = log.read_local_log('alignment', 'min_x')
        self.visible_fov_y_min = log.read_local_log('alignment', 'min_y')
        self.visible_fov_x_max = log.read_local_log('alignment', 'max_x')
        self.visible_fov_y_max = log.read_local_log('alignment', 'max_y')

        self.all_stars = plc.open_dict('all_stars.pickle')
        self.in_fov = self.all_stars['in_fov']
        self.all_stars = self.all_stars['all_stars']

        self.fits = plc.open_fits(find_fits_files(os.path.join(self.reduction_directory, '*'))[0])

        self.targets_indication = self.IntVar(0)
        self.targets_indication_entry = [
            self.Radiobutton(text='      Target           ', variable=self.targets_indication, value=0)]

        self.targets_x_position = [self.DoubleVar(log.read_local_log('photometry', 'target_x_position'))]
        self.targets_x_position_label = [self.Label(textvar=self.targets_x_position[0])]

        self.targets_y_position = [self.DoubleVar(log.read_local_log('photometry', 'target_y_position'))]
        self.targets_y_position_label = [self.Label(textvar=self.targets_y_position[0])]

        self.targets_peak_counts = [self.DoubleVar(0)]
        self.targets_peak_counts_label = [self.Label(textvar=self.targets_peak_counts[0])]

        self.targets_clear_button = [self.close_buttons[0]]

        self.targets_aperture = [self.DoubleVar(log.read_local_log('photometry', 'target_aperture'))]
        self.targets_aperture_entry = [self.Entry(textvar=self.targets_aperture[0], validate='key')]

        self.targets_warning = [self.StringVar('')]
        self.targets_warning_label = [self.Label(textvar=self.targets_warning[0])]

        self.targets_flux = [self.DoubleVar(0)]

        for comparison in range(self.max_comparisons):

            self.targets_indication_entry.append(self.Radiobutton(text='Comparison {0}     '.format(comparison + 1),
                                                                  variable=self.targets_indication,
                                                                  value=comparison + 1))

            self.targets_x_position.append(
                self.DoubleVar(log.read_local_log('photometry', 'comparison_{0}_x_position'.format(comparison + 1))))
            self.targets_x_position_label.append(self.Label(textvar=self.targets_x_position[comparison + 1]))

            self.targets_y_position.append(
                self.DoubleVar(log.read_local_log('photometry', 'comparison_{0}_y_position'.format(comparison + 1))))
            self.targets_y_position_label.append(self.Label(textvar=self.targets_y_position[comparison + 1]))

            self.targets_peak_counts.append(self.DoubleVar(0))
            self.targets_peak_counts_label.append(self.Label(textvar=self.targets_peak_counts[comparison + 1]))

            self.targets_clear_button.append(self.close_buttons[comparison + 1])

            self.targets_aperture.append(
                self.DoubleVar(log.read_local_log('photometry', 'comparison_{0}_aperture'.format(comparison + 1))))
            self.targets_aperture_entry.append(self.Entry(textvar=self.targets_aperture[comparison + 1],
                                                          validate='key'))

            self.targets_warning.append(self.StringVar(''))
            self.targets_warning_label.append(self.Label(textvar=self.targets_warning[comparison + 1]))

            self.targets_flux.append(self.DoubleVar(0))

        for target in range(self.max_targets):

            self.targets_aperture_entry[target]['validatecommand'] = \
                (self.targets_aperture_entry[target].register(test_float_positive_input), '%P', '%d')

            try:
                test = self.targets_aperture[target].get()
            except TclError:
                self.targets_x_position[target].set(0)
                self.targets_y_position[target].set(0)
                self.targets_aperture[target].set(0)

        self.photometry_button = self.Button(text='RUN PHOTOMETRY')
        self.proceed_to_fitting_button = self.Button(text='PROCEED TO FITTING')
        self.return_to_reduction_button = self.Button(text='RETURN TO REDUCTION')

        self.Btn2 = self.Button(text="CHECK SIMBAD", command=openweb_simbad(self.target_ra_dec))

        setup_list = [
            [[self.get_logo(), 0, 1, 8]],
            [],
            [[self.Label(text='Photometry'), 1, 6, 1, 'title']],
            [[self.Label(text="Remember, the best comparison stars need to be:\n"
                              "a) close to your target, b) of similar magnitude to the target,\n"
                              "c) of similar colour to the target,\nd) photometrically stable, i.e. "
                              "not variables!"), 1, 6]],
            [[self.Btn2, 1, 6]],
            [],
            [[self.Label(text='X'), 2], [self.Label(text='Y'), 3], [self.Label(text='Peak'), 4],
             [self.Label(text='Apert. radius'), 6], [self.Label(text='    WARNINGS    '), 7, 2]],
        ]

        for target in range(self.max_targets):

            if target == 1:
                setup_list.append([[self.created_by_label, 0, 1, 3],
                                   [self.targets_indication_entry[target], 1],
                                   [self.targets_x_position_label[target], 2],
                                   [self.targets_y_position_label[target], 3],
                                   [self.targets_peak_counts_label[target], 4],
                                   [self.targets_aperture_entry[target], 6],
                                   [self.targets_clear_button[target], 5],
                                   [self.targets_warning_label[target], 7]])
            elif target == 3:
                setup_list.append([[self.my_profile_button, 0, 1, 3],
                                   [self.targets_indication_entry[target], 1],
                                   [self.targets_x_position_label[target], 2],
                                   [self.targets_y_position_label[target], 3],
                                   [self.targets_peak_counts_label[target], 4],
                                   [self.targets_aperture_entry[target], 6],
                                   [self.targets_clear_button[target], 5],
                                   [self.targets_warning_label[target], 7]])

            elif target == 5:
                setup_list.append([[self.Btn, 0, 1, 3],
                                   [self.targets_indication_entry[target], 1],
                                   [self.targets_x_position_label[target], 2],
                                   [self.targets_y_position_label[target], 3],
                                   [self.targets_peak_counts_label[target], 4],
                                   [self.targets_aperture_entry[target], 6],
                                   [self.targets_clear_button[target], 5],
                                   [self.targets_warning_label[target], 7]])

            else:
                setup_list.append([[self.targets_indication_entry[target], 1],
                                   [self.targets_x_position_label[target], 2],
                                   [self.targets_y_position_label[target], 3],
                                   [self.targets_peak_counts_label[target], 4],
                                   [self.targets_aperture_entry[target], 6],
                                   [self.targets_clear_button[target], 5],
                                   [self.targets_warning_label[target], 7]])

        # setup_list.append([[show_fov_button, 1]])
        setup_list.append([])
        setup_list.append([[self.photometry_button, 2, 3]])
        setup_list.append([[self.proceed_to_fitting_button, 1, 2], [self.return_to_reduction_button, 4, 3]])
        setup_list.append([])

        self.setup_window(setup_list, entries_wd=10)

        # set variables, create and place widgets, main window

        self.plot_black = self.DoubleVar(0)
        self.plot_white = self.DoubleVar(1)

        fov_help_label = Label(frame1, text='- To select a star, double click on it.\n'
                                            '- To zoom in or out use your mouse or touchpad scroll.')

        self.flip_fov_button = Button(frame1, text='Flip FOV')
        self.mirror_fov_button = Button(frame1, text='Mirror FOV')
        self.reset_fov_button = Button(frame1, text='RESET PLOT')
        self.reverse_color_button = Button(frame1, text='INVERT\nBLACK & WHITE')

        self.f = Figure()
        self.f.patch.set_facecolor('white')
        self.ax = self.f.add_subplot(111)
        self.f.subplots_adjust(0.05, 0.08, 1 - 0.005, 1 - 0.13)
        self.ax.set_aspect(aspect="auto")
        self.canvas = self.show_fov_window.FigureCanvasTkAgg(self.f)
        self.canvas.get_tk_widget().pack()
        toolbar = self.show_fov_window.NavigationToolbar2Tk(self.canvas)
        toolbar.pack()

        self.image = self.ax.imshow(self.fits[1].data, origin='lower', extent=(0, len(self.fits[1].data[0]), 0,
                                                                               len(self.fits[1].data)),
                                    cmap=Greys_r,
                                    vmin=(self.fits[1].header[self.mean_key] +
                                          self.frame_low_std * self.fits[1].header[self.std_key]),
                                    vmax=(self.fits[1].header[self.mean_key] +
                                          self.frame_upper_std * self.fits[1].header[self.std_key]))

        x_length = len(self.fits[1].data[0])
        y_length = len(self.fits[1].data)
        self.circles_radius = 0.03 * max(y_length, x_length)

        self.ax.add_patch(mpatches.Rectangle((self.visible_fov_x_min + 1, self.visible_fov_y_min + 1),
                                             self.visible_fov_x_max - self.visible_fov_x_min - 2,
                                             self.visible_fov_y_max - self.visible_fov_y_min - 2,
                                             ec='r', fill=False, label='Available FOV'))

        self.good_comps_boxes1 = []
        self.good_comps_boxes2 = []

        self.targets_box = [mpatches.Circle((self.targets_x_position[0].get(), self.targets_y_position[0].get()),
                                       self.targets_aperture[0].get(),
                                       ec='r', fill=False)]
        for comparison in range(self.max_comparisons):
            box = mpatches.Circle((self.targets_x_position[comparison + 1].get(),
                                   self.targets_y_position[comparison + 1].get()),
                                  self.targets_aperture[comparison + 1].get(),
                                  ec='#07fefc', fill=False)

            if comparison == 0:
                circle1 = mpatches.Circle((-1000, -1000), self.circles_radius, ec='y', fill=False,
                                          label='Stars of similar flux to the target (+/- 40%)')
            else:
                circle1 = mpatches.Circle((-1000, -1000), self.circles_radius, ec='y', fill=False,)
            circle2 = mpatches.Circle((-1000, -1000), 0.75 * self.circles_radius, ec='y', fill=False)

            self.targets_box.append(box)
            self.good_comps_boxes1.append(circle1)
            self.good_comps_boxes2.append(circle2)

        for box in self.targets_box:
            self.ax.add_patch(box)

        for circle1 in self.good_comps_boxes1:
            self.ax.add_patch(circle1)

        for circle2 in self.good_comps_boxes2:
            self.ax.add_patch(circle2)

        self.targets_text = [self.ax.text(self.targets_x_position[0].get(),
                             self.targets_y_position[0].get() - self.targets_aperture[0].get() - 1, 'T',
                             color='r', fontsize=20, va='top')]

        for comparison in range(self.max_comparisons):
            self.targets_text.append(self.ax.text(self.targets_x_position[comparison + 1].get()
                                                  + self.targets_aperture[comparison + 1].get() + 1,
                                                  self.targets_y_position[comparison + 1].get()
                                                  - self.targets_aperture[comparison + 1].get() - 1,
                                                  'C{0}'.format(comparison + 1), color='#07fefc',
                                                  fontsize=20, va='top'))

        self.ax.legend(loc=(0, 1.01))

        self.update_idletasks()
        self.plot_black_entry = Scale(frame1, length=self.root.winfo_reqwidth()/4.0,
                                      from_=max(0, np.min(self.fits[1].data)),
                                      to=(self.fits[1].header[self.mean_key] + 200 * self.fits[1].header[self.std_key]),
                                      resolution=1, variable=self.plot_black, orient=HORIZONTAL)
        self.plot_white_entry = Scale(frame1, length=self.root.winfo_reqwidth()/4.0,
                                      from_=max(0, np.min(self.fits[1].data)),
                                      to=(self.fits[1].header[self.mean_key] + 200 * self.fits[1].header[self.std_key]),
                                      resolution=1, variable=self.plot_white, orient=HORIZONTAL)

        self.plot_black_entry.set(self.fits[1].header[self.mean_key] +
                                  self.frame_low_std * self.fits[1].header[self.std_key])
        self.plot_white_entry.set(self.fits[1].header[self.mean_key] +
                                  self.frame_upper_std * self.fits[1].header[self.std_key])

        self.plot_black.set(self.fits[1].header[self.mean_key] +
                            self.frame_low_std * self.fits[1].header[self.std_key])
        self.plot_white.set(self.fits[1].header[self.mean_key] +
                            self.frame_upper_std * self.fits[1].header[self.std_key])

        self.min_label = StringVar(frame1, value='Min (black pixels)')
        self.max_label = StringVar(frame1, value='Max (white pixels)')

        setup_window(frame1, [
            [[fov_help_label, 0, 6]],
            [],
            [[self.reset_fov_button, 0, 6]],
            [[self.flip_fov_button, 0], [self.mirror_fov_button, 1], [self.plot_black_entry, 2, 2],
             [self.plot_white_entry, 4, 2]],
            [[self.reverse_color_button, 0, 2], [Label(frame1, textvar=self.min_label), 2, 2], [Label(frame1, textvar=self.max_label), 4, 2]],
            [],
        ])

        frame1.pack()

        self.update_window(None)

        # TODO tidy this up

        self.f.canvas.callbacks.connect('button_press_event', self.update_window)
        self.f.canvas.callbacks.connect('scroll_event', self.update_window)

        self.flip_fov_button['command'] = self.flip_fov
        self.reset_fov_button['command'] = self.reset_fov
        self.reverse_color_button['command'] = self.reverse_color
        self.plot_white_entry.bind("<ButtonRelease-1>", self.update_window)
        self.plot_black_entry.bind("<ButtonRelease-1>", self.update_window)

        for num in range(len(self.targets_clear_button)):
            self.targets_clear_button[num]['command'] = self.clear_target(num)

        self.mirror_fov_button['command'] = self.mirror_fov
        self.photometry_button['command'] = self.photometry
        self.proceed_to_fitting_button['command'] = self.proceed_to_fitting
        self.return_to_reduction_button['command'] = self.return_to_reduction

        for target in range(self.max_targets):
            self.targets_aperture_entry[target].bind(sequence='<KeyRelease>', func=self.update_window)

        # TODO tidy this up

        # loop

        self.show_fov_window.show()
        self.loop()

    def update_window(self, event):

        if self.running.get():

            for i_target in range(self.max_targets):
                self.targets_indication_entry[i_target]['state'] = DISABLED
                self.targets_aperture_entry[i_target]['state'] = DISABLED

            self.photometry_button['state'] = DISABLED
            self.proceed_to_fitting_button['state'] = DISABLED
            self.return_to_reduction_button['state'] = DISABLED

        else:

            if isinstance(event, matplotlib.backend_bases.MouseEvent):

                if event.inaxes is None:
                    return None

                elif event.name == 'scroll_event':

                    zoom_factor = 1.2
                    scale_factor = 1.0

                    if event.button == 'up':
                        scale_factor = 1 / zoom_factor
                    elif event.button == 'down':
                        scale_factor = zoom_factor

                    xdata = event.xdata
                    ydata = event.ydata

                    cur_xlim = self.ax.get_xlim()
                    cur_ylim = self.ax.get_ylim()

                    cur_xrange = (cur_xlim[1] - cur_xlim[0])
                    cur_yrange = (cur_ylim[1] - cur_ylim[0])

                    new_xrange = cur_xrange * scale_factor
                    new_yrange = cur_yrange * scale_factor

                    new_xmin = xdata - new_xrange * (xdata - cur_xlim[0]) / cur_xrange
                    new_ymin = ydata - new_yrange * (ydata - cur_ylim[0]) / cur_yrange

                    self.ax.set_xlim([new_xmin, new_xmin + new_xrange])
                    self.ax.set_ylim([new_ymin, new_ymin + new_yrange])

                    new_circles_radius = max(0.03 * max(new_xrange, new_yrange), 5 * self.star_std)
                    new_circles_radius = min(new_circles_radius, self.circles_radius)
                    for nn in range(len(self.good_comps_boxes1)):
                        self.good_comps_boxes1[nn].set_radius(new_circles_radius)
                        self.good_comps_boxes2[nn].set_radius(0.75 * new_circles_radius)

                    self.canvas.draw()

                    return None

                elif event.dblclick:

                    star = plc.find_single_star(
                        self.fits[1].data, event.xdata, event.ydata,
                        mean=self.fits[1].header[self.mean_key], std=self.fits[1].header[self.std_key],
                        burn_limit=self.burn_limit * 0.95, star_std=self.star_std)

                    if not star:
                        showinfo('Star not acceptable.',
                                 'Star could not be located or it is close to saturation.')

                    elif (star[0] < self.visible_fov_x_min or star[0] > self.visible_fov_x_max
                          or star[1] < self.visible_fov_y_min
                          or star[1] > self.visible_fov_y_max):

                        showinfo('Star not acceptable.',
                                 'Star moves outside the FOV later.')

                    else:

                        self.targets_x_position[self.targets_indication.get()].set(round(star[0], 1))
                        self.targets_y_position[self.targets_indication.get()].set(round(star[1], 1))
                        if self.star_psf > 0:
                            self.targets_aperture[self.targets_indication.get()].set(round(3.3 * self.star_psf, 2))
                        else:
                            self.targets_aperture[self.targets_indication.get()].set(round(3.3 * max(star[4], star[5])))
                        self.targets_flux[self.targets_indication.get()].set(2 * np.pi * star[2] * star[4] * star[5])

                else:
                    return None

            if self.plot_white_entry.get() <= self.plot_black_entry.get():
                self.plot_white_entry.set(self.plot_black_entry.get() + 1)
            self.image.set_clim(self.plot_black_entry.get(), self.plot_white_entry.get())

            self.return_to_reduction_button['state'] = NORMAL
            self.photometry_button['state'] = NORMAL

            for i_target in range(self.max_targets):
                self.targets_indication_entry[i_target]['state'] = NORMAL

            try:

                for i_target in range(self.max_targets):

                    if 0 in [self.targets_x_position[i_target].get(), self.targets_y_position[i_target].get()]:

                        self.targets_box[i_target].set_center((-10000, -10000))

                        self.targets_text[i_target].set_x(-10000)
                        self.targets_text[i_target].set_y(-10000)

                        self.targets_aperture_entry[i_target]['state'] = DISABLED

                    else:

                        # for compatibility with older versions

                        if self.targets_flux[i_target].get() == 0:

                            star = plc.find_single_star(
                                self.fits[1].data, self.targets_x_position[i_target].get(),
                                self.targets_y_position[i_target].get(),
                                mean=self.fits[1].header[self.mean_key], std=self.fits[1].header[self.std_key],
                                burn_limit=self.burn_limit * 0.95, star_std=self.star_std)

                            if not star:
                                xx = self.clear_target(i_target)
                                xx()
                                return None

                            self.targets_flux[i_target].set(2 * np.pi * star[2] * star[4] * star[5])

                        try:
                            app = self.targets_aperture[i_target].get()
                            if self.star_psf == 0:
                                self.targets_warning[i_target].set('')
                            elif app < 2 * self.star_psf:
                                self.targets_warning[i_target].set('Apert. too small')
                            else:
                                self.targets_warning[i_target].set('')
                        except TclError:
                            self.targets_warning[i_target].set('')
                            app = 1

                        # for compatibility with older versions

                        if i_target == 0:

                            good_comps = []

                            for j, comp_star in enumerate(self.all_stars):
                                if np.sqrt((comp_star[0] - self.targets_x_position[i_target].get())**2 +
                                           (comp_star[1] - self.targets_y_position[i_target].get())**2) > self.star_std:
                                    if self.in_fov[j]:
                                        if comp_star[-1] < 1.4 * self.targets_flux[i_target].get():
                                            if comp_star[-1] > 0.6 * self.targets_flux[i_target].get():
                                                good_comps.append(comp_star)

                            good_comps = sorted(good_comps,
                                                key=lambda x: np.sqrt((x[0] -
                                                                       self.targets_x_position[i_target].get()) ** 2 +
                                                                      (x[1] -
                                                                       self.targets_y_position[i_target].get()) ** 2))

                            for num in range(self.max_comparisons):

                                if num < len(good_comps):
                                    self.good_comps_boxes1[num].set_center((good_comps[num][0], good_comps[num][1]))
                                    self.good_comps_boxes2[num].set_center((good_comps[num][0], good_comps[num][1]))
                                else:
                                    self.good_comps_boxes1[num].set_center((-1000, -1000))
                                    self.good_comps_boxes2[num].set_center((-1000, -1000))

                        if 0 in [self.targets_x_position[0].get(), self.targets_y_position[0].get()]:
                            for num in range(self.max_comparisons):
                                self.good_comps_boxes1[num].set_center((-1000, -1000))
                                self.good_comps_boxes2[num].set_center((-1000, -1000))

                        if i_target > 0:
                            if 0 not in [self.targets_x_position[0].get(), self.targets_y_position[0].get()]:
                                if self.targets_flux[i_target].get() > 2 * self.targets_flux[0].get():
                                    self.targets_warning[i_target].set(
                                        self.targets_warning[i_target].get() + ' Comp. too bright')
                                elif self.targets_flux[i_target].get() < 0.5 * self.targets_flux[0].get():
                                    self.targets_warning[i_target].set(
                                        self.targets_warning[i_target].get() + ' Comp. too faint')

                        self.targets_box[i_target].set_center((self.targets_x_position[i_target].get(),
                                                               self.targets_y_position[i_target].get()))

                        self.targets_box[i_target].set_radius(app)

                        self.targets_text[i_target].set_x(self.targets_x_position[i_target].get() + app + 1)
                        self.targets_text[i_target].set_y(self.targets_y_position[i_target].get() - app - 1)

                        y1 = int(self.targets_y_position[i_target].get() - app)
                        y2 = y1 + int(2 * app) + 2
                        x1 = int(self.targets_x_position[i_target].get() - app)
                        x2 = x1 + int(2 * app) + 2
                        self.targets_peak_counts[i_target].set(int(np.max(self.fits[1].data[y1:y2, x1:x2])))

                        self.targets_aperture_entry[i_target]['state'] = NORMAL

            except ValueError:
                self.photometry_button['state'] = DISABLED

            for i_target in range(self.max_targets):
                if 0 not in [self.targets_x_position[i_target].get(), self.targets_y_position[i_target].get()]:
                    try:
                        app = self.targets_aperture[i_target].get()
                    except TclError:
                        self.photometry_button['state'] = DISABLED

            if 0 in [self.targets_x_position[0].get(), self.targets_y_position[0].get(),
                     self.targets_x_position[1].get(), self.targets_y_position[1].get()]:
                self.photometry_button['state'] = DISABLED

            if (log.read_local_log('pipeline', 'photometry_complete')
               and len(glob.glob(os.path.join('{0}*'.format(self.photometry_directory),
                                              self.light_curve_aperture_file))) > 0):
                self.proceed_to_fitting_button['state'] = NORMAL

            else:
                self.proceed_to_fitting_button['state'] = DISABLED

        self.canvas.draw()

    # define actions for the different buttons, including calls to the function that updates the window

    def flip_fov(self):
        self.ax.set_ylim(self.ax.get_ylim()[1], self.ax.get_ylim()[0])
        self.canvas.draw()

    def mirror_fov(self):
        self.ax.set_xlim(self.ax.get_xlim()[1], self.ax.get_xlim()[0])
        self.canvas.draw()

    def reverse_color(self):
        if self.image.cmap == Greys_r:
            self.image.cmap = Greys
            self.min_label.set('Min (white pixels)')
            self.max_label.set('Max (black pixels)')
        else:
            self.image.cmap = Greys_r
            self.min_label.set('Min (black pixels)')
            self.max_label.set('Max (white pixels)')
        self.canvas.draw()

    def reset_fov(self):
        self.plot_black_entry.set(self.fits[1].header[self.mean_key] +
                                  self.frame_low_std * self.fits[1].header[self.std_key])
        self.plot_white_entry.set(self.fits[1].header[self.mean_key] +
                                  self.frame_upper_std * self.fits[1].header[self.std_key])

        self.plot_black.set(self.fits[1].header[self.mean_key] +
                            self.frame_low_std * self.fits[1].header[self.std_key])
        self.plot_white.set(self.fits[1].header[self.mean_key] +
                            self.frame_upper_std * self.fits[1].header[self.std_key])

        self.ax.set_xlim(0, len(self.fits[1].data[0]))
        self.ax.set_ylim(0, len(self.fits[1].data))
        self.image.set_clim(self.plot_black_entry.get(), self.plot_white_entry.get())

        new_circles_radius = max(0.03 * max(len(self.fits[1].data[0]), len(self.fits[1].data)), 5 * self.star_std)
        for nn in range(len(self.good_comps_boxes1)):
            self.good_comps_boxes1[nn].set_radius(new_circles_radius)
            self.good_comps_boxes2[nn].set_radius(0.75 * new_circles_radius)

        self.canvas.draw()

    def clear_target(self, num):

        def clear_target_num():

            self.targets_x_position[num].set(0.0)
            self.targets_y_position[num].set(0.0)
            self.targets_aperture[num].set(0.0)
            self.targets_peak_counts[num].set(0)
            self.targets_flux[num].set(0.0)
            self.targets_warning[num].set('')

            for i_target in range(1, self.max_targets - 1):
                if self.targets_x_position[i_target].get() == 0 and self.targets_x_position[i_target + 1].get() != 0:
                    self.targets_x_position[i_target].set(self.targets_x_position[i_target + 1].get())
                    self.targets_y_position[i_target].set(self.targets_y_position[i_target + 1].get())
                    self.targets_aperture[i_target].set(self.targets_aperture[i_target + 1].get())
                    self.targets_peak_counts[i_target].set(self.targets_peak_counts[i_target + 1].get())
                    self.targets_flux[i_target].set(self.targets_flux[i_target + 1].get())
                    self.targets_warning[i_target].set(self.targets_warning[i_target + 1].get())

                    self.targets_x_position[i_target + 1].set(0)
                    self.targets_y_position[i_target + 1].set(0)
                    self.targets_aperture[i_target + 1].set(0)
                    self.targets_peak_counts[i_target + 1].set(0)
                    self.targets_flux[i_target + 1].set(0)
                    self.targets_warning[i_target + 1].set('')

            self.update_window(None)

        return clear_target_num

    def photometry(self):

        self.running.set(True)
        self.update_window(None)

        log.write_local_log('photometry', self.targets_x_position[0].get(), 'target_x_position')
        log.write_local_log('photometry', self.targets_y_position[0].get(), 'target_y_position')
        log.write_local_log('photometry', self.targets_aperture[0].get(), 'target_aperture')
        target_polar = plc.cartesian_to_polar(self.targets_x_position[0].get(), self.targets_y_position[0].get(),
                                              self.fits[1].header[self.align_x0_key],
                                              self.fits[1].header[self.align_y0_key])
        log.write_local_log('photometry', float(target_polar[0]), 'target_r_position')
        log.write_local_log('photometry', float(target_polar[1]), 'target_u_position')

        for i_comparison in range(self.max_comparisons):
            log.write_local_log('photometry', self.targets_x_position[i_comparison + 1].get(),
                                'comparison_{0}_x_position'.format(i_comparison + 1))
            log.write_local_log('photometry', self.targets_y_position[i_comparison + 1].get(),
                                'comparison_{0}_y_position'.format(i_comparison + 1))
            log.write_local_log('photometry', self.targets_aperture[i_comparison + 1].get(),
                                'comparison_{0}_aperture'.format(i_comparison + 1))

            if 0 not in [self.targets_x_position[i_comparison + 1].get(),
                         self.targets_y_position[i_comparison + 1].get()]:

                target_polar = plc.cartesian_to_polar(self.targets_x_position[i_comparison + 1].get(),
                                                      self.targets_y_position[i_comparison + 1].get(),
                                                      self.fits[1].header[self.align_x0_key],
                                                      self.fits[1].header[self.align_y0_key])

            else:

                target_polar = [0, 0]

            log.write_local_log('photometry', float(target_polar[0]),
                                'comparison_{0}_r_position'.format(i_comparison + 1))
            log.write_local_log('photometry', float(target_polar[1]),
                                'comparison_{0}_u_position'.format(i_comparison + 1))
        ax_now = [self.ax.get_xlim(), self.ax.get_ylim()]
        self.ax.set_xlim(0, len(self.fits[1].data[0]))
        self.ax.set_ylim(0, len(self.fits[1].data))
        self.f.savefig(self.fov_figure, dpi=200)
        self.ax.set_xlim(ax_now[0][0], ax_now[0][1])
        self.ax.set_ylim(ax_now[1][0], ax_now[1][1])
        phr_photometry()

        self.running.set(False)
        self.update_window(None)

    def proceed_to_fitting(self):
        self.run.run_from_reduction = False
        self.run.run_from_photometry = False
        self.run.run_from_fitting = True
        self.show_fov_window.close()
        self.close()

    def return_to_reduction(self):
        self.run.run_from_reduction = True
        self.run.run_from_photometry = False
        self.run.run_from_fitting = False
        self.show_fov_window.close()
        self.close()


class FittingWindow(MainWindow):

    def __init__(self, run):

        self.run = run
        MainWindow.__init__(self, 'HOPS - Fitting', position=11)
        self.show_preview_window = SideWindow('FOV', position=10)

        # get variables from log and set as tk variables those to be modified

        self.observation_files = self.StringVar(log.read_local_log('pipeline', 'observation_files'))
        self.exposure_time_key = log.read_local_log('pipeline_keywords', 'exposure_time_key')
        self.light_curve_file = self.StringVar(log.read_local_log('fitting', 'light_curve_file'))
        self.light_curve_aperture_file = log.read_local_log('pipeline', 'light_curve_aperture_file')
        self.light_curve_gauss_file = log.read_local_log('pipeline', 'light_curve_gauss_file')
        self.photometry_directory = log.read_local_log('pipeline', 'photometry_directory')
        self.light_curve_file.set(glob.glob(os.path.join('{0}*'.format(self.photometry_directory), 
                                                         self.light_curve_aperture_file))[-1])
        self.observer = self.StringVar(log.read_local_log('fitting', 'observer'))
        self.observatory = self.StringVar(log.read_local_log('fitting', 'observatory'))
        self.telescope = self.StringVar(log.read_local_log('fitting', 'telescope'))
        self.camera = self.StringVar(log.read_local_log('fitting', 'camera'))
        self.phot_filter = self.StringVar(log.read_local_log('fitting', 'phot_filter'))
        self.scatter = self.DoubleVar(log.read_local_log('fitting', 'scatter'))
        self.iterations = self.IntVar(log.read_local_log('fitting', 'iterations'))
        self.burn = self.IntVar(log.read_local_log('fitting', 'burn'))
        self.manual_planet = self.BooleanVar(log.read_log('fitting', 'manual_planet'))
        self.planet = self.StringVar(log.read_local_log('fitting', 'planet'))
        self.target_ra_dec = self.StringVar(log.read_local_log('fitting', 'target_ra_dec'))
        self.metallicity = self.DoubleVar(log.read_local_log('fitting', 'metallicity'))
        self.temperature = self.DoubleVar(log.read_local_log('fitting', 'temperature'))
        self.logg = self.DoubleVar(log.read_local_log('fitting', 'logg'))
        self.period = self.DoubleVar(log.read_local_log('fitting', 'period'))
        self.mid_time = self.DoubleVar(log.read_local_log('fitting', 'mid_time'))
        self.rp_over_rs = self.DoubleVar(log.read_local_log('fitting', 'rp_over_rs'))
        self.sma_over_rs = self.DoubleVar(log.read_local_log('fitting', 'sma_over_rs'))
        self.inclination = self.DoubleVar(log.read_local_log('fitting', 'inclination'))
        self.eccentricity = self.DoubleVar(log.read_local_log('fitting', 'eccentricity'))
        self.periastron = self.DoubleVar(log.read_local_log('fitting', 'periastron'))
        
        ra_target, dec_target = log.read_local_log('photometry', 'target_ra_dec').split(' ')        
        ecc_planet = plc.find_nearest(plc.Target(plc.Hours(ra_target), plc.Degrees(dec_target)))
        
        self.auto_planet = self.StringVar(ecc_planet.planet.name)
        self.auto_target_ra_dec = self.StringVar('{0} {1}'.format(ecc_planet.star.ra, ecc_planet.star.dec))
        self.auto_metallicity = self.DoubleVar(ecc_planet.star.met)
        self.auto_temperature = self.DoubleVar(ecc_planet.star.teff)
        self.auto_logg = self.DoubleVar(ecc_planet.star.logg)
        self.auto_period = self.DoubleVar(ecc_planet.planet.period)
        self.auto_mid_time = self.DoubleVar(0)
        self.auto_rp_over_rs = self.DoubleVar(ecc_planet.planet.rp_over_rs)
        self.auto_sma_over_rs = self.DoubleVar(ecc_planet.planet.sma_over_rs)
        self.auto_eccentricity = self.DoubleVar(ecc_planet.planet.eccentricity)
        self.auto_inclination = self.DoubleVar(ecc_planet.planet.inclination)
        self.auto_periastron = self.DoubleVar(ecc_planet.planet.periastron)
        target = plc.Target(plc.Hours(ecc_planet.star.ra), plc.Degrees(ecc_planet.star.dec))
        if ecc_planet.planet.time_format in ['BJD_TDB', 'BJD_TT']:
            self.auto_mid_time.set(ecc_planet.planet.mid_time)
        elif ecc_planet.planet.time_format == 'BJD_UTC':
            self.auto_mid_time.set(plc.BJDUTC(ecc_planet.planet.mid_time, target).bjd_tdb(target))
        elif ecc_planet.planet.time_format in ['HJD_TDB', 'HJD_TT']:
            self.auto_mid_time.set(plc.HJDTDB(ecc_planet.planet.mid_time, target).bjd_tdb(target))
        elif ecc_planet.planet.time_format == 'HJD_UTC':
            self.auto_mid_time.set(plc.HJDUTC(ecc_planet.planet.mid_time, target).bjd_tdb(target))
        elif ecc_planet.planet.time_format == 'JD_UTC':
            self.auto_mid_time.set(plc.JD(ecc_planet.planet.mid_time).bjd_tdb(target))

        if self.planet.get() == 'Choose Planet':
            self.planet.set(self.auto_planet.get())
            self.target_ra_dec.set(self.auto_target_ra_dec.get())
            self.metallicity.set(self.auto_metallicity.get())
            self.temperature.set(self.auto_temperature.get())
            self.logg.set(self.auto_logg.get())
            self.period.set(self.auto_period.get())
            self.mid_time.set(self.auto_mid_time.get())
            self.rp_over_rs.set(self.auto_rp_over_rs.get())
            self.sma_over_rs.set(self.auto_sma_over_rs.get())
            self.inclination.set(self.auto_inclination.get())
            self.eccentricity.set(self.auto_eccentricity.get())
            self.periastron.set(self.auto_periastron.get())

        # set progress variables, useful for updating the window

        self.refit = self.BooleanVar(True)

        # create widgets
        combostyle = self.Style()
        combostyle.theme_create('combostyle', parent='alt',
                                settings={'TCombobox': {'configure':
                                                        {'selectbackground': 'white',
                                                         'fieldbackground': 'white',
                                                         'background': 'white'}}})
        combostyle.theme_use('combostyle')

        self.observer_label = self.Label(text='Observer')
        self.observer_entry = self.Entry(textvariable=self.observer)

        self.observatory_label = self.Label(text='Observatory')
        self.observatory_entry = self.Entry(textvariable=self.observatory)

        self.telescope_label = self.Label(text='Telescope')
        self.telescope_entry = self.Entry(textvariable=self.telescope)

        self.camera_label = self.Label(text='Camera')
        self.camera_entry = self.Entry(textvariable=self.camera)

        self.phot_filter_label = self.Label(text='Filter')
        self.phot_filter_entry = self.Combobox(textvariable=self.phot_filter, state='readonly', width=17)
        self.phot_filter_entry['values'] = tuple([ff for ff in filter_map])

        self.light_curve_file_label = self.Label(text='Light-curve file')
        self.light_curve_file_entry = self.Combobox(textvariable=self.light_curve_file, state='readonly', width=55)

        self.scatter_label = self.Label(text='Scatter limit')
        self.scatter_entry = self.Entry(textvariable=self.scatter, validate='key')
        self.scatter_entry['validatecommand'] = (self.scatter_entry.register(test_float_positive_input), '%P', '%d')

        self.iterations_label = self.Label(text='MCMC Iterations')
        self.iterations_entry = self.Entry(textvariable=self.iterations, validate='key')
        self.iterations_entry['validatecommand'] = (self.iterations_entry.register(test_int_positive_non_zero_input), '%P', '%d')

        self.burn_label = self.Label(text='MCMC Burn-in')
        self.burn_entry = self.Entry(textvariable=self.burn, validate='key')
        self.burn_entry['validatecommand'] = (self.burn_entry.register(test_int_positive_non_zero_input), '%P', '%d')

        self.auto_params_entry = self.Radiobutton(text='Use detected planet param.', variable=self.manual_planet, value=False)
        self.manual_params_entry = self.Radiobutton(text='Enter param. manually', variable=self.manual_planet, value=True)

        self.planet_label = self.Label(text='Planet')
        self.auto_planet_entry = self.Entry(textvariable=self.auto_planet, state=DISABLED)
        self.planet_entry = self.Entry(textvariable=self.planet)

        self.target_ra_dec_label = self.Label(text='Planet RA DEC\n(hh:mm:ss +/-dd:mm:ss)')
        self.auto_target_ra_dec_entry = self.Entry(textvariable=self.auto_target_ra_dec, state=DISABLED)
        self.target_ra_dec_entry = self.Entry(textvariable=self.target_ra_dec)
        self.target_ra_dec_test = self.Label(text=' ')

        self.metallicity_label = self.Label(text='M* [Fe/H, dex]')
        self.auto_metallicity_entry = self.Entry(textvariable=self.auto_metallicity, state=DISABLED)
        self.metallicity_entry = self.Entry(textvariable=self.metallicity, validate='key')
        self.metallicity_entry['validatecommand'] = (self.metallicity_entry.register(test_float_input), '%P', '%d')

        self.temperature_label = self.Label(text='T* [K]')
        self.auto_temperature_entry = self.Entry(textvariable=self.auto_temperature, state=DISABLED)
        self.temperature_entry = self.Entry(textvariable=self.temperature, validate='key')
        self.temperature_entry['validatecommand'] = (self.temperature_entry.register(test_float_positive_input), '%P', '%d')

        self.logg_label = self.Label(text='log(g*) [cm/s^2]')
        self.auto_logg_entry = self.Entry(textvariable=self.auto_logg, state=DISABLED)
        self.logg_entry = self.Entry(textvariable=self.logg, validate='key')
        self.logg_entry['validatecommand'] = (self.logg_entry.register(test_float_positive_input), '%P', '%d')

        self.period_label = self.Label(text='Period [days]')
        self.auto_period_entry = self.Entry(textvariable=self.auto_period, state=DISABLED)
        self.period_entry = self.Entry(textvariable=self.period, validate='key')
        self.period_entry['validatecommand'] = (self.period_entry.register(test_float_positive_input), '%P', '%d')

        self.mid_time_label = self.Label(text='Mid-time [BJD_TDB]')
        self.auto_mid_time_entry = self.Entry(textvariable=self.auto_mid_time, state=DISABLED)
        self.mid_time_entry = self.Entry(textvariable=self.mid_time, validate='key')
        self.mid_time_entry['validatecommand'] = (self.mid_time_entry.register(test_float_positive_input), '%P', '%d')

        self.rp_over_rs_label = self.Label(text='Rp/Rs')
        self.auto_rp_over_rs_entry = self.Entry(textvariable=self.auto_rp_over_rs, state=DISABLED)
        self.rp_over_rs_entry = self.Entry(textvariable=self.rp_over_rs, validate='key')
        self.rp_over_rs_entry['validatecommand'] = (self.rp_over_rs_entry.register(test_float_positive_input), '%P', '%d')

        self.sma_over_rs_label = self.Label(text='a/Rs')
        self.auto_sma_over_rs_entry = self.Entry(textvariable=self.auto_sma_over_rs, state=DISABLED)
        self.sma_over_rs_entry = self.Entry(textvariable=self.sma_over_rs, validate='key')
        self.sma_over_rs_entry['validatecommand'] = (self.sma_over_rs_entry.register(test_float_positive_input), '%P', '%d')

        self.inclination_label = self.Label(text='Inclination [deg]')
        self.auto_inclination_entry = self.Entry(textvariable=self.auto_inclination, state=DISABLED)
        self.inclination_entry = self.Entry(textvariable=self.inclination, validate='key')
        self.inclination_entry['validatecommand'] = (self.inclination_entry.register(test_float_positive_input), '%P', '%d')

        self.eccentricity_label = self.Label(text='Eccentricity')
        self.auto_eccentricity_entry = self.Entry(textvariable=self.auto_eccentricity, state=DISABLED)
        self.eccentricity_entry = self.Entry(textvariable=self.eccentricity, validate='key')
        self.eccentricity_entry['validatecommand'] = (self.eccentricity_entry.register(test_float_positive_input), '%P', '%d')

        self.periastron_label = self.Label(text='Periastron [deg]')
        self.auto_periastron_entry = self.Entry(textvariable=self.auto_periastron, state=DISABLED)
        self.periastron_entry = self.Entry(textvariable=self.periastron, validate='key')
        self.periastron_entry['validatecommand'] = (self.periastron_entry.register(test_float_positive_input), '%P', '%d')

        self.show_preview_button = self.Button(text='Show Preview')

        self.return_to_photometry_button = self.Button(text='RETURN TO PHOTOMETRY')

        self.return_to_reduction_button = self.Button(text='RETURN TO REDUCTION')

        self.fitting_button = self.Button(text='RUN FITTING', bg='green')

        self.my_profile_button = self.Button(text='MY PROFILE')

        # define additional windows

        funit = 1.0
        fcol = 7
        frow = 5
        fbottom = 0.11
        fright = 0.05
        self.fsmain = 10
        self.fsbig = 15
        fig = Figure(figsize=(funit * fcol / (1 - fright), funit * frow / (1 - fbottom)))
        self.canvas = self.show_preview_window.FigureCanvasTkAgg(fig)
        self.canvas.get_tk_widget().pack()
        self.show_preview_window.NavigationToolbar2Tk(self.canvas)
        try:
            gs = gridspec.GridSpec(frow, fcol, fig, 0, fbottom, 1.0 - fright, 1.0, 0.0, 0.0)
        except TypeError:
            gs = gridspec.GridSpec(frow, fcol, 0, fbottom, 1.0 - fright, 1.0, 0.0, 0.0)

        logo_ax = fig.add_subplot(gs[0, 0])
        logo_ax.imshow(log.holomon_logo_jpg)
        logo_ax.spines['top'].set_visible(False)
        logo_ax.spines['bottom'].set_visible(False)
        logo_ax.spines['left'].set_visible(False)
        logo_ax.spines['right'].set_visible(False)
        logo_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

        self.ax1 = fig.add_subplot(gs[1:4, 1:])
        self.ax2 = fig.add_subplot(gs[4, 1:])
        fig.text(0.04, fbottom + 2.5 * (1 - fbottom) / frow, 'relative flux (de-trended)', fontsize=self.fsbig, va='center',
                 ha='center', rotation='vertical')
        fig.text(0.04, fbottom + 0.5 * (1 - fbottom) / frow, 'residuals', fontsize=self.fsbig, va='center',
                 ha='center', rotation='vertical')

        self.title1 = fig.text(0.5, 0.94, '', fontsize=24, va='center', ha='center')
        self.title2 = fig.text(0.97, 0.97, '', fontsize=self.fsmain, va='top', ha='right')
        self.title3 = fig.text((1 - fright) / fcol, 1 - (1 - fbottom) / frow, '', fontsize=self.fsmain, ha='left', va='bottom')

        # connect widgets to functions

        self.light_curve_file_entry.bind('<<ComboboxSelected>>', self.update_window)
        self.planet_entry.bind(sequence='<KeyRelease>', func=self.update_window)
        self.auto_params_entry['command'] = self.update_window
        self.manual_params_entry['command'] = self.update_window
        self.camera_entry.bind(sequence='<KeyRelease>', func=self.update_window_no_refit)
        self.telescope_entry.bind(sequence='<KeyRelease>', func=self.update_window_no_refit)
        self.observatory_entry.bind(sequence='<KeyRelease>', func=self.update_window_no_refit)
        self.observer_entry.bind(sequence='<KeyRelease>', func=self.update_window_no_refit)
        self.scatter_entry.bind(sequence='<KeyRelease>', func=self.update_window)
        self.iterations_entry.bind(sequence='<KeyRelease>', func=self.update_window_no_refit)
        self.burn_entry.bind(sequence='<KeyRelease>', func=self.update_window_no_refit)
        self.phot_filter_entry.bind('<<ComboboxSelected>>', self.update_window)

        self.metallicity_entry.bind(sequence='<KeyRelease>', func=self.update_window_no_refit_if_auto)
        self.temperature_entry.bind(sequence='<KeyRelease>', func=self.update_window_no_refit_if_auto)
        self.logg_entry.bind(sequence='<KeyRelease>', func=self.update_window_no_refit_if_auto)
        self.period_entry.bind(sequence='<KeyRelease>', func=self.update_window_no_refit_if_auto)
        self.mid_time_entry.bind(sequence='<KeyRelease>', func=self.update_window_no_refit_if_auto)
        self.rp_over_rs_entry.bind(sequence='<KeyRelease>', func=self.update_window_no_refit_if_auto)
        self.sma_over_rs_entry.bind(sequence='<KeyRelease>', func=self.update_window_no_refit_if_auto)
        self.inclination_entry.bind(sequence='<KeyRelease>', func=self.update_window_no_refit_if_auto)
        self.eccentricity_entry.bind(sequence='<KeyRelease>', func=self.update_window_no_refit_if_auto)
        self.periastron_entry.bind(sequence='<KeyRelease>', func=self.update_window_no_refit_if_auto)
        self.target_ra_dec_entry.bind(sequence='<KeyRelease>', func=self.update_window_no_refit_if_auto)
        self.show_preview_button['command'] = self.show_preview_window.show
        self.return_to_photometry_button['command'] = self.return_to_photometry
        self.return_to_reduction_button['command'] = self.return_to_reduction
        self.fitting_button['command'] = self.fitting
        self.my_profile_button['command'] = self.my_profile_window.show

        # setup window

        self.setup_window([
            [[self.get_logo(), 0, 1, 8]],
            [],
            [[self.Label(text='Fitting'), 1, 5, 1, 'title']],
            [],
            [[self.light_curve_file_label, 1], [self.light_curve_file_entry, 2, 4, 1]],
            [],
            [[self.auto_params_entry, 4], [self.manual_params_entry, 5]],
            [[self.planet_label, 3], [self.auto_planet_entry, 4], [self.planet_entry, 5]],
            [[self.created_by_label, 0, 1, 3], [self.scatter_label, 1], [self.scatter_entry, 2], [self.target_ra_dec_label, 3],
             [self.auto_target_ra_dec_entry, 4], [self.target_ra_dec_entry, 5]],
            [[self.target_ra_dec_test, 5]],
            [[self.iterations_label, 1], [self.iterations_entry, 2], [self.period_label, 3], [self.auto_period_entry, 4],
             [self.period_entry, 5]],
            [[self.my_profile_button, 0], [self.burn_label, 1], [self.burn_entry, 2], [self.mid_time_label, 3], [self.auto_mid_time_entry, 4],
             [self.mid_time_entry, 5]],
            [[self.rp_over_rs_label, 3], [self.auto_rp_over_rs_entry, 4], [self.rp_over_rs_entry, 5]],
            [[self.Btn, 0, 1, 2], [self.phot_filter_label, 1], [self.phot_filter_entry, 2], [self.sma_over_rs_label, 3],
             [self.auto_sma_over_rs_entry, 4], [self.sma_over_rs_entry, 5]],
            [[self.camera_label, 1], [self.camera_entry, 2], [self.inclination_label, 3], [self.auto_inclination_entry, 4],
             [self.inclination_entry, 5]],
            [[self.telescope_label, 1], [self.telescope_entry, 2], [self.eccentricity_label, 3], [self.auto_eccentricity_entry, 4],
             [self.eccentricity_entry, 5]],
            [[self.observatory_label, 1], [self.observatory_entry, 2], [self.periastron_label, 3], [self.auto_periastron_entry, 4],
             [self.periastron_entry, 5]],
            [[self.observer_label, 1], [self.observer_entry, 2], [self.metallicity_label, 3], [self.auto_metallicity_entry, 4],
             [self.metallicity_entry, 5]],
            [[self.temperature_label, 3], [self.auto_temperature_entry, 4], [self.temperature_entry, 5]],
            [[self.logg_label, 3], [self.auto_logg_entry, 4], [self.logg_entry, 5]],
            [[self.fitting_button, 1, 4]],
            [[self.return_to_photometry_button, 1, 2], [self.return_to_reduction_button, 3, 2]],
            []
        ])

        self.update_window(None)

        # finalise and show  window

        self.show_preview_window.show()
        self.loop()

    def update_window(self, *entry):

        if not entry:
            pass

        self.phot_filter_entry.selection_clear()
        self.light_curve_file_entry.selection_clear()

        all_files = (glob.glob(os.path.join('{0}*'.format(self.photometry_directory), self.light_curve_aperture_file)) +
                     glob.glob(os.path.join('{0}*'.format(self.photometry_directory), self.light_curve_gauss_file)))
        all_files.sort()
        self.light_curve_file_entry['values'] = tuple(all_files)

        if self.running.get():

            self.light_curve_file_entry['state'] = DISABLED
            self.scatter_entry['state'] = DISABLED
            self.iterations_entry['state'] = DISABLED
            self.burn_entry['state'] = DISABLED
            self.metallicity_entry['state'] = DISABLED
            self.temperature_entry['state'] = DISABLED
            self.logg_entry['state'] = DISABLED
            self.phot_filter_entry['state'] = DISABLED
            self.period_entry['state'] = DISABLED
            self.mid_time_entry['state'] = DISABLED
            self.rp_over_rs_entry['state'] = DISABLED
            self.sma_over_rs_entry['state'] = DISABLED
            self.inclination_entry['state'] = DISABLED
            self.eccentricity_entry['state'] = DISABLED
            self.periastron_entry['state'] = DISABLED
            self.target_ra_dec_entry['state'] = DISABLED
            self.observer_entry['state'] = DISABLED
            self.observatory_entry['state'] = DISABLED
            self.telescope_entry['state'] = DISABLED
            self.camera_entry['state'] = DISABLED
            self.planet_entry['state'] = DISABLED
            self.manual_params_entry['state'] = DISABLED
            self.return_to_photometry_button['state'] = DISABLED
            self.return_to_reduction_button['state'] = DISABLED
            self.fitting_button['state'] = DISABLED
            self.my_profile_button['state'] = DISABLED
            self.show_preview_button['state'] = DISABLED

        elif not os.path.isfile(self.light_curve_file.get()):

            self.light_curve_file_entry['state'] = NORMAL
            self.scatter_entry['state'] = DISABLED
            self.iterations_entry['state'] = DISABLED
            self.burn_entry['state'] = DISABLED
            self.metallicity_entry['state'] = DISABLED
            self.temperature_entry['state'] = DISABLED
            self.logg_entry['state'] = DISABLED
            self.phot_filter_entry['state'] = DISABLED
            self.period_entry['state'] = DISABLED
            self.mid_time_entry['state'] = DISABLED
            self.rp_over_rs_entry['state'] = DISABLED
            self.sma_over_rs_entry['state'] = DISABLED
            self.inclination_entry['state'] = DISABLED
            self.eccentricity_entry['state'] = DISABLED
            self.periastron_entry['state'] = DISABLED
            self.target_ra_dec_entry['state'] = DISABLED
            self.observer_entry['state'] = DISABLED
            self.observatory_entry['state'] = DISABLED
            self.telescope_entry['state'] = DISABLED
            self.camera_entry['state'] = DISABLED
            self.planet_entry['state'] = DISABLED
            self.manual_params_entry['state'] = DISABLED
            self.return_to_photometry_button['state'] = DISABLED
            self.return_to_reduction_button['state'] = DISABLED
            self.fitting_button['state'] = DISABLED
            self.my_profile_button['state'] = NORMAL
            self.show_preview_button['state'] = DISABLED

        else:

            self.light_curve_file_entry['state'] = 'readonly'
            self.scatter_entry['state'] = NORMAL
            self.iterations_entry['state'] = NORMAL
            self.burn_entry['state'] = NORMAL
            self.phot_filter_entry['state'] = 'readonly'
            self.observer_entry['state'] = NORMAL
            self.observatory_entry['state'] = NORMAL
            self.telescope_entry['state'] = NORMAL
            self.camera_entry['state'] = NORMAL
            self.manual_params_entry['state'] = NORMAL
            self.my_profile_button['state'] = NORMAL
            self.show_preview_button['state'] = NORMAL
            self.planet_entry['state'] = NORMAL
            self.target_ra_dec_entry['state'] = NORMAL
            self.metallicity_entry['state'] = NORMAL
            self.temperature_entry['state'] = NORMAL
            self.logg_entry['state'] = NORMAL
            self.period_entry['state'] = NORMAL
            self.mid_time_entry['state'] = NORMAL
            self.rp_over_rs_entry['state'] = NORMAL
            self.sma_over_rs_entry['state'] = NORMAL
            self.inclination_entry['state'] = NORMAL
            self.eccentricity_entry['state'] = NORMAL
            self.periastron_entry['state'] = NORMAL
            self.target_ra_dec_entry['state'] = NORMAL

            if self.phot_filter.get() == 'default':
                for key in log.read_local_log_profile('filter_key').split(','):
                    check_filter = test_fits_keyword(self.observation_files.get(), key)
                    if check_filter[0]:
                        if check_filter[2] in filter_map:
                            self.phot_filter.set(check_filter[2])
                            break
                if self.phot_filter.get() == 'default':
                    if log.read_local_log_profile('filter') in filter_map:
                        self.phot_filter.set(log.read_local_log_profile('filter'))
                    else:
                        self.phot_filter.set(' ')

            if self.telescope.get() == 'default':
                for key in log.read_local_log_profile('telescope_key').split(','):
                    check_telescope = test_fits_keyword(self.observation_files.get(), key)
                    if check_telescope[0]:
                        self.telescope.set(check_telescope[2])
                        break
                if self.telescope.get() == 'default':
                    self.telescope.set(log.read_local_log_profile('telescope'))

            if self.camera.get() == 'default':
                for key in log.read_local_log_profile('camera_key').split(','):
                    check_camera = test_fits_keyword(self.observation_files.get(), key)
                    if check_camera[0]:
                        self.camera.set(check_camera[2])
                        break
                if self.camera.get() == 'default':
                    self.camera.set(log.read_local_log_profile('camera'))

            if self.observer.get() == 'default':
                for key in log.read_local_log_profile('observer_key').split(','):
                    check_observer = test_fits_keyword(self.observation_files.get(), key)
                    if check_observer[0]:
                        self.observer.set(check_observer[2])
                        break
                if self.observer.get() == 'default':
                    self.observer.set(log.read_local_log_profile('observer'))

            if self.observatory.get() == 'default':
                for key in log.read_local_log_profile('observatory_key').split(','):
                    check_observatory = test_fits_keyword(self.observation_files.get(), key)
                    if check_observatory[0]:
                        self.observatory.set(check_observatory[2])
                        break
                if self.observatory.get() == 'default':
                    self.observatory.set(log.read_local_log_profile('observatory'))

            enable_buttons = True

            if not os.path.isfile(self.light_curve_file.get()):
                enable_buttons = False

            check_ra_dec = test_coordinates(self.target_ra_dec_entry.get(), single_line=True)
            self.target_ra_dec_test.configure(text=check_ra_dec[1])

            if not check_ra_dec[0]:
                enable_buttons = False

            if self.phot_filter.get() not in filter_map:
                enable_buttons = False

            for input_entry in [self.scatter_entry, self.iterations_entry, self.burn_entry,
                                self.metallicity_entry, self.temperature_entry, self.logg_entry, self.period_entry, self.mid_time_entry,
                                self.rp_over_rs_entry, self.sma_over_rs_entry, self.inclination_entry, self.eccentricity_entry,
                                self.periastron_entry, self.target_ra_dec_entry]:

                if len(str(input_entry.get())) == 0:
                    enable_buttons = False

            if enable_buttons:
                self.fitting_button['state'] = NORMAL
                self.show_preview_button['state'] = NORMAL

            else:
                self.fitting_button['state'] = DISABLED
                self.show_preview_button['state'] = DISABLED

            self.return_to_photometry_button['state'] = NORMAL
            self.return_to_reduction_button['state'] = NORMAL

        if self.manual_planet.get():
            planet_to_plot = self.planet.get()
            target_ra_dec_to_plot = self.target_ra_dec.get()
            metallicity_to_plot = self.metallicity.get()
            temperature_to_plot = self.temperature.get()
            logg_to_plot = self.logg.get()
            period_to_plot = self.period.get()
            mid_time_to_plot = self.mid_time.get()
            rp_over_rs_to_plot = self.rp_over_rs.get()
            sma_over_rs_to_plot = self.sma_over_rs.get()
            inclination_to_plot = self.inclination.get()
            eccentricity_to_plot = self.eccentricity.get()
            periastron_to_plot = self.periastron.get()
        else:
            planet_to_plot = self.auto_planet.get()
            target_ra_dec_to_plot = self.auto_target_ra_dec.get()
            metallicity_to_plot = self.auto_metallicity.get()
            temperature_to_plot = self.auto_temperature.get()
            logg_to_plot = self.auto_logg.get()
            period_to_plot = self.auto_period.get()
            mid_time_to_plot = self.auto_mid_time.get()
            rp_over_rs_to_plot = self.auto_rp_over_rs.get()
            sma_over_rs_to_plot = self.auto_sma_over_rs.get()
            inclination_to_plot = self.auto_inclination.get()
            eccentricity_to_plot = self.auto_eccentricity.get()
            periastron_to_plot = self.auto_periastron.get()

        light_curve = np.loadtxt(self.light_curve_file.get(), unpack=True)

        date = plc.JD(light_curve[0][0]).utc.isoformat()[:16].replace('T', ' ')
        obs_duration = round(24 * (light_curve[0][-1] - light_curve[0][0]), 1)
        exp_time = test_fits_keyword(self.observation_files.get(), self.exposure_time_key)[2]

        self.title1.set_text('{0}{1}{2}'.format('$\mathbf{', planet_to_plot, '}$'))
        self.title2.set_text('{0} (UT)\nDur: {1}h / Exp: {2}s\nFilter: {3}'.format(date, obs_duration, exp_time,
                                                                                  self.phot_filter.get()))
        self.title3.set_text('\n\n{0}\n{1}'.format(
            self.observer.get(), '{0} / {1} / {2}'.format(self.observatory.get(), self.telescope.get(), self.camera.get())))

        if self.refit.get():

            self.ax1.cla()
            self.ax2.cla()
            try:

                # filter out outliers

                light_curve_0 = light_curve[0]
                light_curve_1 = light_curve[1]

                light_curve_0 = light_curve_0[np.where(~np.isnan(light_curve_1))]
                light_curve_1 = light_curve_1[np.where(~np.isnan(light_curve_1))]

                moving_average = []
                for i in range(-10, 11):
                    moving_average.append(np.roll(light_curve_1, i))

                median = np.median(moving_average, 0)
                med = np.median([np.abs(ff - median) for ff in moving_average], 0)

                flag = np.where((np.abs(light_curve_1 - median) < self.scatter.get() * med))[0]
                outliers = np.where((np.abs(light_curve_1 - median) >= self.scatter.get() * med))[0]

                light_curve_0_outliers = light_curve_0[outliers]
                light_curve_1_outliers = light_curve_1[outliers]
                light_curve_0 = light_curve_0[flag]
                light_curve_1 = light_curve_1[flag]

                # fix timing

                ra_dec_string = target_ra_dec_to_plot.replace(':', ' ').split(' ')
                target = plc.Target(plc.Hours(*ra_dec_string[:3]), plc.Degrees(*ra_dec_string[3:]))
                light_curve_0 = np.array([plc.JD(ff + 0.5 * exp_time / 60.0 / 60.0 / 24.0).bjd_tdb(target) for ff in light_curve_0])

                # predictions

                limb_darkening_coefficients = plc.clablimb('claret', logg_to_plot, max(4000, temperature_to_plot),
                                                           metallicity_to_plot, filter_map[self.phot_filter.get()])

                predicted_mid_time = (mid_time_to_plot +
                                      round((np.mean(light_curve_0) - mid_time_to_plot) / period_to_plot) * period_to_plot)

                predicted_transit_model = plc.transit_integrated('claret', limb_darkening_coefficients, rp_over_rs_to_plot,
                                                      period_to_plot, sma_over_rs_to_plot, eccentricity_to_plot,
                                                      inclination_to_plot, periastron_to_plot, predicted_mid_time,
                                                      light_curve_0, float(exp_time), max(1, int(float(exp_time) / 10)))

                # define models

                def mcmc_f(time_array, detrend_zero, detrend_one, detrend_two, model_rp_over_rs, model_mid_time):

                    data_delta_t = time_array - light_curve_0[0]

                    detrend = detrend_zero * (1 + detrend_one * data_delta_t +
                                              detrend_two * data_delta_t * data_delta_t)
                    transit_model = plc.transit_integrated('claret', limb_darkening_coefficients, model_rp_over_rs,
                                                period_to_plot, sma_over_rs_to_plot, eccentricity_to_plot,
                                                inclination_to_plot, periastron_to_plot,
                                                predicted_mid_time + model_mid_time,
                                                time_array, float(exp_time), max(1, int(float(exp_time) / 10)))

                    return detrend * transit_model

                def independent_f(time_array, detrend_zero, detrend_one, detrend_two, model_rp_over_rs, model_mid_time):

                    data_delta_t = time_array - light_curve_0[0]

                    detrend = detrend_zero * (1 + detrend_one * data_delta_t +
                                              detrend_two * data_delta_t * data_delta_t)
                    transit_model = plc.transit_integrated('claret', limb_darkening_coefficients, model_rp_over_rs, period_to_plot,
                                                sma_over_rs_to_plot, eccentricity_to_plot, inclination_to_plot,
                                                periastron_to_plot, predicted_mid_time + model_mid_time,
                                                time_array, float(exp_time), max(1, int(float(exp_time) / 10)))

                    return detrend, transit_model

                # set noise level

                sigma = np.array([np.roll(light_curve_1, ff) for ff in range(-10, 10)])
                sigma = np.std(sigma, 0)

                popt, pcov = curve_fit(mcmc_f, light_curve_0, light_curve_1,
                                       p0=[np.mean(light_curve_1), 1, -1, rp_over_rs_to_plot, 0],
                                       sigma=sigma, maxfev=10000)

                fit_detrend, fit_transit_model = independent_f(light_curve_0, *popt)

                test = []
                for i in range(-int(len(light_curve_0) / 2), int(len(light_curve_0) / 2)):
                    test.append([np.sum((light_curve_1 / fit_detrend - np.roll(fit_transit_model, i)) ** 2), i])
                test.sort()

                popt, pcov = curve_fit(mcmc_f, light_curve_0, light_curve_1,
                                       p0=[popt[0], popt[1], popt[2], popt[3],
                                           popt[4] + (test[0][1]) * exp_time / 60.0 / 60.0 / 24.0],
                                       sigma=sigma, maxfev=10000)

                residuals = light_curve_1 - mcmc_f(light_curve_0, *popt)

                sigma = np.array([np.roll(residuals, ff) for ff in range(-10, 10)])
                sigma = np.std(sigma, 0)

                # results

                popt, pcov = curve_fit(mcmc_f, light_curve_0, light_curve_1, p0=popt, sigma=sigma, maxfev=10000)
                new_mid_time = predicted_mid_time + popt[-1]

                fit_detrend, fit_transit_model = independent_f(light_curve_0, *popt)
                phase = np.array((light_curve_0 - new_mid_time) / period_to_plot)
                detrended_data = light_curve_1 / fit_detrend
                residuals = detrended_data - fit_transit_model
                std_res = np.std(residuals)
                res_autocorr = np.correlate(residuals, residuals, mode='full')
                res_autocorr = round(np.max(np.abs(
                    res_autocorr[res_autocorr.size // 2:] / res_autocorr[res_autocorr.size // 2:][0])[1:]), 2)

                fit_detrend_outliers, fit_transit_model_outliers = independent_f(light_curve_0_outliers, *popt)
                phase_outliers = np.array((light_curve_0_outliers - new_mid_time) / period_to_plot)
                detrended_data_outliers = light_curve_1_outliers / fit_detrend_outliers

                # plot

                self.ax1.plot(phase, detrended_data, 'ko', ms=2, label='De-trended data')
                self.ax1.plot(phase, fit_transit_model, 'r-', label='Best-fit model (T$_0$ = {0}, R$_p$/R$_s$ = {1})'.format(
                    round(new_mid_time, 5), round(popt[3], 5)))
                self.ax1.plot(phase, predicted_transit_model, 'c-', label='Expected model (T$_0$ = {0}, R$_p$/R$_s$ = {1})'.format(
                    round(predicted_mid_time, 5), round(rp_over_rs_to_plot, 5)))

                data_ymin = min(detrended_data) - 3 * std_res
                data_ymax = max(detrended_data) + 2 * std_res
                self.ax1.set_yticks(self.ax1.get_yticks()[np.where(self.ax1.get_yticks() > data_ymin)])
                ymin, ymax = data_ymax - 1.15 * (data_ymax - data_ymin), data_ymax
                self.ax1.set_ylim(ymin, ymax)
                x_max = max(np.abs(phase) + 0.05 * (max(phase) - min(phase)))

                self.ax1.set_xlim(-x_max, x_max)
                self.ax1.tick_params(labelbottom=False, labelsize=self.fsmain)

                self.ax1.plot(phase_outliers, detrended_data_outliers, 'ro', ms=2, label='Filtered outliers')
                self.ax1.legend(loc=3)

                self.ax1.text(x_max - 0.4*x_max, ymin + 0.1 * (ymax - ymin),
                         'O-C: {0} min'.format(round((new_mid_time - predicted_mid_time) * 24 * 60, 1)))

                self.ax2.plot(phase, residuals, 'ko', ms=2)
                self.ax2.plot(phase, np.zeros_like(phase), 'r-')

                self.ax2.set_ylim(- 5 * std_res, 5 * std_res)

                self.ax2.set_xlabel('phase', fontsize=self.fsbig)

                self.ax2.set_xlim(-x_max, x_max)
                self.ax2.tick_params(labelsize=self.fsmain)

                self.ax2.text(self.ax2.get_xlim()[0] + 0.02 * (self.ax2.get_xlim()[-1] - self.ax2.get_xlim()[0]),
                         self.ax2.get_ylim()[1] - 0.15 * (self.ax2.get_ylim()[-1] - self.ax2.get_ylim()[0]),
                         r'STD = %.1f $‰$' % round((std_res * 1000), 1),
                         fontsize=self.fsmain)
                self.ax2.text(self.ax2.get_xlim()[0] + 0.02 * (self.ax2.get_xlim()[-1] - self.ax2.get_xlim()[0]),
                         self.ax2.get_ylim()[0] + 0.07 * (self.ax2.get_ylim()[-1] - self.ax2.get_ylim()[0]),
                         r'AutoCorr = %.1f' %res_autocorr,
                         fontsize=self.fsmain)

            except:
                pass

        self.canvas.draw()

    # define actions for the different buttons, including calls to the function that updates the window

    def return_to_photometry(self):

        if self.manual_planet.get():
            planet_to_plot = self.planet.get()
            target_ra_dec_to_plot = self.target_ra_dec.get()
            metallicity_to_plot = self.metallicity.get()
            temperature_to_plot = self.temperature.get()
            logg_to_plot = self.logg.get()
            period_to_plot = self.period.get()
            mid_time_to_plot = self.mid_time.get()
            rp_over_rs_to_plot = self.rp_over_rs.get()
            sma_over_rs_to_plot = self.sma_over_rs.get()
            inclination_to_plot = self.inclination.get()
            eccentricity_to_plot = self.eccentricity.get()
            periastron_to_plot = self.periastron.get()
        else:
            planet_to_plot = self.auto_planet.get()
            target_ra_dec_to_plot = self.auto_target_ra_dec.get()
            metallicity_to_plot = self.auto_metallicity.get()
            temperature_to_plot = self.auto_temperature.get()
            logg_to_plot = self.auto_logg.get()
            period_to_plot = self.auto_period.get()
            mid_time_to_plot = self.auto_mid_time.get()
            rp_over_rs_to_plot = self.auto_rp_over_rs.get()
            sma_over_rs_to_plot = self.auto_sma_over_rs.get()
            inclination_to_plot = self.auto_inclination.get()
            eccentricity_to_plot = self.auto_eccentricity.get()
            periastron_to_plot = self.auto_periastron.get()

        log.write_local_log('fitting', self.light_curve_file.get(), 'light_curve_file')
        log.write_local_log('fitting', self.scatter.get(), 'scatter')
        log.write_local_log('fitting', self.iterations.get(), 'iterations')
        log.write_local_log('fitting', self.burn.get(), 'burn')
        log.write_local_log('fitting', self.observer.get(), 'observer')
        log.write_local_log('fitting', self.observatory.get(), 'observatory')
        log.write_local_log('fitting', self.telescope.get(), 'telescope')
        log.write_local_log('fitting', self.camera.get(), 'camera')
        log.write_local_log('fitting', self.phot_filter.get(), 'phot_filter')
        log.write_local_log('fitting', self.manual_planet.get(), 'manual_planet')
        log.write_local_log('fitting', planet_to_plot, 'planet')
        log.write_local_log('fitting', target_ra_dec_to_plot, 'target_ra_dec')
        log.write_local_log('fitting', metallicity_to_plot, 'metallicity')
        log.write_local_log('fitting', temperature_to_plot, 'temperature')
        log.write_local_log('fitting', logg_to_plot, 'logg')
        log.write_local_log('fitting', period_to_plot, 'period')
        log.write_local_log('fitting', mid_time_to_plot, 'mid_time')
        log.write_local_log('fitting', rp_over_rs_to_plot, 'rp_over_rs')
        log.write_local_log('fitting', sma_over_rs_to_plot, 'sma_over_rs')
        log.write_local_log('fitting', inclination_to_plot, 'inclination')
        log.write_local_log('fitting', eccentricity_to_plot, 'eccentricity')
        log.write_local_log('fitting', periastron_to_plot, 'periastron')

        self.run.run_from_reduction = False
        self.run.run_from_photometry = True
        self.run.run_from_fitting = False

        self.show_preview_window.close()
        self.my_profile_window.close()
        self.close()

    def return_to_reduction(self):

        if self.manual_planet.get():
            planet_to_plot = self.planet.get()
            target_ra_dec_to_plot = self.target_ra_dec.get()
            metallicity_to_plot = self.metallicity.get()
            temperature_to_plot = self.temperature.get()
            logg_to_plot = self.logg.get()
            period_to_plot = self.period.get()
            mid_time_to_plot = self.mid_time.get()
            rp_over_rs_to_plot = self.rp_over_rs.get()
            sma_over_rs_to_plot = self.sma_over_rs.get()
            inclination_to_plot = self.inclination.get()
            eccentricity_to_plot = self.eccentricity.get()
            periastron_to_plot = self.periastron.get()
        else:
            planet_to_plot = self.auto_planet.get()
            target_ra_dec_to_plot = self.auto_target_ra_dec.get()
            metallicity_to_plot = self.auto_metallicity.get()
            temperature_to_plot = self.auto_temperature.get()
            logg_to_plot = self.auto_logg.get()
            period_to_plot = self.auto_period.get()
            mid_time_to_plot = self.auto_mid_time.get()
            rp_over_rs_to_plot = self.auto_rp_over_rs.get()
            sma_over_rs_to_plot = self.auto_sma_over_rs.get()
            inclination_to_plot = self.auto_inclination.get()
            eccentricity_to_plot = self.auto_eccentricity.get()
            periastron_to_plot = self.auto_periastron.get()

        log.write_local_log('fitting', self.light_curve_file.get(), 'light_curve_file')
        log.write_local_log('fitting', self.scatter.get(), 'scatter')
        log.write_local_log('fitting', self.iterations.get(), 'iterations')
        log.write_local_log('fitting', self.burn.get(), 'burn')
        log.write_local_log('fitting', self.observer.get(), 'observer')
        log.write_local_log('fitting', self.observatory.get(), 'observatory')
        log.write_local_log('fitting', self.telescope.get(), 'telescope')
        log.write_local_log('fitting', self.camera.get(), 'camera')
        log.write_local_log('fitting', self.phot_filter.get(), 'phot_filter')
        log.write_local_log('fitting', self.manual_planet.get(), 'manual_planet')
        log.write_local_log('fitting', planet_to_plot, 'planet')
        log.write_local_log('fitting', target_ra_dec_to_plot, 'target_ra_dec')
        log.write_local_log('fitting', metallicity_to_plot, 'metallicity')
        log.write_local_log('fitting', temperature_to_plot, 'temperature')
        log.write_local_log('fitting', logg_to_plot, 'logg')
        log.write_local_log('fitting', period_to_plot, 'period')
        log.write_local_log('fitting', mid_time_to_plot, 'mid_time')
        log.write_local_log('fitting', rp_over_rs_to_plot, 'rp_over_rs')
        log.write_local_log('fitting', sma_over_rs_to_plot, 'sma_over_rs')
        log.write_local_log('fitting', inclination_to_plot, 'inclination')
        log.write_local_log('fitting', eccentricity_to_plot, 'eccentricity')
        log.write_local_log('fitting', periastron_to_plot, 'periastron')

        self.run.run_from_reduction = True
        self.run.run_from_photometry = False
        self.run.run_from_fitting = False
        self.show_preview_window.close()
        self.my_profile_window.close()
        self.close()

    def fitting(self):

        if self.manual_planet.get():
            planet_to_plot = self.planet.get()
            target_ra_dec_to_plot = self.target_ra_dec.get()
            metallicity_to_plot = self.metallicity.get()
            temperature_to_plot = self.temperature.get()
            logg_to_plot = self.logg.get()
            period_to_plot = self.period.get()
            mid_time_to_plot = self.mid_time.get()
            rp_over_rs_to_plot = self.rp_over_rs.get()
            sma_over_rs_to_plot = self.sma_over_rs.get()
            inclination_to_plot = self.inclination.get()
            eccentricity_to_plot = self.eccentricity.get()
            periastron_to_plot = self.periastron.get()
        else:
            planet_to_plot = self.auto_planet.get()
            target_ra_dec_to_plot = self.auto_target_ra_dec.get()
            metallicity_to_plot = self.auto_metallicity.get()
            temperature_to_plot = self.auto_temperature.get()
            logg_to_plot = self.auto_logg.get()
            period_to_plot = self.auto_period.get()
            mid_time_to_plot = self.auto_mid_time.get()
            rp_over_rs_to_plot = self.auto_rp_over_rs.get()
            sma_over_rs_to_plot = self.auto_sma_over_rs.get()
            inclination_to_plot = self.auto_inclination.get()
            eccentricity_to_plot = self.auto_eccentricity.get()
            periastron_to_plot = self.auto_periastron.get()

        self.running.set(True)
        self.update_window(None)

        log.write_local_log('fitting', self.light_curve_file.get(), 'light_curve_file')
        log.write_local_log('fitting', self.scatter.get(), 'scatter')
        log.write_local_log('fitting', self.iterations.get(), 'iterations')
        log.write_local_log('fitting', self.burn.get(), 'burn')
        log.write_local_log('fitting', self.observer.get(), 'observer')
        log.write_local_log('fitting', self.observatory.get(), 'observatory')
        log.write_local_log('fitting', self.telescope.get(), 'telescope')
        log.write_local_log('fitting', self.camera.get(), 'camera')
        log.write_local_log('fitting', self.phot_filter.get(), 'phot_filter')
        log.write_local_log('fitting', self.manual_planet.get(), 'manual_planet')
        log.write_local_log('fitting', planet_to_plot, 'planet')
        log.write_local_log('fitting', target_ra_dec_to_plot, 'target_ra_dec')
        log.write_local_log('fitting', metallicity_to_plot, 'metallicity')
        log.write_local_log('fitting', temperature_to_plot, 'temperature')
        log.write_local_log('fitting', logg_to_plot, 'logg')
        log.write_local_log('fitting', period_to_plot, 'period')
        log.write_local_log('fitting', mid_time_to_plot, 'mid_time')
        log.write_local_log('fitting', rp_over_rs_to_plot, 'rp_over_rs')
        log.write_local_log('fitting', sma_over_rs_to_plot, 'sma_over_rs')
        log.write_local_log('fitting', inclination_to_plot, 'inclination')
        log.write_local_log('fitting', eccentricity_to_plot, 'eccentricity')
        log.write_local_log('fitting', periastron_to_plot, 'periastron')

        ftr_fitting()

        self.running.set(False)
        self.update_window(None)

    def update_window_no_refit(self, *entry):
        self.refit.set(False)
        self.update_window()
        self.refit.set(True)

    def update_window_no_refit_if_auto(self, *entry):
        print(self.manual_planet.get())
        if self.manual_planet.get():
            self.refit.set(True)
        else:
            self.refit.set(False)
        self.update_window()
        self.refit.set(True)


class HOPS:

    def __init__(self):
        self.run_from_reduction = True
        self.run_from_photometry = False
        self.run_from_fitting = False
        self.exit = False

        self.location = os.path.abspath(os.path.dirname(__file__))

        # find current version

        self.current_version = '0.0.0'
        for i in open(os.path.join(self.location, '__init__.py')):
            if len(i.split('__version__')) > 1:
                self.current_version = i.split()[-1][1:-1]

        # find latest version
        self.latest_version = '0.0.0'
        self.latest_version_message = ''

        try:
            for i in urlopen(
                    'https://raw.githubusercontent.com/ExoWorldsSpies/hops/master/hops/__init__.py').readlines():
                if len(str(i).split('__version__')) > 1:
                    self.latest_version = str(i).split()[-1][1:-4]
                if len(str(i).split('__message__')) > 1:
                    message = str(i).split('__message__ = ')[-1][1:-4]
                    self.latest_version_message = message.replace('\\\\n', '\n')
        except:
            pass

        c1 = int(self.current_version.split('.')[0]) * 100 * 100 * 100
        c2 = int(self.current_version.split('.')[1]) * 100 * 100
        c3 = int(self.current_version.split('.')[2]) * 100
        v1 = int(self.latest_version.split('.')[0]) * 100 * 100 * 100
        v2 = int(self.latest_version.split('.')[1]) * 100 * 100
        v3 = int(self.latest_version.split('.')[2]) * 100

        if v1 + v2 + v3 > c1 + c2 + c3:
            self.new_avail = '\nv{0} now available'.format(self.latest_version)
        else:
            self.new_avail = ''

    def check_for_update(self):

        if self.new_avail != '':
            showinfo('Update available',
                     'There is a newer version ({0}) of the code available!\n\n{1}\n\nDownload and install it from:'
                     '\nhttps://www.exoworldsspies.com/en/software'.format(self.latest_version,
                                                                           self.latest_version_message))

    def run(self):

        while not self.exit:

            if self.run_from_reduction:
                ReductionWindow(self)

            if self.run_from_photometry:
                PhotometryWindow(self)

            if self.run_from_fitting:
                FittingWindow(self)


def run_app():
    print('Loading... Please wait for the main window to appear.')

    hops = HOPS()
    hops.run()
