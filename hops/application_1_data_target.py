
import os
import numpy as np
import exoclock
import hops.pylightcurve41 as plc

from astropy.io import fits as pf

from hops.hops_tools.fits import *
from hops.hops_tools.image_analysis import bin_frame

from .application_windows import MainWindow, AddOnWindow


filter_map = {'Clear': 'clear',
              'Luminance': 'luminance',
              'U': 'JOHNSON_U',
              'B': 'JOHNSON_B',
              'V': 'JOHNSON_V',
              'R': 'COUSINS_R',
              'I': 'COUSINS_I',
              'H': '2mass_h',
              'J': '2mass_j',
              'K': '2mass_ks',
              'Astrodon ExoPlanet-BB': 'exoplanets_bb',
              'u\'': 'sdss_u',
              'g\'': 'sdss_g',
              'r\'': 'sdss_r',
              'z\'': 'sdss_z',
              'i\'': 'sdss_i',
              }

filter_translations = {
    'Clear': ['clear', 'None', 'Clear'],
    'Luminance': ['luminance', 'Luminance', 'lum', 'Lum', 'L'],
    'U': ['U', 'Uj', 'u'],
    'B': ['B', 'Bj', 'b'],
    'V': ['V', 'Vj', 'v'],
    'R': ['R', 'Rc', 'r'],
    'I': ['I', 'Ic', 'i'],
    'H': ['H', 'h'],
    'J': ['J', 'j'],
    'K': ['K', 'k', 'Ks', 'KS'],
    'Astrodon ExoPlanet-BB': ['exoplanets_bb', 'exoplanets'],
    'u\'': ['up', 'u\''],
    'g\'': ['gp', 'g\''],
    'r\'': ['rp', 'r\''],
    'z\'': ['zp', 'z\''],
    'i\'': ['ip', 'i\''],
}

__location__ = os.path.abspath(os.path.dirname(__file__))


class DataTargetWindow(MainWindow):

    def __init__(self, log):

        MainWindow.__init__(self, log, name='HOPS - Data & Target', position=2)

        # extra windows

        self.content_window = AddOnWindow(self, name='Files list', sizex=3, sizey=3, position=1)
        self.header_window = AddOnWindow(self, name='Header keywords list', sizex=3, sizey=3, position=7)
        self.target_window = AddOnWindow(self, name='Select/Change target')
        self.location_window = AddOnWindow(self, name='Select/Change location')
        self.advanced_settings_window = AddOnWindow(self, name='Advanced settings', position=2)
        self.observer_information_window = AddOnWindow(self, name='Observer information', position=2)

        # set variables, create and place widgets

        # main window

        self.directory_test = self.Label(text=self.log.get_param('directory_short'))

        self.observation_files = self.Entry(value=self.log.get_param('observation_files'), instance=str,
                                            command=self.update_observation_files)
        self.observation_files_test = self.Label(text=' ')
        self.science_files = 0
        self.science_name = ''
        self.science_header = []
        self.science_data = None

        self.show_files_button = self.Button(text='Show files', command=self.content_window.show)

        self.bias_files = self.Entry(value=self.log.get_param('bias_files'), instance=str,
                                     command=self.update_bias_files)
        self.bias_files_test = self.Label(text=' ')

        self.dark_files = self.Entry(value=self.log.get_param('dark_files'), instance=str,
                                     command=self.update_dark_files)
        self.dark_files_test = self.Label(text=' ')

        self.dark_flat_files = self.Entry(value=self.log.get_param('darkf_files'), instance=str,
                                     command=self.update_dark_flat_files)
        self.dark_flat_files_test = self.Label(text=' ')

        self.flat_files = self.Entry(value=self.log.get_param('flat_files'), instance=str,
                                     command=self.update_flat_files)
        self.flat_files_test = self.Label(text=' ')

        self.show_header_button = self.Button(text='Show header', command=self.header_window.show)

        self.exposure_time_key = self.Entry(value=self.log.get_param('exposure_time_key'), instance=str,
                                            command=self.update_exposure_time_key)
        self.exposure_time_key_test = self.Label(text=' ')

        self.observation_date_key = self.Entry(value=self.log.get_param('observation_date_key'), instance=str,
                                               command=self.update_observation_date_key)
        self.observation_date_key_test = self.Label(text=' ')

        self.observation_time_key = self.Entry(value=self.log.get_param('observation_time_key'), instance=str,
                                               command=self.update_observation_time_key)
        self.observation_time_key_test = self.Label(text=' ')

        self.time_stamp = self.DropDown(initial=self.log.get_param('time_stamp'),
                                        instance=str,
                                        options=['exposure start', 'mid-exposure', 'exposure end'])

        self.select_target_button = self.Button(text='Select/Change target', command=self.target_window.show)
        self.target_ra_dec = self.Label(text=self.log.get_param('target_ra_dec'))
        self.target_name = self.Label(text=self.log.get_param('target_name'))
        self.target_ra_dec_test = self.Label(text='')

        self.select_location_button = self.Button(text='Select/Change location', command=self.location_window.show)
        self.location = self.Label(text=self.log.get_param('location'))
        self.location_test = self.Label(text='')

        filters = ['No filter chosen'] + list(filter_map.keys())
        self.filter = self.DropDown(initial=self.log.get_param('filter'),
                                        instance=str,
                                        options=filters, command=self.check_filter)
        self.filter_test = self.Label(text=' ')

        self.show_advanced_settings_button = self.Button(text='Advanced settings',
                                                         command=[self.update_preview, self.advanced_settings_window.show])
        self.show_observer_information_button = self.Button(text='Observer\'s information (optional)',
                                                            command=self.observer_information_window.show)

        self.save_and_return_button = self.Button(text='SAVE OPTIONS &\nRETURN TO MAIN MENU',
                                                  command=self.save_and_return)
        self.save_and_proceed_button = self.Button(text='SAVE OPTIONS &\nPROCEED',
                                                   command=self.save_and_proceed,
                                                   bg='green', highlightbackground='green')

        self.setup_window([
            [],
            [
                [self.Label(text='Directory:'), 1],
                [self.directory_test, 2, 2],
                [self.Button(text='Select/Change directory', command=self.change_directory), 5, 2],
            ],
            [],
            [
                [self.Label(text='Name identifiers:'), 1, 3],
                [self.Label(text='Header information:'), 4, 3]
            ],
            [
                [self.Label(text='Observation files'), 1], [self.observation_files, 2], [self.observation_files_test, 3],
                [self.Label(text='Exposure time key'), 4], [self.exposure_time_key, 5], [self.exposure_time_key_test, 6],
            ],
            [
                [self.Label(text='Bias files'), 1], [self.bias_files, 2], [self.bias_files_test, 3],
                [self.Label(text='Observation date key'), 4], [self.observation_date_key, 5], [self.observation_date_key_test, 6]
            ],
            [
                [self.Label(text='Dark files'), 1], [self.dark_files, 2], [self.dark_files_test, 3],
                [self.Label(text='Observation time key'), 4], [self.observation_time_key, 5], [self.observation_time_key_test, 6]
            ],
            [
                [self.Label(text='Dark-flat files'), 1], [self.dark_flat_files, 2], [self.dark_flat_files_test, 3],
                [self.Label(text='Time-stamp'), 4], [self.time_stamp, 5], [self.Label(text='   OK   '), 6],
            ],
            [
                [self.Label(text='Flat files'), 1], [self.flat_files, 2], [self.flat_files_test, 3],
                [self.Label(text='Filter'), 4], [self.filter, 5], [self.filter_test, 6]
            ],
            [
                [self.show_files_button, 1, 3],
                [self.show_header_button, 4, 3]
            ],
            [],
            [
                [self.Label(text='Location:'), 1], [self.location, 2], [self.location_test, 3],
                [self.Label(text='Target:'), 4], [self.target_ra_dec, 5], [self.target_ra_dec_test, 6],
            ],
            [
                [self.select_location_button, 1,3],  [self.select_target_button, 4, 3],
            ],
            [
                [self.target_name, 4, 3]
            ],
            [],
            [
                [self.show_advanced_settings_button, 1, 3], [self.show_observer_information_button, 4, 3]
            ],
            [],
            [
                [self.save_and_return_button, 1, 3],
                [self.save_and_proceed_button, 4, 3]
            ],
            []

        ], entries_wd=self.log.entries_width)

        # files window

        self.content_list = self.content_window.ListDisplay()

        self.content_window.setup_window([
            [[self.content_list, 0]]
        ])

        # headers window

        self.header_list = self.header_window.ListDisplay()

        self.header_window.setup_window([
            [[self.header_list, 0]]
        ])

        # target window

        self.target_ra_dec_choice = self.target_window.IntVar(self.log.get_param('target_ra_dec_choice'))
        self.target_ra_dec_choice_0 = self.target_window.Radiobutton(
            text='Use the RA/DEC found in the file\'s header:',
            variable=self.target_ra_dec_choice, value=0, command=self.update_ra_dec)
        self.target_ra_dec_choice_1 = self.target_window.Radiobutton(
            text='Provide the name of the target:',
            variable=self.target_ra_dec_choice, value=1, command=self.update_ra_dec)
        self.target_ra_dec_choice_2 = self.target_window.Radiobutton(
            text='Provide the RA/DEC of the target\n(hh:mm:ss +/-dd:mm:ss):',
            variable=self.target_ra_dec_choice, value=2, command=self.update_ra_dec)
        self.auto_target_ra_dec = self.target_window.Label(text=self.log.get_param('auto_target_ra_dec'))
        self.auto_target_ra_dec_check = self.target_window.Label(text='')
        self.simbad_target_name = self.target_window.Entry(value=self.log.get_param('simbad_target_name'), command=self.update_ra_dec)
        self.simbad_target_name_check = self.target_window.Label(text='')
        self.manual_target_ra_dec = self.target_window.Entry(value=self.log.get_param('manual_target_ra_dec'), command=self.update_ra_dec)
        self.target_ra_dec_2 = self.target_window.Label(text=self.log.get_param('target_ra_dec'))
        self.target_name_2 = self.target_window.Label(text=self.log.get_param('target_name'))

        self.target_window.setup_window([
            [],
            [[self.target_window.Label(text='How would you like to choose your target? Select one of the three options:'), 0, 5]],
            [],
            [[self.target_ra_dec_choice_0, 1], [self.auto_target_ra_dec, 2, 2]],
            [[self.target_ra_dec_choice_1, 1], [self.simbad_target_name, 2, 2]],
            [[self.target_ra_dec_choice_2, 1], [self.manual_target_ra_dec, 2, 2]],
            [],
            [[self.target_window.Label(text='Target RA/DEC: '), 0], [self.target_ra_dec_2, 1, 3]],
            [[self.target_window.Label(text='Target Name: '), 0], [self.target_name_2, 1, 3]],
            [],
            [[self.target_window.Button(text='  Cancel  ', command=self.target_window.hide), 2],
             [self.target_window.Button(text='  Choose  ', command=self.choose_target), 3]],
            []

        ], entries_wd=self.log.entries_width)

        # location window

        self.location_choice = self.location_window.IntVar(self.log.get_param('location_choice'))
        self.location_choice_0 = self.location_window.Radiobutton(
            text='Use the LAT/LONG found in the file\'s header:',
            variable=self.location_choice, value=0, command=self.update_location)
        self.location_choice_1 = self.location_window.Radiobutton(
            text='Use the LAT/LONG found in your profile:',
            variable=self.location_choice, value=1, command=self.update_location)
        self.location_choice_2 = self.location_window.Radiobutton(
            text='Provide the LAT/LONG of the location:\n(+/-dd:mm:ss +/-dd:mm:ss)\npositive LAT is NORTH\npositive LONG is EAST',
            variable=self.location_choice, value=2, command=self.update_location)
        self.auto_location = self.location_window.Label(text=self.log.get_param('auto_location'))
        self.profile_location = self.location_window.Label(text=self.log.get_param('profile_location'))
        self.manual_location = self.location_window.Entry(value=self.log.get_param('manual_location'), command=self.update_location)

        self.location_window.setup_window([
            [],
            [[self.location_window.Label(text='How would you like to choose your location? Select one of the two options:'), 0, 5]],
            [],
            [[self.location_choice_0, 1], [self.auto_location, 2, 2]],
            [[self.location_choice_1, 1], [self.profile_location, 2, 2]],
            [[self.location_choice_2, 1], [self.manual_location, 2, 2]],
            [],
            [[self.location_window.Button(text='  Cancel  ', command=self.location_window.hide), 2],
             [self.location_window.Button(text='  Choose  ', command=self.choose_location), 3]],
            []

        ], entries_wd=self.log.entries_width)

        # advanced settings window

        self.faint_target_mode = self.advanced_settings_window.CheckButton(
            text='UltraShortExposure mode (use for occultations)',
            initial=self.log.get_param('faint_target_mode'))

        self.moving_target_mode = self.advanced_settings_window.CheckButton(
            text='MovingTarget mode (use for asteroids).',
            initial=self.log.get_param('moving_target_mode')
        )

        self.rotating_field_mode = self.advanced_settings_window.CheckButton(
            text='RotatingField mode (use for AltAz mounts without rotator).',
            initial=self.log.get_param('rotating_field_mode')
        )

        self.colour_camera_mode = self.advanced_settings_window.CheckButton(
            text='ColourCamera mode (use for colour cameras).',
            initial=self.log.get_param('colour_camera_mode')
        )

        self.bin_fits = self.advanced_settings_window.DropDown(initial=self.log.get_param('bin_fits'),
                                                               options=[1, 2, 3, 4],
                                                               instance=int, command=self.update_preview_final)
        self.crop_edge_pixels = self.advanced_settings_window.DropDown(initial=self.log.get_param('crop_edge_pixels'),
                                                                       options=np.int_(np.arange(0, 50, 1)),
                                                                       instance=int, command=self.update_preview_final)

        self.crop_x1 = self.advanced_settings_window.Entry(value=self.log.get_param('crop_x1'), instance=int)
        self.crop_x2 = self.advanced_settings_window.Entry(value=self.log.get_param('crop_x2'), instance=int)
        self.crop_y1 = self.advanced_settings_window.Entry(value=self.log.get_param('crop_y1'), instance=int)
        self.crop_y2 = self.advanced_settings_window.Entry(value=self.log.get_param('crop_y2'), instance=int)

        self.initial_figure = self.advanced_settings_window.FitsWindow(
            show_controls=True, show_axes=True,
            subplots_adjust=(0.07, 0.99, 0.07, 0.99)
        )
        self.final_figure = self.advanced_settings_window.FitsWindow(
            show_controls=True, show_axes=True,
            subplots_adjust=(0.07, 0.99, 0.07, 0.99)
        )

        self.initial_figure_size = self.advanced_settings_window.Label(text='Image before reduction: 0 MB')
        self.final_figure_size = self.advanced_settings_window.Label(text='Image after reduction: 0 MB')

        self.advanced_settings_window.setup_window([
            [
                [self.advanced_settings_window.Label(
                    text='To crop the original image (on the left panel),\n'
                         'zoom in to the area you want to select and then press "Select area".'
                ), 0],
                [self.advanced_settings_window.Button(text='Select area', command=self.crop_to_selected_area), 1],
                [self.advanced_settings_window.Label(text='Binning (use 2 for coloured cameras)'), 3],
                [self.bin_fits, 4],
            ],
            [
                [self.advanced_settings_window.Label(
                    text='Remove pixels from the edge of the image,\n'
                         'useful for cameras with GPS data in the first rows/columns.'
                ), 0],
                [self.crop_edge_pixels, 1],
                [self.faint_target_mode, 3, 2]
            ],
            [
                [self.moving_target_mode, 3, 2]
            ],
            [
                [self.rotating_field_mode, 3, 2]
            ],
            [
                [self.colour_camera_mode, 3, 2]
            ],
            [],
            [[self.initial_figure_size, 0, 3], [self.final_figure_size, 3, 2]],
            [[self.initial_figure, 0, 3], [self.final_figure, 3, 2]],
            [],
        ])

        # observer information window

        self.observer = self.observer_information_window.Entry(value=self.log.get_param('observer'))
        self.observatory = self.observer_information_window.Entry(value=self.log.get_param('observatory'))
        self.telescope = self.observer_information_window.Entry(value=self.log.get_param('telescope'))
        self.camera = self.observer_information_window.Entry(value=self.log.get_param('camera'))

        self.observer_information_window.setup_window([
            [],
            [],
            [
                [self.observer_information_window.Label(text='Observer Name'), 1],
                [self.observer, 2, 2],
                [self.observer_information_window.Label(text='Observatory Name'), 4],
                [self.observatory, 5, 2]
            ],
            [
                [self.observer_information_window.Label(text='Telescope Name'), 1],
                [self.telescope, 2, 2],
                [self.observer_information_window.Label(text='Camera Name'), 4], [self.camera, 5, 2],
            ],
            [],
            []
        ])

    def check_dir(self):

        self.disable()

        dir_ok = False

        if os.path.isdir(self.log.get_param('directory')):
            try:
                self.log.load_local_log()
                dir_ok = True
            except:
                pass

        if not dir_ok:
            new_directory = self.askdirectory()
            if new_directory == '':
                self.close()
            else:
                try:
                    os.chdir(new_directory)
                    self.log.set_param('directory', new_directory)
                    self.log.set_param('directory_short', os.path.split(new_directory)[1])
                    self.log.load_main_log()
                    self.log.load_local_log_profile()
                    try:
                        self.log.load_local_log()
                    except:
                        self.log.initiate_local_log()
                    self.log.save_local_log_user()
                    dir_ok = True
                except:
                    self.showinfo('Not valid directory', 'Not valid directory')
                    os.chdir(self.log.__home__)
                    self.log.set_param('directory', 'Choose Directory')
                    self.log.set_param('directory_short', 'Choose Directory')
                    self.log.save_local_log_user()
                    self.close()

        if dir_ok:
            self.activate()
            self.update_directory()

    # define functions

    def change_directory(self, *event):

        current_directory = self.log.get_param('directory')
        new_directory = self.askdirectory()
        if new_directory != '':
            try:
                os.chdir(new_directory)
                self.log.set_param('directory', new_directory)
                self.log.set_param('directory_short', os.path.split(new_directory)[1])
                self.log.load_main_log()
                self.log.load_local_log_profile()
                try:
                    self.log.load_local_log()
                except:
                    self.log.initiate_local_log()
                self.log.save_local_log_user()
                self.update_directory()
            except:
                self.showinfo('Not valid directory', 'Not valid directory')
                os.chdir(current_directory)
                self.log.set_param('directory', current_directory)
                self.log.set_param('directory_short', os.path.split(current_directory)[1])
                self.log.load_main_log()
                self.log.load_local_log_profile()
                self.log.initiate_local_log()
                self.log.save_local_log_user()

    def update_directory(self, *event):

        content_list = ['  List of files in your directory:', '  ']

        xx = find_fits_files('*')

        for ii in xx:
            content_list.append('  {0}'.format(str(ii).split(os.sep)[-1]))

        self.content_list.update_list(content_list)

        self.directory_test.set(self.log.get_param('directory_short'))
        self.observation_files.set(self.log.get_param('observation_files'))
        self.bias_files.set(self.log.get_param('bias_files'))
        self.dark_files.set(self.log.get_param('dark_files'))
        self.dark_flat_files.set(self.log.get_param('darkf_files'))
        self.flat_files.set(self.log.get_param('flat_files'))
        self.bin_fits.set(self.log.get_param('bin_fits'))
        self.crop_x1.set(self.log.get_param('crop_x1'))
        self.crop_x2.set(self.log.get_param('crop_x2'))
        self.crop_y1.set(self.log.get_param('crop_y1'))
        self.crop_y2.set(self.log.get_param('crop_y2'))
        self.faint_target_mode.set(self.log.get_param('faint_target_mode'))
        self.moving_target_mode.set(self.log.get_param('moving_target_mode'))
        self.rotating_field_mode.set(self.log.get_param('rotating_field_mode'))
        self.colour_camera_mode.set(self.log.get_param('colour_camera_mode'))
        self.target_ra_dec_choice.set(self.log.get_param('target_ra_dec_choice'))
        self.target_ra_dec.set(self.log.get_param('target_ra_dec'))
        self.target_name.set(self.log.get_param('target_name'))
        self.auto_target_ra_dec.set(self.log.get_param('auto_target_ra_dec'))
        self.manual_target_ra_dec.set(self.log.get_param('manual_target_ra_dec'))
        self.location_choice.set(self.log.get_param('location_choice'))
        self.location.set(self.log.get_param('location'))
        self.auto_location.set(self.log.get_param('auto_location'))
        self.profile_location.set(self.log.get_param('profile_location'))
        self.manual_location.set(self.log.get_param('manual_location'))
        self.simbad_target_name.set(self.log.get_param('simbad_target_name'))
        self.time_stamp.set(self.log.get_param('time_stamp'))
        self.exposure_time_key.set(self.log.get_param('exposure_time_key'))
        self.observation_date_key.set(self.log.get_param('observation_date_key'))
        self.observation_time_key.set(self.log.get_param('observation_time_key'))
        self.filter.set(self.log.get_param('filter'))

        self.update_observation_files()
        self.update_bias_files()
        self.update_dark_files()
        self.update_dark_flat_files()
        self.update_flat_files()

    def update_observation_files(self, *event):

        check = find_fits_files(self.observation_files.get())

        self.science_files = len(check)

        if self.science_files == 0:

            self.observation_files_test.set('0 files found\nyou cannot proceed')
            self.science_name = ''
            self.science_header = []
            self.science_data = None

            header_list = ['  Keywords:      Values:', '  ']
            self.header_list.update_list(header_list)

            self.update_exposure_time_key()
            self.update_observation_date_key()
            self.update_observation_time_key()
            self.update_ra_dec_options()
            self.choose_target()
            self.update_location_options()
            self.choose_location()
            self.update_observing_info()

        else:

            self.observation_files_test.set('{0} files - OK'.format(len(check)))

            if check[0] != self.science_name:

                self.science_name = check[0]
                self.science_data, self.science_header = get_fits_data_and_header(check[0])

                header_list = ['  Keywords:      Values:', '  ']
                for ii in self.science_header:
                    if ii != '':
                        header_list.append('  {0}{1}{2}'.format(str(ii[:10]), ' ' * (15 - len(str(ii[:10]))),
                                                                str(self.science_header[ii])))
                self.header_list.update_list(header_list)

                self.update_exposure_time_key()
                self.update_observation_date_key()
                self.update_observation_time_key()
                self.update_ra_dec_options()
                self.choose_target()
                self.update_location_options()
                self.choose_location()
                self.update_observing_info()

    def update_preview(self):

        if self.crop_x2.get() == 0:
           self.crop_x2.set(len(self.science_data[0]))
        if self.crop_y2.get() == 0:
           self.crop_y2.set(len(self.science_data))
        self.initial_figure_size.set(
            'Image before reduction: {0} MB'.format(
                round(os.stat(self.science_name).st_size/1024.0/1024.0, 2)
            )
        )
        self.final_figure_size.set(
            'Image after reduction: {0} MB'.format(
                round(os.stat(self.science_name).st_size/1024.0/1024.0, 2)
            )
        )
        self.initial_figure.load_fits(self.science_data, self.science_header, input_name=self.science_name)
        self.initial_figure.adjust_size()
        self.final_figure.load_fits(self.science_data, self.science_header,
                                    input_name='This is only an example of the reduction output.  '
                                               'Select the area to crop on the left image.')
        self.final_figure.adjust_size()
        self.update_preview_final()

    def update_preview_final(self):

        data_frame = np.ones_like(self.science_data) * self.science_data

        data_frame = data_frame[self.crop_y1.get(): self.crop_y2.get()]
        data_frame = data_frame[:, self.crop_x1.get(): self.crop_x2.get()]

        if self.crop_edge_pixels.get() > 0:
            data_frame = data_frame[self.crop_edge_pixels.get():-self.crop_edge_pixels.get(),
                         self.crop_edge_pixels.get():-self.crop_edge_pixels.get()]

        if self.bin_fits.get() > 1:
            data_frame = bin_frame(data_frame, self.bin_fits.get())
        try:
            os.remove('.test_size.fits')
        except:
            pass

        hdu = pf.CompImageHDU(data=np.array(data_frame, dtype=np.int32))
        plc.save_fits(pf.HDUList([pf.PrimaryHDU(), hdu]), '.test_size.fits')

        self.final_figure_size.set(
            'Image after reduction: {0} MB'.format(round(os.stat('.test_size.fits').st_size/1024.0/1024.0, 2))
        )

        os.remove('.test_size.fits')

        self.final_figure.load_fits(data_frame, self.science_header,
                                    input_name='This is only an example of the reduction output.  '
                                                    'Select the area to crop on the left image.')
        # self.final_figure.draw()

    def update_bias_files(self, *event):

        check = find_fits_files(self.bias_files.get())
        self.bias_files_test.set('{0} files - OK'.format(len(check)))

    def update_dark_files(self, *event):

        check = find_fits_files(self.dark_files.get())
        self.dark_files_test.set('{0} files - OK'.format(len(check)))

    def update_dark_flat_files(self, *event):

        check = find_fits_files(self.dark_flat_files.get())
        self.dark_flat_files_test.set('{0} files - OK'.format(len(check)))

    def update_flat_files(self, *event):

        check = find_fits_files(self.flat_files.get())
        self.flat_files_test.set('{0} files - OK'.format(len(check)))

    def update_ra_dec_options(self, *event):

        self.auto_target_ra_dec.set('Not found')
        self.target_ra_dec_choice_0['state'] = self.DISABLED

        ra = None
        for key in self.log.get_param('target_ra_key').split(','):
            if key in self.science_header:
                ra = self.science_header[key]
                break

        dec = None
        for key in self.log.get_param('target_dec_key').split(','):
            if key in self.science_header:
                dec = self.science_header[key]
                break

        if ra and dec:
            try:
                if isinstance(ra, str):
                    target = exoclock.FixedTarget(exoclock.Hours(ra.replace(',', '.').replace(' ', ':')),
                                                  exoclock.Degrees(dec.replace(',', '.').replace(' ', ':')))
                    self.auto_target_ra_dec.set(target.coord())
                    self.target_ra_dec_choice_0['state'] = self.NORMAL
                elif isinstance(ra, float):
                    target = exoclock.FixedTarget(exoclock.Degrees(ra), exoclock.Degrees(dec))
                    self.auto_target_ra_dec.set(target.coord())
                    self.target_ra_dec_choice_0['state'] = self.NORMAL
            except:
                pass

        self.update_ra_dec()

    def update_ra_dec(self, *event):

        if self.target_ra_dec_choice.get() == 0:

            self.simbad_target_name.disable()
            self.simbad_target_name.set('')
            self.manual_target_ra_dec.disable()
            self.manual_target_ra_dec.set('')

            try:
                target = exoclock.FixedTarget(
                    exoclock.Hours(self.auto_target_ra_dec.get().split(' ')[0]),
                    exoclock.Degrees(self.auto_target_ra_dec.get().split(' ')[1])
                )

                try:
                    nearest = exoclock.simbad_search_by_coordinates(target.ra, target.dec,
                                                                    radius=exoclock.Degrees(0.25))
                    nearest._notes = ''
                    if len(nearest.all_names) > 0:
                        test = list(set(nearest.all_names).intersection(list(exoclock.exoclock_data.ecc()['hosts'])))
                        if len(test) > 0:
                            nearest._notes = '/ Host of: ' + ', '.join(exoclock.exoclock_data.ecc()['hosts'][test[0]])

                    self.target_ra_dec_2.set(nearest.coord())
                    if nearest._notes:
                        self.target_name_2.set(nearest.name + ' ' + nearest._notes)
                    else:
                        self.target_name_2.set(nearest.name)
                except:
                    self.target_ra_dec_2.set(self.auto_target_ra_dec.get())
                    self.target_name_2.set('Name not resolved')
            except:
                self.target_ra_dec_2.set('Coordinates not found')
                self.target_name_2.set('Name not resolved')

        elif self.target_ra_dec_choice.get() == 1:

            self.simbad_target_name.activate()
            self.simbad_target_name.widget.focus()
            self.manual_target_ra_dec.disable()
            self.manual_target_ra_dec.set('')

            try:
                nearest = exoclock.simbad_search_by_name(self.simbad_target_name.get())
                nearest._notes = ''
                if len(nearest.all_names) > 0:
                    test = list(set(nearest.all_names).intersection(list(exoclock.exoclock_data.ecc()['hosts'])))
                    if len(test) > 0:
                        nearest._notes = '/ Host of: ' + ', '.join(exoclock.exoclock_data.ecc()['hosts'][test[0]])

                self.target_ra_dec_2.set(nearest.coord())
                if nearest._notes:
                    self.target_name_2.set(nearest.name + ' ' + nearest._notes)
                else:
                    self.target_name_2.set(nearest.name)
            except:
                self.target_ra_dec_2.set('Coordinates not found')
                self.target_name_2.set('Name not resolved')

        else:

            self.simbad_target_name.disable()
            self.simbad_target_name.set('')
            self.manual_target_ra_dec.activate()
            self.manual_target_ra_dec.widget.focus()

            try:
                if len(self.manual_target_ra_dec.get().split(':')) == 5:
                    target = exoclock.FixedTarget(
                        exoclock.Hours(self.manual_target_ra_dec.get().split(' ')[0]),
                        exoclock.Degrees(self.manual_target_ra_dec.get().split(' ')[1])
                    )
                    nearest = exoclock.simbad_search_by_coordinates(target.ra, target.dec,
                                                                    radius=exoclock.Degrees(0.25))
                    nearest._notes = ''
                    if len(nearest.all_names) > 0:
                        test = list(set(nearest.all_names).intersection(list(exoclock.exoclock_data.ecc()['hosts'])))
                        if len(test) > 0:
                            nearest._notes = '/ Host of: ' + ', '.join(exoclock.exoclock_data.ecc()['hosts'][test[0]])

                    self.target_ra_dec_2.set(nearest.coord())
                    if nearest._notes:
                        self.target_name_2.set(nearest.name + ' ' + nearest._notes)
                    else:
                        self.target_name_2.set(nearest.name)
                else:
                    self.target_ra_dec_2.set(self.manual_target_ra_dec.get())
                    self.target_name_2.set(' ')
            except:
                self.target_ra_dec_2.set(self.manual_target_ra_dec.get())
                self.target_name_2.set(' ')

    def choose_target(self, *event):

        self.target_ra_dec.set(self.target_ra_dec_2.get())
        self.target_name.set(self.target_name_2.get())

        try:

            _ = exoclock.FixedTarget(exoclock.Hours(self.target_ra_dec.get().split(' ')[0]),
                                     exoclock.Degrees(self.target_ra_dec.get().split(' ')[1]))
            self.target_ra_dec_test.set('   OK   ')

        except:
            self.target_ra_dec_test.set('Wrong coordinates\nyou cannot proceed')

        self.update_save_button()
        self.target_window.hide()

    def update_location_options(self, *event):

        self.auto_location.set('Not found')
        self.location_choice_0['state'] = self.DISABLED

        lat = None
        for key in self.log.get_param('observatory_latitude_key').split(','):
            if key in self.science_header:
                lat = self.science_header[key]
                break

        long = None
        for key in self.log.get_param('observatory_longitude_key').split(','):
            if key in self.science_header:
                long = self.science_header[key]
                break

        if lat and long:
            try:
                if isinstance(lat, str):
                    observatory = exoclock.Observatory(exoclock.Degrees(lat.replace(',', '.')),
                                                  exoclock.Degrees(long.replace(',', '.')))
                    self.auto_location.set(observatory.coord())
                    self.location_choice_0['state'] = self.NORMAL
                elif isinstance(lat, float):
                    observatory = exoclock.Observatory(exoclock.Degrees(lat), exoclock.Degrees(long))
                    self.auto_location.set(observatory.coord())
                    self.location_choice_0['state'] = self.NORMAL
            except:
                pass

        self.profile_location.set('Not found')
        self.location_choice_1['state'] = self.DISABLED
        try:
            observatory = exoclock.Observatory(exoclock.Degrees(self.log.get_param('observatory_lat')),
                                          exoclock.Degrees(self.log.get_param('observatory_long')))
            self.profile_location.set(observatory.coord())
            self.location_choice_1['state'] = self.NORMAL
        except:
            pass

        self.update_location()

    def update_location(self, *event):

        if self.location_choice.get() == 0:

            self.manual_location.disable()
            self.location.set(self.auto_location.get())

        elif self.location_choice.get() == 1:

            self.manual_location.disable()
            self.location.set(self.profile_location.get())

        else:

            self.manual_location.activate()
            self.manual_location.widget.focus()
            self.location.set(self.manual_location.get())

    def choose_location(self, *event):

        try:

            observatory = exoclock.Observatory(exoclock.Degrees(self.location.get().split(' ')[0]),
                                          exoclock.Degrees(self.location.get().split(' ')[1]))
            self.location_test.set('   OK   ')

        except:
            self.location_test.set('Wrong coordinates\nyou cannot proceed')

        self.update_save_button()
        self.location_window.hide()

    def update_exposure_time_key(self, *event):

        if self.exposure_time_key.get() in self.science_header:
            self.exposure_time_key_test.set('   OK   ')
        else:
            self.exposure_time_key_test.set('Keyword not found\nyou cannot proceed')

        self.update_save_button()

    def update_observation_date_key(self, *event):
        if self.observation_date_key.get() in self.science_header:
            self.observation_date_key_test.set('   OK   ')

            if len(self.science_header[self.observation_date_key.get()].split('T')) == 2:
                self.observation_time_key.set(self.observation_date_key.get())
                self.observation_time_key.disable()
            else:
                self.observation_time_key.activate()

            self.update_observation_time_key()

        else:
            self.observation_date_key_test.set('Keyword not found\nyou cannot proceed')

        self.update_save_button()

    def update_observation_time_key(self, *event):

        if self.observation_time_key.get() in self.science_header:
            self.observation_time_key_test.set('   OK   ')
        else:
            self.observation_time_key_test.set('Keyword not found\nyou cannot proceed')

        self.update_save_button()

    def update_observing_info(self, *event):

        if self.telescope.get() == 'default':
            for key in self.log.get_param('telescope_key').split(','):
                if key in self.science_header:
                    self.telescope.set(self.science_header[key])
                    break
            if self.telescope.get() == 'default':
                self.telescope.set(self.log.get_param('telescope'))

        if self.camera.get() == 'default':
            for key in self.log.get_param('camera_key').split(','):
                if key in self.science_header:
                    self.camera.set(self.science_header[key])
                    break
            if self.camera.get() == 'default':
                self.camera.set(self.log.get_param('camera'))

        if self.observer.get() == 'default':
            for key in self.log.get_param('observer_key').split(','):
                if key in self.science_header:
                    self.observer.set(self.science_header[key])
                    break
            if self.observer.get() == 'default':
                self.observer.set(self.log.get_param('observer'))

        if self.observatory.get() == 'default':
            for key in self.log.get_param('observatory_key').split(','):
                if key in self.science_header:
                    self.observatory.set(self.science_header[key])
                    break
            if self.observatory.get() == 'default':
                self.observatory.set(self.log.get_param('observatory'))

        self.check_filter()

    def check_filter(self):

        if self.filter.get() not in filter_map:

            for filter in filter_translations:
                if self.filter.get() in filter_translations[filter]:
                    self.filter.set(filter)

        if self.filter.get() not in filter_map:

            for key in self.log.get_param('filter_key').split(','):
                if key in self.science_header:
                    self.filter.set(self.science_header[key])
                    break

            for filter in filter_translations:
                if self.filter.get() in filter_translations[filter]:
                    self.filter.set(filter)

        if self.filter.get() not in filter_map:

            self.filter.set(self.log.get_param('filter'))

            for filter in filter_translations:
                if self.filter.get() in filter_translations[filter]:
                    self.filter.set(filter)

        if self.filter.get() not in filter_map:
            self.filter.set('No filter chosen')
            self.filter_test.set('Filter not valid\nyou cannot proceed')
        else:
            self.filter_test.set('   OK   ')

        self.update_save_button()

    def update_save_button(self, *event):

        if (self.science_files > 0 and
                'OK' in self.location_test.get() and
                'OK' in self.target_ra_dec_test.get() and
                'OK' in self.exposure_time_key_test.get() and
                'OK' in self.observation_date_key_test.get() and
                'OK' in self.observation_time_key_test.get() and
                'OK' in self.filter_test.get()
        ):
            self.save_and_return_button.activate()
            self.save_and_proceed_button.activate()

        else:
            self.save_and_return_button.disable()
            self.save_and_proceed_button.disable()

    def crop_to_selected_area(self):

        x1, x2 = self.initial_figure.ax.get_xlim()
        y1, y2 = self.initial_figure.ax.get_ylim()

        self.crop_x1.set(int(max(0, x1)))
        self.crop_x2.set(int(min(x2, len(self.science_data[0]))))
        self.crop_y1.set(int(max(0, y1)))
        self.crop_y2.set(int(min(y2, len(self.science_data))))

        self.update_preview_final()

    def save(self):

        self.log.set_param('observation_files', self.observation_files.get())
        self.log.set_param('bias_files', self.bias_files.get())
        self.log.set_param('dark_files', self.dark_files.get())
        self.log.set_param('darkf_files', self.dark_flat_files.get())
        self.log.set_param('flat_files', self.flat_files.get())
        self.log.set_param('bin_fits', self.bin_fits.get())
        self.log.set_param('crop_x1', self.crop_x1.get())
        self.log.set_param('crop_x2', self.crop_x2.get())
        self.log.set_param('crop_y1', self.crop_y1.get())
        self.log.set_param('crop_y2', self.crop_y2.get())
        self.log.set_param('faint_target_mode', self.faint_target_mode.get())
        if self.faint_target_mode.get():
            self.log.set_param('centroids_snr', 2)
            self.log.set_param('psf_guess', 1)
            self.log.set_param('stars_snr', 2)
        self.log.set_param('moving_target_mode', self.moving_target_mode.get())
        self.log.set_param('rotating_field_mode', self.rotating_field_mode.get())
        self.log.set_param('colour_camera_mode', self.colour_camera_mode.get())
        self.log.set_param('crop_edge_pixels', self.crop_edge_pixels.get())
        self.log.set_param('target_ra_dec_choice', self.target_ra_dec_choice.get())
        self.log.set_param('auto_target_ra_dec', self.auto_target_ra_dec.get())
        self.log.set_param('manual_target_ra_dec', self.manual_target_ra_dec.get())
        self.log.set_param('simbad_target_name', self.simbad_target_name.get())
        self.log.set_param('target_ra_dec', self.target_ra_dec.get())
        self.log.set_param('target_name', self.target_name.get())
        self.log.set_param('location_choice', self.location_choice.get())
        self.log.set_param('auto_location', self.auto_location.get())
        self.log.set_param('manual_location', self.manual_location.get())
        self.log.set_param('location', self.location.get())
        self.log.set_param('time_stamp', self.time_stamp.get())
        self.log.set_param('exposure_time_key', self.exposure_time_key.get())
        self.log.set_param('observation_date_key', self.observation_date_key.get())
        self.log.set_param('observation_time_key', self.observation_time_key.get())
        self.log.set_param('observer', self.observer.get())
        self.log.set_param('observatory', self.observatory.get())
        self.log.set_param('telescope', self.telescope.get())
        self.log.set_param('camera', self.camera.get())
        self.log.set_param('filter', self.filter.get())

        self.log.set_param('data_target_complete', True)
        self.log.set_param('data_target_version', self.log.version)

        self.log.save_local_log_user()
        self.log.save_local_log()

    def save_and_return(self):

        self.save()
        self.log.set_param('proceed', False)
        self.close()

    def save_and_proceed(self):

        self.save()
        self.log.set_param('proceed', True)
        self.close()
