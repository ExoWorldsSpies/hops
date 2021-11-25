
import os

import hops.pylightcurve3 as plc

from hops.hops_tools.fits import find_fits_files, get_fits_data
from .application_windows import MainWindow, AddOnWindow

from astroquery.simbad import Simbad
import astropy.coordinates as coord
import astropy.units as u
from urllib.request import urlopen
import yaml
import warnings


filter_map = {'Clear': 'V', 'Luminance': 'V',
              'U': 'U', 'B': 'B', 'V': 'V', 'R': 'R', 'I': 'I', 'H': 'H', 'J': 'J', 'K': 'K',
              'Astrodon ExoPlanet-BB': 'R',
              'UV': 'U', 'Rc': 'R', 'Ic': 'I', 'Re': 'R', 'Ie': 'I',
              }

__location__ = os.path.abspath(os.path.dirname(__file__))
ecc_stars = yaml.load(open(os.path.join(__location__ , 'stars.yaml'), 'r'), Loader=yaml.SafeLoader)
for star in ecc_stars:
    ecc_stars[star]['planets'] = []

ecc_planets = yaml.load(open(os.path.join(__location__ , 'planets.yaml'), 'r'), Loader=yaml.SafeLoader)
for planet in ecc_planets:
    ecc_stars[planet[:-1]]['planets'].append(planet)

def flat_name(name):

    flat_name_list = [
        [' ', ''],
        ['-', ''],
    ]

    name = name.lower()

    for char in flat_name_list:
        name = name.replace(char[0], char[1])

    return name

ecc_stars_flat_names = {flat_name(ff): ecc_stars[ff]['simbad_id'] for ff in ecc_stars}

ecc_stars = {ecc_stars[ff]['simbad_id']: 'Host of {0}'.format(','.join(ecc_stars[ff]['planets'])) for ff in ecc_stars}

class DataTargetWindow(MainWindow):

    def __init__(self, log):

        MainWindow.__init__(self, log, name='HOPS - Data & Target', position=2)

        # extra windows

        self.content_window = AddOnWindow(self, name='Files list', sizex=3, sizey=3, position=1)
        self.header_window = AddOnWindow(self, name='Header keywords list', sizex=3, sizey=3, position=7)
        self.target_window = AddOnWindow(self, name='Select Target')

        # set variables, create and place widgets

        # main window

        self.directory_test = self.Label(text=self.log.get_param('directory_short'))

        self.observation_files = self.Entry(value=self.log.get_param('observation_files'), instance=str,
                                            command=self.update_observation_files)
        self.observation_files_test = self.Label(text=' ')
        self.science_files = 0
        self.science_header = []

        self.show_files_button = self.Button(text='Show files', command=self.content_window.show)

        self.bias_files = self.Entry(value=self.log.get_param('bias_files'), instance=str,
                                     command=self.update_bias_files)
        self.bias_files_test = self.Label(text=' ')

        self.dark_files = self.Entry(value=self.log.get_param('dark_files'), instance=str,
                                     command=self.update_dark_files)
        self.dark_files_test = self.Label(text=' ')

        self.flat_files = self.Entry(value=self.log.get_param('flat_files'), instance=str,
                                     command=self.update_flat_files)
        self.flat_files_test = self.Label(text=' ')

        try:
            initial_bin = self.log.get_param('bin_fits')
            initial_bin = max(1, initial_bin)
            initial_bin = min(initial_bin, 4)
        except:
            print('WARNING the bin_fits parameter should be 1, 2, 3 or 4. Re-setting it to 1.')
            initial_bin = 1

        self.bin_fits = self.DropDown(initial=initial_bin, options=[1, 2, 3, 4], instance=int)

        try:
            initial_crop = self.log.get_param('crop_fits')
            initial_crop = max(0.2, initial_crop)
            initial_crop = min(initial_crop, 1)
        except:
            print('WARNING the crop_fits parameter should be 1, 0.8, 0.6, 0.4, or 0.2. Re-setting it to 1.')
            initial_crop = 1

        self.crop_fits = self.DropDown(initial=initial_crop, options=[1, 0.8, 0.6, 0.4, 0.2], instance=float)

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
                                        options=['exposure start', 'mid-exposure', 'exposure end'])

        self.select_target_button = self.Button(text='CHOOSE TARGET', command=self.update_ra_dec_options)

        self.target_ra_dec = self.Label(text=self.log.get_param('target_ra_dec'))
        self.target_name = self.Label(text=self.log.get_param('target_name'))
        self.target_ra_dec_test = self.Label(text='')

        self.observer = self.Entry(value=self.log.get_param('observer'))
        self.observatory = self.Entry(value=self.log.get_param('observatory'))
        self.telescope = self.Entry(value=self.log.get_param('telescope'))
        self.camera = self.Entry(value=self.log.get_param('camera'))

        filters = ['default'] + list(filter_map.keys())
        if self.log.get_param('filter') not in filters:
            self.log.set_param('filter', 'default')

        self.filter = self.DropDown(initial=self.log.get_param('filter'), options=filters, command=self.check_filter)
        self.filter_test = self.Label(text=' ')

        self.save_and_return_button = self.Button(text='SAVE OPTIONS & RETURN TO MAIN MENU', command=self.save_and_return)
        self.save_and_proceed_button = self.Button(text='SAVE OPTIONS & PROCEED', command=self.save_and_proceed,
                                                   bg='green', highlightbackground='green')

        self.setup_window([
            [],

            [[self.Button(text='CHOOSE DIRECTORY', command=self.update_directory), 1], [self.directory_test, 2, 2]],

            [[self.show_files_button, 2]],

            [[self.Label(text='Name identifier for observation files'), 1],
             [self.observation_files, 2], [self.observation_files_test, 3]],

            [[self.Label(text='Name identifier for bias files'), 1],
             [self.bias_files, 2], [self.bias_files_test, 3]],

            [[self.Label(text='Name identifier for dark files'), 1],
             [self.dark_files, 2], [self.dark_files_test, 3]
             ],

            [[self.Label(text='Name identifier for flat files'), 1],
             [self.flat_files, 2], [self.flat_files_test, 3],
             ],

            [[self.Label(text='Bin fits files (reduced only)'), 1], [self.bin_fits, 2]],
            [[self.Label(text='Subframe (reduced only)'), 1], [self.crop_fits, 2]],

            [[self.show_header_button, 2]],

            [[self.Label(text='Exposure time header keyword'), 1],
             [self.exposure_time_key, 2], [self.exposure_time_key_test, 3]
             ],

            [[self.Label(text='Observation date header keyword\n(no JD, HJD, BJD)'), 1],
             [self.observation_date_key, 2], [self.observation_date_key_test, 3],
             ],

            [[self.Label(text='Observation time header keyword'), 1],
             [self.observation_time_key, 2], [self.observation_time_key_test, 3]],

            [[self.Label(text='Time-stamp\n(which time is saved in your fits files?)'), 1],
             [self.time_stamp, 2]],

            [],

            [[self.select_target_button, 1], [self.target_ra_dec, 2], [self.target_ra_dec_test, 3]],
            [[self.target_name, 2]],

            [[self.Label(text='Observer'), 1], [self.observer, 2], [self.Label(text='OK'), 3]],
            [[self.Label(text='Observatory'), 1], [self.observatory, 2], [self.Label(text='OK'), 3]],
            [[self.Label(text='Telescope'), 1], [self.telescope, 2], [self.Label(text='OK'), 3]],
            [[self.Label(text='Camera'), 1], [self.camera, 2], [self.Label(text='OK'), 3]],
            [[self.Label(text='Filter'), 1], [self.filter, 2], [self.filter_test, 3]],

            [],
            [[self.Button(text='RETURN TO MAIN MENU', command=self.close), 1, 6]],
            [[self.save_and_return_button, 1, 6]],
            [[self.save_and_proceed_button, 1, 6]],
            []

        ], entries_wd=self.log.entries_width)

        # files window

        self.content_list = self.content_window.ListDisplay()

        self.content_window.setup_window([
            [[self.content_list, 0]]
        ])

        content_list = ['  List of files in your directory:', '  ']

        xx = find_fits_files('*')

        for ii in xx:
            content_list.append('  {0}'.format(str(ii).split(os.sep)[-1]))

        self.content_list.update_list(content_list)

        # headers window

        self.header_list = self.header_window.ListDisplay()

        self.header_window.setup_window([
            [[self.header_list, 0]]
        ])

        self.stars = plc.open_dict(os.path.join(plc.databases.oec, 'stars_catalogue.pickle'))
        self.stars = {self.stars[ff]['simbad_id']: self.stars[ff]['planets'][0] for ff in self.stars}

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
            [[self.target_ra_dec_choice_0, 1], [self.auto_target_ra_dec, 2, 2], [self.auto_target_ra_dec_check, 4]],
            [[self.target_ra_dec_choice_1, 1], [self.simbad_target_name, 2, 2], [self.simbad_target_name_check, 4]],
            [[self.target_ra_dec_choice_2, 1], [self.manual_target_ra_dec, 2, 2]],
            [],
            [[self.target_window.Label(text='Target RA/DEC: '), 0], [self.target_ra_dec_2, 1, 3]],
            [[self.target_window.Label(text='Target Name: '), 0], [self.target_name_2, 1, 3]],
            [],
            [[self.target_window.Button(text='  Cancel  ', command=self.target_window.hide), 2],
             [self.target_window.Button(text='  Choose  ', command=self.choose_target), 3]],
            []

        ], entries_wd=self.log.entries_width)

        self.update_observation_files()
        self.update_bias_files()
        self.update_dark_files()
        self.update_flat_files()

    # define functions

    def update_directory(self, *event):

        new_directory = self.askdirectory()

        if len(new_directory) > 0:

            self.log.set_param('directory', new_directory)
            self.log.set_param('directory_short', os.path.split(new_directory)[1])

            self.log.load_main_log()
            self.log.load_local_log_profile()
            try:
                os.chdir(self.log.get_param('directory'))
                self.log.load_local_log()
            except:
                try:
                    os.chdir(self.log.get_param('directory'))
                    self.log.initiate_local_log()
                except:
                    os.chdir(self.log.__home__)
                    self.log.set_param('directory', 'Choose Directory')
                    self.log.set_param('directory_short', 'Choose Directory')

            self.log.save_local_log_user()

            content_list = ['  List of files in your directory:', '  ']

            xx = find_fits_files('*')

            for ii in xx:
                content_list.append('  {0}'.format(str(ii).split(os.sep)[-1]))

            self.content_list.update_list(content_list)

            self.directory_test.set(self.log.get_param('directory_short'))
            self.observation_files.set(self.log.get_param('observation_files'))
            self.bias_files.set(self.log.get_param('bias_files'))
            self.dark_files.set(self.log.get_param('dark_files'))
            self.flat_files.set(self.log.get_param('flat_files'))
            self.bin_fits.set(self.log.get_param('bin_fits'))
            self.crop_fits.set(self.log.get_param('crop_fits'))
            self.target_ra_dec_choice.set(self.log.get_param('target_ra_dec_choice'))
            self.target_ra_dec.set(self.log.get_param('target_ra_dec'))
            self.target_name.set(self.log.get_param('target_name'))
            self.auto_target_ra_dec.set(self.log.get_param('auto_target_ra_dec'))
            self.manual_target_ra_dec.set(self.log.get_param('manual_target_ra_dec'))
            self.simbad_target_name.set(self.log.get_param('simbad_target_name'))
            self.time_stamp.set(self.log.get_param('time_stamp'))
            self.exposure_time_key.set(self.log.get_param('exposure_time_key'))
            self.observation_date_key.set(self.log.get_param('observation_date_key'))
            self.observation_time_key.set(self.log.get_param('observation_time_key'))

            self.update_observation_files()
            self.update_bias_files()
            self.update_dark_files()
            self.update_flat_files()

        self.show()

    def update_observation_files(self, *event):

        check = find_fits_files(self.observation_files.get())
        self.science_files = len(check)

        if self.science_files == 0:

            self.observation_files_test.set('{0} files found\nyou cannot proceed'.format(len(check)))
            self.science_header = []
            header_list = ['  Keywords:      Values:', '  ']

        else:
            self.observation_files_test.set('{0} files found - OK'.format(len(check)))
            self.science_header = []
            header_list = ['  Keywords:      Values:', '  ']

            self.science_header = get_fits_data(check[0])[0].header

            for ii in self.science_header:

                if ii != '':
                    header_list.append('  {0}{1}{2}'.format(str(ii[:10]), ' ' * (15 - len(str(ii[:10]))),
                                                            str(self.science_header[ii])))

        self.header_list.update_list(header_list)

        self.update_exposure_time_key()
        self.update_observation_date_key()
        self.update_observation_time_key()
        self.update_ra_dec()
        self.choose_target()
        self.update_observing_info()
        self.update_save_button()

    #     and then update the header keywards

    def update_bias_files(self, *event):

        check = find_fits_files(self.bias_files.get())
        self.bias_files_test.set('{0} files found - OK'.format(len(check)))

        self.update_save_button()

    def update_dark_files(self, *event):

        check = find_fits_files(self.dark_files.get())
        self.dark_files_test.set('{0} files found - OK'.format(len(check)))

        self.update_save_button()

    def update_flat_files(self, *event):

        check = find_fits_files(self.flat_files.get())
        self.flat_files_test.set('{0} files found - OK'.format(len(check)))

        self.update_save_button()

    def update_ra_dec_options(self, *event):

        self.auto_target_ra_dec.set('----------')
        self.target_ra_dec_choice_0['state'] = self.DISABLED
        auto_found = False

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
                    target = plc.Target(plc.Hours(ra.replace(',', '.')),
                                        plc.Degrees(dec.replace(',', '.')))
                    self.auto_target_ra_dec.set(target.coord)
                    self.target_ra_dec_choice_0['state'] = self.NORMAL
                    auto_found = True
                elif isinstance(ra, float):
                    target = plc.Target(plc.Degrees(ra), plc.Degrees(dec))
                    self.auto_target_ra_dec.set(target.coord)
                    self.target_ra_dec_choice_0['state'] = self.NORMAL
                    auto_found = True
            except:
                pass

        if not auto_found:
            self.auto_target_ra_dec_check.set('Not found')
        else:
            self.auto_target_ra_dec_check.set('')

        if not auto_found and self.target_ra_dec_choice.get() == 0:
            self.target_ra_dec_choice.set(1)

        try:
            urlopen('http://www.google.com', timeout=2)
            self.target_ra_dec_choice_1['state'] = self.NORMAL
            self.simbad_target_name.activate()
            auto_found = True
            self.simbad_target_name_check.set('')
        except:
            self.target_ra_dec_choice_1['state'] = self.DISABLED
            self.simbad_target_name.disable()
            auto_found = False
            self.simbad_target_name_check.set('No connection')

        if not auto_found and self.target_ra_dec_choice.get() == 1:
            self.target_ra_dec_choice.set(2)

        self.update_ra_dec()
        self.target_window.show()

    def update_ra_dec(self, *event):

        if self.target_ra_dec_choice.get() == 0:

            self.simbad_target_name.disable()
            self.simbad_target_name.set('')
            self.manual_target_ra_dec.disable()
            self.manual_target_ra_dec.set('')

            try:
                result_table = Simbad.query_region(coord.SkyCoord(self.auto_target_ra_dec.get(), unit=(u.hourangle, u.deg),
                                                                  frame='icrs'), radius='0d5m0s')[0]
                self.target_ra_dec_2.set('{0} {1}'.format(str(result_table['RA']).replace(' ', ':'),
                                                        str(result_table['DEC']).replace(' ', ':')))

                try:
                    xx = result_table['MAIN_ID'].decode("utf-8")
                except:
                    xx = result_table['MAIN_ID']
                self.target_name_2.set(xx)
            except:
                self.target_ra_dec_2.set(self.auto_target_ra_dec.get())
                self.target_name_2.set(' ')

        elif self.target_ra_dec_choice.get() == 1:

            self.simbad_target_name.activate()
            self.simbad_target_name.widget.focus()
            self.manual_target_ra_dec.disable()
            self.manual_target_ra_dec.set('')

            try:

                flat_input_name = flat_name(self.simbad_target_name.get())

                if flat_input_name in ecc_stars_flat_names:
                    name = ecc_stars_flat_names[flat_input_name]
                elif flat_input_name[:-1] in ecc_stars_flat_names:
                    name = ecc_stars_flat_names[flat_input_name[:-1]]

                else:
                    name = self.simbad_target_name.get().lower()
                    name = name.replace('hd-', 'hd ')
                    name = name.replace('2mass-', '2mass ')
                    name = name.replace('gj-', 'gj ')
                    name = name.replace('hip-', 'hip ')
                    name = name.replace('lhs-', 'lhs ')
                    name = name.replace('tyc-', 'tyc ')
                    name = name.replace('kic-', 'kic ')
                    name = name.replace('tic-', 'tic ')
                    name = name.replace('toi-', 'toi ')
                    name = name.replace('qatar-', 'qatar ')

                    if 'bd' in name:
                        name = name[:3] + name[3:].replace('-', ' ')

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    result_table = Simbad.query_object(name)[0]

                self.target_ra_dec_2.set('{0} {1}'.format(str(result_table['RA']).replace(' ', ':'), str(result_table['DEC']).replace(' ', ':')))
                main_id = result_table['MAIN_ID']
                while '  ' in main_id:
                    main_id = main_id.replace('  ', ' ')
                print(main_id)
                self.target_name_2.set(main_id)
            except:
                self.target_ra_dec_2.set('Coordinates not found')
                self.target_name_2.set(' ')

        else:

            self.simbad_target_name.disable()
            self.simbad_target_name.set('')
            self.manual_target_ra_dec.activate()
            self.manual_target_ra_dec.widget.focus()

            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    result_table = Simbad.query_region(coord.SkyCoord(self.manual_target_ra_dec.get(), unit=(u.hourangle, u.deg),
                                                                  frame='icrs'), radius='0d5m0s')[0]
                self.target_ra_dec_2.set('{0} {1}'.format(str(result_table['RA']).replace(' ', ':'),
                                                        str(result_table['DEC']).replace(' ', ':')))
                main_id = result_table['MAIN_ID']
                while '  ' in main_id:
                    main_id = main_id.replace('  ', ' ')
                self.target_name_2.set(main_id)
            except:
                self.target_ra_dec_2.set(self.manual_target_ra_dec.get())
                self.target_name_2.set(' ')

        print(self.target_name_2.get())
        print(self.target_name_2.get() in ecc_stars)
        if self.target_name_2.get() in ecc_stars:
            self.target_name_2.set('{0} - {1}'.format(self.target_name_2.get(), ecc_stars[self.target_name_2.get()]))

    def choose_target(self, *event):

        self.target_ra_dec.set(self.target_ra_dec_2.get())
        self.target_name.set(self.target_name_2.get())

        try:

            ra_dec_string = self.target_ra_dec.get().split(' ')[0].split(':') + self.target_ra_dec.get().split(' ')[1].split(':')
            target = plc.Target(plc.Hours(ra_dec_string[0], ra_dec_string[1], ra_dec_string[2]),
                                plc.Degrees(ra_dec_string[3], ra_dec_string[4], ra_dec_string[5]))
            self.target_ra_dec_test.set('Coordinates accepted - OK')

        except:
            self.target_ra_dec_test.set('Wrong coordinates\nyou cannot proceed')

        self.update_save_button()
        self.target_window.hide()

    def update_exposure_time_key(self, *event):

        if self.exposure_time_key.get() in self.science_header:
            self.exposure_time_key_test.set('Keyword found - OK')
        else:
            self.exposure_time_key_test.set('Keyword not found\nyou cannot proceed')

        self.update_save_button()

    def update_observation_date_key(self, *event):
        if self.observation_date_key.get() in self.science_header:
            self.observation_date_key_test.set('Keyword found - OK')

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
            self.observation_time_key_test.set('Keyword found - OK')
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

        if self.filter.get() == 'default':
            self.filter_test.set('Filter not valid\nyou cannot proceed')
        else:
            self.filter_test.set('OK')

        self.update_save_button()

    def update_save_button(self, *event):

        if (self.science_files > 0 and 'OK' in self.target_ra_dec_test.get() and
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

    def save(self):

        self.log.set_param('observation_files', self.observation_files.get())
        self.log.set_param('bias_files', self.bias_files.get())
        self.log.set_param('dark_files', self.dark_files.get())
        self.log.set_param('flat_files', self.flat_files.get())
        self.log.set_param('bin_fits', self.bin_fits.get())
        self.log.set_param('crop_fits', self.crop_fits.get())
        self.log.set_param('target_ra_dec_choice', self.target_ra_dec_choice.get())
        self.log.set_param('auto_target_ra_dec', self.auto_target_ra_dec.get())
        self.log.set_param('manual_target_ra_dec', self.manual_target_ra_dec.get())
        self.log.set_param('simbad_target_name', self.simbad_target_name.get())
        self.log.set_param('target_ra_dec', self.target_ra_dec.get())
        self.log.set_param('target_name', self.target_name.get())
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
