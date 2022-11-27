
__all__=['HOPS']

import os
import shutil
import yaml
import matplotlib.image as mpimg

from urllib.request import urlopen
from tkinter.messagebox import showinfo

from hops import __version__


class HOPSLog:

    def __init__(self):

        # check version

        self.__location__ = os.path.split(os.path.abspath(os.path.dirname(__file__)))[0]
        self.__home__ = os.path.expanduser('~')
        self.version = __version__

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

        c1 = int(self.version.split('.')[0]) * 100 * 100 * 100
        c2 = int(self.version.split('.')[1]) * 100 * 100
        c3 = int(self.version.split('.')[2]) * 100
        v1 = int(self.latest_version.split('.')[0]) * 100 * 100 * 100
        v2 = int(self.latest_version.split('.')[1]) * 100 * 100
        v3 = int(self.latest_version.split('.')[2]) * 100

        if v1 + v2 + v3 > c1 + c2 + c3:
            self.new_avail = 'v{0} now available'.format(self.latest_version)
        else:
            self.new_avail = ''

        # file paths

        self.files = {
            'main_log': self.code_path('log.yaml'),
            'main_log_profile': self.code_path('log_profile.yaml'),
            'main_log_user': self.code_path('log_user.yaml'),
            'local_log': 'log.yaml',
            'local_log_profile': self.home_path('hops_log_profile.yaml'),
            'local_log_user': self.home_path('hops_log_user.yaml'),
            'holomon_logo': self.code_path('holomon.gif'),
            'close_icon': self.code_path('close.gif'),
            'holomon_logo_jpg': self.code_path('logo.jpg'),
            'reduction_output_description': self.code_path('reduction_output_description.txt'),
            'photometry_output_description': self.code_path('photometry_output_description.txt'),
            'fitting_output_description': self.code_path('fitting_output_description.txt')
        }

        self.logo_jpg = mpimg.imread(self.files['holomon_logo'])

        # parameters

        self.params = {}
        self.load_main_log()
        self.load_main_log_user()
        self.load_main_log_profile()
        self.load_local_log_user()
        self.load_local_log_profile()

        try:
            os.chdir(self.get_param('directory'))
            self.load_local_log()
        except:
            self.set_param('directory', 'Choose Directory')
            self.set_param('directory_short', 'Choose Directory')

        self.main_log_profile = self.open_yaml(self.files['main_log_profile'])
        self.local_log_profile = self.open_yaml(self.files['local_log_profile'])

        # fixed parameters

        self.reduction_directory = 'REDUCED_DATA'
        self.reduction_prefix = 'out_'
        self.reduction_trash_directory = 'REDUCED_DATA_TRASH'
        self.all_frames = 'all_frames.pickle'
        self.all_stars = 'all_stars.pickle'
        self.photometry_directory_base = 'PHOTOMETRY'
        self.fov_figure = 'FOV.pdf'
        self.lc_figure = 'RESULTS.pdf'
        self.photometry_file = 'PHOTOMETRY.txt'
        self.light_curve_aperture_file = 'PHOTOMETRY_APERTURE.txt'
        self.light_curve_gauss_file = 'PHOTOMETRY_GAUSS.txt'
        self.fitting_directory_base = 'FITTING'

        # windows

        self.software_name = 'HOlomon Photometric Software'
        if self.new_avail != '':
            self.updates = "UPDATES & USER MANUAL - " + self.new_avail
        else:
            self.updates = "UPDATES & USER MANUAL"

        self.button_height = 2
        self.button_width = 30
        self.entries_width = 20
        self.entries_bd = 3

        self.title_font = ['times', 17, 'bold']
        self.main_font = ['times', 15]
        self.names_font = ['Courier', 15]
        self.button_font = ['times', 15, 'bold']

        self.frame_low_std = -3
        self.frame_upper_std = 20

        # process - keywords

        self.align_star_area_key = 'ALIGN_SA'
        self.align_u0_key = 'ALIGN_U0'
        self.align_x0_key = 'ALIGN_X0'
        self.align_y0_key = 'ALIGN_Y0'
        self.mean_key = 'MEAN'
        self.psf_key = 'PSF'
        self.std_key = 'STD'
        self.skip_key = 'SKIP'
        self.time_key = 'HOPSJD'
        self.airmass_key= 'AIRMASS'

        # process - values

        self.bin_to = 10000

    def get_param(self, parameter):
        return self.params[parameter]

    def set_param(self, parameter, value):
        self.params[parameter] = value

    # log files operations

    def load_main_log(self):

        main_log = self.open_yaml(self.files['main_log'])

        for i in main_log:
            self.params[i] = main_log[i]

    def load_main_log_user(self):

        main_log_user = self.open_yaml(self.files['main_log_user'])

        for i in main_log_user:
            self.params[i] = main_log_user[i]

    def load_main_log_profile(self):

        main_log_profile = self.open_yaml(self.files['main_log_profile'])

        for i in main_log_profile:
            self.params[i] = main_log_profile[i]

    def load_local_log_user(self):

        try:
            local_log_user = self.open_yaml(self.files['local_log_user'])
            main_log_user = self.open_yaml(self.files['main_log_user'])

            for i in local_log_user:
                if isinstance(local_log_user[i], dict):
                    for j in local_log_user[i]:
                        if j in main_log_user:
                            self.set_param(j, local_log_user[i][j])

                elif i in main_log_user:
                    self.set_param(i, local_log_user[i])
        except:
            pass

        self.save_local_log_user()

    def load_local_log_profile(self):

        try:
            local_log_profile = self.open_yaml(self.files['local_log_profile'])
            main_log_profile = self.open_yaml(self.files['main_log_profile'])

            for i in local_log_profile:
                if isinstance(local_log_profile[i], dict):
                    for j in local_log_profile[i]:
                        if j in main_log_profile:
                            self.set_param(j, local_log_profile[i][j])

                elif i in main_log_profile:
                    self.set_param(i, local_log_profile[i])
        except:
            pass

        self.save_local_log_profile()

    def load_local_log(self):

        local_log_path = os.path.join(self.get_param('directory'), self.files['local_log'])

        local_log = self.open_yaml(local_log_path)
        main_log = self.open_yaml(self.files['main_log'])

        for i in local_log:
            if isinstance(local_log[i], dict):
                for j in local_log[i]:
                    if j in main_log:
                        self.set_param(j, local_log[i][j])

            elif i in main_log:
                self.set_param(i, local_log[i])

    def initiate_local_log(self):

        local_log_path = os.path.join(self.get_param('directory'), self.files['local_log'])

        try:
            test = self.open_yaml(local_log_path)
            self.load_local_log()
        except:
            self.save_local_log()

    def save_local_log(self):

        local_log_path = os.path.join(self.get_param('directory'), self.files['local_log'])

        local_log = {}
        main_log = self.open_yaml(self.files['main_log'])

        for i in self.params:
            if i in main_log:
                local_log[i] = self.get_param(i)

        self.save_yaml(local_log, local_log_path)

    def export_local_log(self, path):

        local_log_path = os.path.join(self.get_param('directory'), self.files['local_log'])

        self.save_local_log()
        shutil.copy(local_log_path, path)

    def save_local_log_profile(self, local_log_profile=None):

        if local_log_profile:
            pass
        else:
            local_log_profile = {}
            main_log_profile = self.open_yaml(self.files['main_log_profile'])

            for i in self.params:
                if i in main_log_profile:
                    local_log_profile[i] = self.get_param(i)

        self.save_yaml(local_log_profile, self.files['local_log_profile'])
        self.local_log_profile = self.open_yaml(self.files['local_log_profile'])

    def save_local_log_user(self):

        local_log_user = {}
        main_log_user = self.open_yaml(self.files['main_log_user'])

        for i in self.params:
            if i in main_log_user:
                local_log_user[i] = self.get_param(i)

        self.save_yaml(local_log_user, self.files['local_log_user'])

    def read_local_log_user(self, keyword):
        x = self.open_yaml(self.files['local_log_user'])
        return x[keyword]

    def write_local_log_user(self, keyword, value):
        x = self.open_yaml(self.files['local_log_user'])
        x[keyword] = value
        self.save_yaml(x, self.files['local_log_user'])

    def read_log_profile(self, keyword):
        x = self.open_yaml(self.files['log_profile'])
        if not x[keyword]:
            return ''
        else:
            test = str(x[keyword])
            if test[0] == ' ':
                test = test[1:]
            return test

    def read_local_log_profile(self, keyword):
        try:
            x = self.open_yaml(self.files['local_log_profile'])
            test = str(x[keyword])
            if len(str(test)) > 0:
                if test[0] == ' ':
                    test = test[1:]
            return test
        except KeyError:
            return self.read_log_profile(keyword)

    # operations

    def code_path(self, path):

        return os.path.join(self.__location__, 'hops', 'files', path)

    def home_path(self, path):

        return os.path.join(self.__home__, path)

    def open_yaml(self, path):
        return yaml.load(open(path, 'r'), Loader=yaml.SafeLoader)

    def save_yaml(self, data, path):
        yaml.dump(data, open(path, 'w'), default_flow_style=False)

    def check_for_updates(self):

        if self.new_avail != '':

            showinfo('Update available',
                     'There is a newer version ({0}) of the code available!\n\n{1}\n\nDownload and install it from:'
                     '\nhttps://www.exoworldsspies.com/en/software'.format(self.latest_version,
                                                                           self.latest_version_message))

    def is_version_new(self, version1, version2):

        if not version1:
            return False

        c1 = int(version1.split('.')[0]) * 100 * 100 * 100
        c2 = int(version1.split('.')[1]) * 100 * 100
        c3 = int(version1.split('.')[2]) * 100
        v1 = int(version2.split('.')[0]) * 100 * 100 * 100
        v2 = int(version2.split('.')[1]) * 100 * 100
        v3 = int(version2.split('.')[2]) * 100

        if v1 + v2 + v3 >= c1 + c2 + c3:
            return False
        else:
            return True


