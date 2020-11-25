
import datetime
import os
import shutil
import yaml
import matplotlib.image as mpimg


class HOPSLog:

    def __init__(self):

        self.__location__ = os.path.split(os.path.abspath(os.path.dirname(__file__)))[0]
        self.__home__ = os.path.expanduser('~')

        self.files = {
            'log': self.code_path('log.yaml'),
            'log_profile': self.code_path('log_profile.yaml'),
            'log_user': self.code_path('log_user.yaml'),
            'local_log': 'log.yaml',
            'local_log_profile': self.home_path('hops_log_profile.yaml'),
            'local_log_user': self.home_path('hops_log_user.yaml'),
        }

        try:
            test = self.open_yaml(self.files['local_log_user'])
        except:
            shutil.copy(self.files['log_user'], self.files['local_log_user'])

        try:
            test = self.open_yaml(self.files['local_log_profile'])
        except:
            shutil.copy(self.files['log_profile'], self.files['local_log_profile'])

        self.holomon_logo = self.code_path('holomon.gif')
        self.close_icon = self.code_path('close.gif')
        self.holomon_logo_jpg = mpimg.imread(self.code_path('logo.jpg'))

        self.reduction_output_description = self.code_path('reduction_output_description.txt')
        self.photometry_output_description = self.code_path('photometry_output_description.txt')
        self.fitting_output_description = self.code_path('fitting_output_description.txt')

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

    def write_local_log_profile(self, keyword, value):
        x = self.open_yaml(self.files['local_log_profile'])
        x[keyword] = value

        self.save_yaml(x, self.files['local_log_profile'])

    # main log

    def read_log(self, keyword, keyword2=None):
        x = self.open_yaml(self.files['log'])
        try:
            if keyword2:
                return x[keyword][keyword2]
            else:
                return x[keyword]
        except KeyError:
            return ''

    def read_local_log(self, keyword, keyword2=None):
        try:
            x = self.open_yaml(self.files['local_log'])
            if keyword2:
                return x[keyword][keyword2]
            else:
                return x[keyword]
        except KeyError:
            value = self.read_log(keyword, keyword2)
            self.write_local_log(keyword, value, keyword2)
            return value

    def write_local_log(self, keyword, value, keyword2=None):
        x = self.open_yaml(self.files['local_log'])
        if keyword2:
            x[keyword][keyword2] = value
        else:
            x[keyword] = value

        x['last_changed'] = datetime.datetime.now().isoformat()

        self.save_yaml(x, self.files['local_log'])

    def initiate_local_log_file(self):
        try:
            test = self.open_yaml(self.files['local_log'])
        except:
            shutil.copy(self.files['log'], self.files['local_log'])
            self.write_local_log('pipeline', self.read_local_log_profile('observation_files'), 'observation_files')
            self.write_local_log('reduction', self.read_local_log_profile('bias_files'), 'bias_files')
            self.write_local_log('reduction', self.read_local_log_profile('dark_files'), 'dark_files')
            self.write_local_log('reduction', self.read_local_log_profile('flat_files'), 'flat_files')
            self.write_local_log('reduction', self.read_local_log_profile('bin_fits'), 'bin_fits')
            self.write_local_log('pipeline_keywords', self.read_local_log_profile('exposure_time_key'), 'exposure_time_key')
            self.write_local_log('pipeline_keywords', self.read_local_log_profile('observation_date_key'), 'observation_date_key')
            self.write_local_log('pipeline_keywords', self.read_local_log_profile('observation_time_key'), 'observation_time_key')


    def copy_local_log_file(self, directory):
        shutil.copy(self.files['local_log'], os.path.join(directory, 'log.yaml'))

    # operations

    def code_path(self, path):

        return os.path.join(self.__location__, 'files', path)

    def home_path(self, path):

        return os.path.join(self.__home__, path)

    def open_yaml(self, path):
        return yaml.load(open(path, 'r'), Loader=yaml.SafeLoader)

    def save_yaml(self, data, path):
        yaml.dump(data, open(path, 'w'), default_flow_style=False)


log = HOPSLog()
