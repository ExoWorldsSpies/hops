
import os
import glob
import time
import shutil
import datetime
import numpy as np
from scipy.interpolate import interp1d

from pylightcurve.processes.files import open_dict, open_yaml, save_dict, download, open_dict_online
from pylightcurve import __version__

try:
    import zipfile
    download_zip = True
except:
    download_zip = False

databases_file = '__databases__.pickle'
package_name = 'pylightcurve4'
github_link = 'https://github.com/ucl-exoplanets/pylightcurve/raw/master/pylightcurve/__databases__.pickle?raw=true'


class PlcData:

    def __init__(self, _reset=False, _test=False):

        self.package_name = package_name
        self.version = '.'.join(__version__.split('.')[:2])

        self.build_in_databases_file_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), databases_file)

        self.databases_directory_path = os.path.join(os.path.abspath(os.path.expanduser('~')),
                                                     '.{0}'.format(self.package_name))

        self.databases_file_path = os.path.join(self.databases_directory_path, databases_file)
        self.databases_file_path_new = os.path.join(self.databases_directory_path, databases_file + '_new')

        # initiate databases

        if not os.path.isdir(self.databases_directory_path):
            os.mkdir(self.databases_directory_path)

        test = False

        if test:

            shutil.copy(self.build_in_databases_file_path, self.databases_file_path)

        else:

            if not os.path.isfile(self.databases_file_path):
                shutil.copy(self.build_in_databases_file_path, self.databases_file_path)

            # check for updates in the databases (identified on github)

            test_online_db = open_dict_online(github_link)
            test_local_db = open_dict(self.databases_file_path)
            if test_online_db and test_online_db != test_local_db:
                save_dict(test_online_db, self.databases_file_path)

        # load databases

        self.databases = open_dict(self.databases_file_path)

        self.exotethys_loaded = self._setup_database('exotethys')
        self.ephemerides_loaded = self._setup_database('ephemerides')
        self.photometry_loaded = self._setup_database('photometry')
        self.catalogues_loaded = self._setup_database('catalogues')
        self.ut_loaded = self._setup_database('ut')

        self.leap_seconds_data = None
        self.earth_rotation_data = None
        self.barycenter_data = None
        self.sun_data = None
        self.moon_data = None
        self.ecc_data = None
        self.all_filters_data = None

    def exotethys(self):
        return self.exotethys_loaded

    def ephemeris(self):
        return self.ephemerides_loaded

    def photometry(self):
        return self.photometry_loaded

    def catalogues(self):
        return self.catalogues_loaded

    def ut(self):
        return self.ut_loaded

    def all_filters(self):

        if not self.all_filters_data:
            self.all_filters_data = [os.path.split(ff)[1].split('.')[0]
                                     for ff in glob.glob(os.path.join(self.photometry(), '*'))]

        return self.all_filters_data

    def ecc(self):

        if not self.ecc_data:
            stars = open_dict(os.path.join(self.catalogues(), 'ecc_stars.pickle'))
            planets = open_dict(os.path.join(self.catalogues(), 'ecc_planets.pickle'))

            hosts = {stars[ff]['simbad_id']: stars[ff]['planets'] for ff in stars}

            def _flat_name(name):

                flat_name_list = [
                    [' ', ''],
                    ['-', ''],
                    ['cancri', 'cnc'],
                    ['hatp10', 'wasp11'],
                    ['wasp40', 'hatp27'],
                    ['wasp51', 'hatp30'],
                    ['wasp86', 'kelt12'],
                    ['kelt22', 'wasp173'],
                ]

                name = name.lower()

                for char in flat_name_list:
                    name = name.replace(char[0], char[1])

                return name

            flats = {_flat_name(ff): stars[ff]['simbad_id'] for ff in stars}
            for planet in planets:
                star = planets[planet]['star']
                flats[_flat_name(planet)] = star
                flats[_flat_name(star)] = star

            self.ecc_data = {'stars': stars, 'planets': planets, 'hosts': hosts, 'flats': flats}

        return self.ecc_data

    def barycentre(self, utc_jd):

        if not self.barycenter_data:
            self.barycenter_data = open_dict(os.path.join(self.ephemeris(), 'bjd_dict.pickle'))

        bjd_dict = self.barycenter_data[int(utc_jd)]

        ssb_t = bjd_dict['t']
        ssb_ra = bjd_dict['ra']
        ssb_dec = bjd_dict['dec']
        ssb_d = bjd_dict['d']
        ssb_dt = bjd_dict['dt']

        ssb_ra = interp1d(ssb_t, ssb_ra, kind='cubic')(utc_jd)
        ssb_dec = interp1d(ssb_t, ssb_dec, kind='cubic')(utc_jd)
        ssb_d = interp1d(ssb_t, ssb_d, kind='cubic')(utc_jd)
        ssb_dt = interp1d(ssb_t, ssb_dt, kind='cubic')(utc_jd)

        return ssb_ra, ssb_dec, ssb_d, ssb_dt

    def heliocentre(self, utc_jd):

        if not self.sun_data:
            self.sun_data = open_dict(os.path.join(self.ephemeris(), 'hjd_dict.pickle'))

        hjd_dict = self.sun_data[int(utc_jd)]

        ssb_t = hjd_dict['t']
        ssb_ra = hjd_dict['ra']
        ssb_dec = hjd_dict['dec']
        ssb_d = hjd_dict['d']
        ssb_dt = hjd_dict['dt']

        ssb_ra = interp1d(ssb_t, ssb_ra, kind='cubic')(utc_jd)
        ssb_dec = interp1d(ssb_t, ssb_dec, kind='cubic')(utc_jd)
        ssb_d = interp1d(ssb_t, ssb_d, kind='cubic')(utc_jd)
        ssb_dt = interp1d(ssb_t, ssb_dt, kind='cubic')(utc_jd)

        return ssb_ra, ssb_dec, ssb_d, ssb_dt

    def mooncentre(self, utc_jd):

        if not self.moon_data:
            self.moon_data = open_dict(os.path.join(self.ephemeris(), 'moon_dict.pickle'))

        moon_dict = self.moon_data[int(utc_jd)]

        ssb_t = moon_dict['t']
        ssb_ra = moon_dict['ra']
        ssb_dec = moon_dict['dec']
        ssb_d = moon_dict['d']
        ssb_dt = moon_dict['dt']

        ssb_ra = interp1d(ssb_t, ssb_ra, kind='cubic')(utc_jd)
        ssb_dec = interp1d(ssb_t, ssb_dec, kind='cubic')(utc_jd)
        ssb_d = interp1d(ssb_t, ssb_d, kind='cubic')(utc_jd)
        ssb_dt = interp1d(ssb_t, ssb_dt, kind='cubic')(utc_jd)

        return ssb_ra, ssb_dec, ssb_d, ssb_dt

    def leap_seconds(self, utc):

        if not self.leap_seconds_data:
            self.leap_seconds_data = [[datetime.datetime(int(line[3]), int(line[2]), int(line[1])), line[4]] for line in
                                      np.loadtxt(open(os.path.join(self.ut(), 'leap_seconds.txt')))]

        ls = self.leap_seconds_data[-1][1]

        if utc < self.leap_seconds_data[0][0]:
            print('Conversion to TT is not valid before {0}'.format(self.leap_seconds_data[0][0].isoformat()))
            ls = self.leap_seconds_data[0][1]
        else:
            for check in range(1, len(self.leap_seconds_data)):
                if utc < self.leap_seconds_data[check][0]:
                    ls = self.leap_seconds_data[check - 1][1]
                    break

        return ls

    def earth_rotation(self, utc_jd):

        if not self.earth_rotation_data:

            earth_rotation_data_x = []
            earth_rotation_data_y = []
            for ff in open(os.path.join(self.ut(), 'earth_rotation.txt')).readlines():
                if ff[58:68].replace(' ', '') != '':
                    earth_rotation_data_x.append(float(ff[7:12]) + 2400000.5)
                    earth_rotation_data_y.append(float(ff[58:68].replace(' ', '')))

            self.earth_rotation_data = interp1d(earth_rotation_data_x, earth_rotation_data_y, kind='cubic')

        return float(self.earth_rotation_data(utc_jd))

    def _setup_database(self, database_name):

        print('Checking {0} database...'.format(database_name))

        # define paths

        database_directory_path = os.path.join(self.databases_directory_path, database_name)
        database_file_path = os.path.join(self.databases_directory_path, database_name + '.pickle')
        database_link_file_path = os.path.join(self.databases_directory_path, database_name + '_link.txt')
        database_file_path_new = os.path.join(self.databases_directory_path, database_name + '_new.pickle')
        database_file_path_old = os.path.join(self.databases_directory_path, database_name + '_old.pickle')
        last_update_file_path = os.path.join(self.databases_directory_path, '{0}_last_update.txt'.format(database_name))

        # define paths

        # check if everything exists, if not reset database

        if not os.path.isdir(database_directory_path) or not os.path.isfile(database_file_path) or not os.path.isfile(database_link_file_path):

            try:
                shutil.rmtree(database_directory_path)
            except:
                pass

            try:
                os.remove(database_file_path)
            except:
                pass

            try:
                os.remove(database_file_path_old)
            except:
                pass

            try:
                os.remove(database_file_path_new)
            except:
                pass

            try:
                os.remove(database_link_file_path)
            except:
                pass

            try:
                os.remove(last_update_file_path)
            except:
                pass

            os.mkdir(database_directory_path)

            if not download(self.databases[self.version][database_name], database_file_path):
                print('\n{0} features cannot be used.'.format(database_name))
                return False
            else:
                shutil.copy(database_file_path, database_file_path_old)
                w = open(database_link_file_path, 'w')
                w.write(self.databases[self.version][database_name])
                w.close()

                try:
                    new_database = open_dict(database_file_path)
                    download(new_database['zipfile'], database_directory_path + '.zip')
                    new_database = zipfile.ZipFile(database_directory_path + '.zip', 'r')
                    here = os.path.abspath('.')
                    os.chdir(self.databases_directory_path)
                    new_database.extractall()
                    os.chdir(here)
                    os.remove(database_directory_path + '.zip')
                except:
                    pass

        # check if everything exists, if not reset database

        # download database if there is an update

        if self.databases[self.version][database_name] != open(database_link_file_path).read():

            if not download(self.databases[self.version][database_name], database_file_path_new):
                pass
            else:
                shutil.move(database_file_path, database_file_path_old)
                shutil.move(database_file_path_new, database_file_path)

                w = open(database_link_file_path, 'w')
                w.write(self.databases[self.version][database_name])
                w.close()

        # download database if there is an update

        # check all files in database, remove files that need to be updated

        current_database = open_dict(database_file_path_old)
        new_database = open_dict(database_file_path)

        for dbx_file in current_database['files']:

            if dbx_file not in new_database['files']:
                try:
                    os.remove(os.path.join(self.databases_directory_path,
                                           new_database['files'][dbx_file]['local_path']))
                except:
                    pass
            elif new_database['files'][dbx_file]['link'] != current_database['files'][dbx_file]['link']:
                try:
                    os.remove(os.path.join(self.databases_directory_path,
                                           new_database['files'][dbx_file]['local_path']))
                except:
                    pass

        # check for updates, remove files that need to be updated

        # download missing files

        final_check = True

        for dbx_file in new_database['files']:
            if not os.path.isfile(os.path.join(self.databases_directory_path,
                                               new_database['files'][dbx_file]['local_path'])):
                try:
                    os.remove(last_update_file_path)
                except:
                    pass
                if not download(new_database['files'][dbx_file]['link'],
                                os.path.join(self.databases_directory_path,
                                             new_database['files'][dbx_file]['local_path'])):
                    final_check = False

        # download missing files

        # update files from external links

        frequency = new_database['frequency']
        if frequency:

            try:
                last_update_date = int(open(last_update_file_path).read())
            except:
                last_update_date = 0

            today = int(time.strftime('%y%m%d'))

            if today >= last_update_date + frequency:

                for dbx_file in new_database['files']:
                    if 'external_link' in new_database['files'][dbx_file]:
                        print('\tUpdating: ', dbx_file)
                        if not download(new_database['files'][dbx_file]['external_link'],
                                        os.path.join(self.databases_directory_path,
                                                     new_database['files'][dbx_file]['local_path'])):
                            final_check = False

                w = open(last_update_file_path, 'w')
                w.write(time.strftime('%y%m%d'))
                w.close()

        # update files from external links

        if not final_check:
            print('\n{0} features cannot be used.'.format(database_name))
            return False
        else:
            if current_database != new_database:
                shutil.copy(database_file_path, database_file_path_old)
            return database_directory_path


plc_data = PlcData()
