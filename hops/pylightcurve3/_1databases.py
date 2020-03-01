from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .tools_databases import *
from .tools_files import *

from scipy.interpolate import interp1d


class Databases:

    def __init__(self, update=False):

        print('\n')

        self.clablimb = Database('clablimb', 'pylightcurve', __file__, date_to_update='200228',
                                 update=update).path
        self.ephemeris = Database('ephemeris', 'pylightcurve', __file__, date_to_update='200228',
                                  expire_date='200628', frequencey=30, update=update).path
        self.oec = Database('oec', 'pylightcurve', __file__, date_to_update='200228', frequencey=30,
                            update=update).path

        print('\n')

    def phoenix(self):
        return Database('phoenixplc3', 'pylightcurve', __file__, date_to_update='181213', ask_size='3.5GB').path


databases = Databases()


class Data:

    def __init__(self):

        self.leap_seconds_data = None
        self.earth_rotation_data = None
        self.tt_j2000 = datetime.datetime(2000, 1, 1, 12, 0, 0, 0)
        self.bjd_dict = None
        self.hjd_dict = None

    def leap_seconds(self, utc):

        if not self.leap_seconds_data:
            self.leap_seconds_data = [[datetime.datetime(int(line[3]), int(line[2]), int(line[1])), line[4]] for line in
                                      np.loadtxt(open(os.path.join(databases.ephemeris, 'leap_seconds.txt')))]

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

            self.earth_rotation_data = {}
            for ff in open(os.path.join(databases.ephemeris, 'earth_rotation.txt')).readlines():
                if ff[58:68].replace(' ', '') != '':
                    self.earth_rotation_data[int(int(ff[7:12]) + 2400000.5)] = float(ff[58:68].replace(' ', ''))

        return self.earth_rotation_data[int(utc_jd)]

    def barycentre(self, utc_jd):

        if not self.bjd_dict:
            self.bjd_dict = open_dict(glob.glob(os.path.join(databases.ephemeris, 'bjd_dict.pickle''*'))[0])

        bjd_dict = self.bjd_dict[int(utc_jd)]

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

        if not self.hjd_dict:
            self.hjd_dict = open_dict(glob.glob(os.path.join(databases.ephemeris, 'hjd_dict.pickle''*'))[0])

        hjd_dict = self.hjd_dict[int(utc_jd)]

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

plc_data = Data()
