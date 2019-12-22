from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import datetime
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from functools import lru_cache

from ._0errors import *
from ._1databases import databases, plc_data

# Angles


def _break_min_sec(value):

    a = int(value)
    b = (float(value) - a) * 60.0
    c = int(b)
    d = (b - c) * 60.0
    return a, c, d


def _collapse_min_sec(degrees_hours, minutes, seconds):

        try:
            degrees_hours = float(degrees_hours)
        except:

            if isinstance(degrees_hours, str):

                if minutes != 0 or seconds != 0:
                    raise PyLCInputError('Not valid angle format')
                else:
                    try:
                        degrees_hours, minutes, seconds = degrees_hours.replace(':', ' ').split()
                    except:
                        raise PyLCInputError('Not valid angle format')

        try:

            degrees_hours = float(degrees_hours)
            minutes = float(minutes)
            seconds = float(seconds)

            if degrees_hours >= 0:
                sign = 1.0
            else:
                sign = -1.0

        except:
            raise PyLCInputError('Not valid angle format')

        if minutes < 0 or seconds < 0:
            raise PyLCInputError('Not valid angle format. You should use positive values for minutes and seconds.')

        return sign * (abs(degrees_hours) + minutes / 60.0 + seconds / 3600.0)


class _DMS:

    def __init__(self, degrees):
            self.d, self.m, self.s = _break_min_sec(degrees)

    def __str__(self):
        return 'dms(d = {0}, m = {1}, s = {2})'.format(self.d, self.m, self.s)

    def __repr__(self):
        return self.__str__()


class _HMS:

    def __init__(self, hours):
        self.h, self.m, self.s = _break_min_sec(hours)

    def __str__(self):
        return 'hms(h = {0}, m = {1}, s = {2})'.format(self.h, self.m, self.s)

    def __repr__(self):
        return self.__str__()


class _Angle:

    def __init__(self, angle):

        if angle < 0:
            angle += 360 * (int(abs(angle) / 360.0) + 1.0)
        angle -= 360.0 * int(angle / 360.0)

        self.deg = angle
        if angle > 180:
            self.deg_pm = - 360 + angle
        else:
            self.deg_pm = angle
        self.hours = angle / 15.0
        self.rad = angle * np.pi / 180

        self.dms = _DMS(self.deg)
        self.hms = _HMS(self.hours)

        self.sin = np.sin(self.rad)
        self.cos = np.cos(self.rad)
        self.tan = np.tan(self.rad)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return 'plc.Angle(dms: {0}:{1}:{2}, hms: {3}:{4}:{5}, rad:{6})'.format(
            self.dms.d, self.dms.m, self.dms.s, self.hms.h, self.hms.m, self.hms.s, self.rad)


class Degrees(_Angle):

    def __init__(self, degrees, minutes=0.0, seconds=0.0):

        angle = _collapse_min_sec(degrees, minutes, seconds)

        _Angle.__init__(self, angle)


class Hours(_Angle):

    def __init__(self, hours, minutes=0.0, seconds=0.0):

        hours = _collapse_min_sec(hours, minutes, seconds)

        _Angle.__init__(self, hours * 15.0)


class Rad(_Angle):

    def __init__(self, rad):

        try:
            rad = float(rad)
        except:
            raise PyLCInputError('Not valid angle format')

        _Angle.__init__(self, rad * 180.0 / np.pi)


def is_angle(obsject):
    if isinstance(obsject, Degrees) or isinstance(obsject, Hours) or isinstance(obsject, Rad):
        return True
    else:
        return False



# Observatory

class _Horizon:

    def __init__(self, horizon):

        self.horizon_value = None
        self.horizon_function = None

        try:
            self.horizon_value = float(horizon)
        except:
            pass

        if self.horizon_value is None:

            try:

                if isinstance(horizon, str):
                    horizon_list = []
                    for horizon_line in horizon.split('\n'):
                        if horizon_line != '':
                            horizon_list.append([Degrees(horizon_line.split()[0]), Degrees(horizon_line.split()[1])])

                elif isinstance(horizon, list):
                    horizon_list = [[Degrees(ff[0]), Degrees(ff[1])] for ff in horizon]

                elif isinstance(horizon, np.ndarray):
                    horizon_list = [[Degrees(ff[0]), Degrees(ff[1])] for ff in horizon]

                else:
                    raise PyLCInputError('Not valid horizon format')

                horizon_list = [[ff[0].deg, ff[1].deg] for ff in horizon_list]
                horizon_list.append([360.0, horizon_list[0][0]])

                horizon_list = np.swapaxes(np.array(horizon_list, dtype=float), 0, 1)
                self.horizon_function = interp1d(horizon_list[0], horizon_list[1])

            except:
                raise PyLCInputError('Not valid horizon format')

    def horizon(self, azimuth):

        if is_angle(azimuth):

            if self.horizon_value is None:
                return Degrees(self.horizon_function(azimuth.deg))

            else:
                return Degrees(self.horizon_value)

        else:
            raise PyLCInputError('Azimuth should be an Angle object')


class Observatory:

    def __init__(self, latitude, longitude, time_zone=0, horizon=0):

        if is_angle(latitude):
            self.latitude = latitude
        else:
            raise PyLCInputError('Latitude should be an Angle object')

        if 270 > self.latitude.deg > 90:
            raise PyLCInputError('Latitude must within -90, 90 degrees')

        if is_angle(longitude):
            self.longitude = longitude
        else:
            raise PyLCInputError('Longitude should be an Angle object')

        try:
            self.time_zone = float(time_zone)
            self.time_zone_shift = datetime.timedelta(days=(self.time_zone / 24))
        except:
            raise PyLCInputError('Not valid time zone')

        if self.time_zone > 12 or self.time_zone < -12:
            raise PyLCInputError('Time zone must within -12, 12 hours')

        self.horizon = _Horizon(horizon).horizon

    def __str__(self):
        return 'plc.Observatory(latitude = {0} deg, longitude = {1} deg, ' \
               'time zone = {2} hours)'.format(self.latitude.deg, self.longitude.deg, self.time_zone)

    def __repr__(self):
        return self.__str__()


# # Target

class Target:

    def __init__(self, ra, dec):

        if is_angle(dec):
            self.dec = dec
        else:
            raise PyLCInputError('Declination should be an Angle object')

        if 270 > self.dec.deg > 90:
            raise PyLCInputError('Declination must within -90, 90 degrees')

        if is_angle(ra):
            self.ra = ra
        else:
            raise PyLCInputError('Right Ascension should be an Angle object')

        if self.dec.deg > 180:
            sign = '-'
            dec_print = Degrees(360 - self.dec.deg)
        else:
            sign = '+'
            dec_print = self.dec

        self.coord = '{0}:{1}:{2}.{3} {4}{5}:{6}:{7}.{8}'.format(str(self.ra.hms.h).zfill(2),
                                                                 str(self.ra.hms.m).zfill(2),
                                                                 str(int(self.ra.hms.s)).zfill(2),
                                                                 str(round(self.ra.hms.s - int(self.ra.hms.s), 4))[2:],
                                                                 sign,
                                                                 str(dec_print.dms.d).zfill(2),
                                                                 str(dec_print.dms.m).zfill(2),
                                                                 str(int(dec_print.dms.s)).zfill(2),
                                                                 str(round(dec_print.dms.s - int(dec_print.dms.s), 4))[2:])

    def __str__(self):
        return 'plc.Target(RA = {0} deg, DEC = {1} deg'.format(self.ra.deg, self.dec.deg_pm)

    def __repr__(self):
        return self.__str__()


# Times

class _Time:

    def __init__(self, date, history=[]):

        self.utc = date

        # jd
        time_dif = self.utc - datetime.datetime(1600, 1, 1, 0, 0, 0, 0)
        self.jd = (time_dif.days + (time_dif.seconds + time_dif.microseconds / 1000000.0) / 60.0 / 60.0 / 24.0 +
                   2305447.5)

        # mjd
        self.mjd = self.jd - 2400000.5

        self.history = history
        self.history.append('Final UTC: {0}'.format(self.utc.isoformat()))

    @lru_cache()
    def tt(self):
        return self.utc + datetime.timedelta(seconds=(32.184 + plc_data.leap_seconds(self.utc)))

    @lru_cache()
    def ut1(self):
        return self.utc + datetime.timedelta(seconds=plc_data.earth_rotation(self.jd))

    @lru_cache()
    def gmst(self):

        tu = (self.ut1() - plc_data.tt_j2000).total_seconds() / 60.0 / 60.0 / 24.0

        t = (self.tt() - plc_data.tt_j2000).total_seconds() / 60.0 / 60.0 / 24.0 / 36525.0

        gmst = (1296000 * (0.7790572732640 + 1.00273781191135448 * tu)
                + 0.014506
                + 4612.15739966 * t
                + 1.39667721 * (t ** 2)
                - 0.00009344 * (t ** 3)
                + 0.00001882 * (t ** 4)
                )

        return Degrees(0, 0, gmst)

    @lru_cache()
    def get_barycetre_data(self):
        return plc_data.barycentre(self.jd)

    @lru_cache()
    def get_heliocetre_data(self):
        return plc_data.heliocentre(self.jd)

    @lru_cache()
    def get_sun(self):
        ssb_ra, ssb_dec, ssb_d, ssb_dt = self.get_heliocetre_data()
        return Target(Rad(ssb_ra), Rad(ssb_dec))

    def bjd_tdb(self, target):

        if not isinstance(target, Target):
            raise PyLCInputError('A plc.Target object needs to be passed.')

        ssb_ra, ssb_dec, ssb_d, ssb_dt = self.get_barycetre_data()

        a = ssb_d / 60.0 / 24.0
        b = target.dec.sin * np.sin(ssb_dec)
        c = target.dec.cos * np.cos(ssb_dec) * np.cos(target.ra.rad - ssb_ra)

        return self.jd - a * (b + c) + ssb_dt / 60.0 / 60.0 / 24.0

    def bjd_utc(self, target):

        if not isinstance(target, Target):
            raise PyLCInputError('A plc.Target object needs to be passed.')

        ssb_ra, ssb_dec, ssb_d, ssb_dt = self.get_barycetre_data()

        a = ssb_d / 60.0 / 24.0
        b = target.dec.sin * np.sin(ssb_dec)
        c = target.dec.cos * np.cos(ssb_dec) * np.cos(target.ra.rad - ssb_ra)

        return self.jd - a * (b + c)

    def hjd_tdb(self, target):

        if not isinstance(target, Target):
            raise PyLCInputError('A plc.Target object needs to be passed.')

        ssb_ra, ssb_dec, ssb_d, ssb_dt = self.get_heliocetre_data()

        a = ssb_d / 60.0 / 24.0
        b = target.dec.sin * np.sin(ssb_dec)
        c = target.dec.cos * np.cos(ssb_dec) * np.cos(target.ra.rad - ssb_ra)

        return self.jd - a * (b + c) + ssb_dt / 60.0 / 60.0 / 24.0

    def hjd_utc(self, target):

        if not isinstance(target, Target):
            raise PyLCInputError('A plc.Target object needs to be passed.')

        ssb_ra, ssb_dec, ssb_d, ssb_dt = self.get_heliocetre_data()

        a = ssb_d / 60.0 / 24.0
        b = target.dec.sin * np.sin(ssb_dec)
        c = target.dec.cos * np.cos(ssb_dec) * np.cos(target.ra.rad - ssb_ra)

        return self.jd - a * (b + c)

    def lt(self, observatory):

        if not isinstance(observatory, Observatory):
            raise PyLCInputError('A plc.Observatory object needs to be passed.')

        return self.utc + observatory.time_zone_shift

    def lst(self, observatory):

        if not isinstance(observatory, Observatory):
            raise PyLCInputError('A plc.Observatory object needs to be passed.')

        return Degrees(self.gmst().deg + observatory.longitude.deg)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return 'plc.Time({0} UTC)'.format(self.utc.isoformat().replace('T', ' '))


class NOW(_Time):
    def __init__(self):
        _Time.__init__(self, datetime.datetime.now())


class UTC(_Time):

    def __init__(self, date, history=[]):

        if isinstance(date, datetime.datetime):

            utc = date

        elif isinstance(date, str):

            try:

                date = date.replace('T', ' ')

                year, month, day = list(map(float, date.split(' ')[0].split('-')))
                hour, minutes, seconds = date.split(' ')[1].split(':')

                time_zone = 0
                if '+' in seconds or '-' in seconds:
                    if '+' in seconds:
                        time_zone_sign = 1.0
                        seconds, time_zone = seconds.split('+')
                    else:
                        time_zone_sign = -1.0
                        seconds, time_zone = seconds.split('-')

                    time_zone = time_zone.replace(':', '')
                    time_zone_hours = float(time_zone[:2])
                    time_zone_minutes = float(time_zone[2:])
                    time_zone = time_zone_sign * (time_zone_hours + time_zone_minutes / 60)

                hour, minutes, seconds = [float(ff) for ff in [hour, minutes, seconds]]

            except:
                raise PyLCInputError('Not valid date format')

            if (int(year) != year or int(month) != month or int(day) != day or int(hour) != hour or
                    int(minutes) != minutes):
                raise PyLCInputError('Not valid date format')

            try:
                year = int(year)
                month = int(month)
                day = int(day)
                hour = int(hour)
                minutes = int(minutes)
                microsec = int((seconds - int(seconds)) * 1000000)
                seconds = int(seconds)

                utc = (datetime.datetime(year, month, day, hour, minutes, seconds, microsec) -
                       datetime.timedelta(days=time_zone / 24))

            except:
                raise PyLCInputError('Not valid date format')

        else:
            raise PyLCInputError('Not valid date format')

        self.history = history

        _Time.__init__(self, utc, history=self.history)


class JD(_Time):

    def __init__(self, date, history=[]):

        try:
            jd = float(date)
        except:
            raise PyLCInputError('Not valid date format')

        self.history = history
        self.history.append('Converted to UTC from JD: {0}'.format(jd))

        _Time.__init__(self, datetime.datetime(1600, 1, 1, 0, 0, 0, 0) + datetime.timedelta(
                    days=(jd - 2305447.5)), history=self.history)


class MJD(_Time):

    def __init__(self, date):

        try:
            mjd = float(date)
        except:
            raise PyLCInputError('Not valid date format')

        self.history = []
        self.history.append('Converted to UTC from MJD: {0}'.format(mjd))

        _Time.__init__(self, datetime.datetime(1600, 1, 1, 0, 0, 0, 0) + datetime.timedelta(
            days=(mjd - 2305447.5 + 2400000.5)), history=self.history)


class BJDTDB(JD):

    def __init__(self, date, target):

        try:
            bjd_tdb = float(date)
        except:
            raise PyLCInputError('Not valid date format')

        if not isinstance(target, Target):
            raise PyLCInputError('A plc.Target object needs to be passed.')

        self.history = []
        self.history.append('Converted to UTC from BJD_TDB: {0}, with target at RA: {1} deg, DEC: {2} deg'.format(
            bjd_tdb, target.ra.deg, target.dec.deg))

        def _jd_utc_to_bjd_tdb(jd):

            ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.barycentre(jd)

            a = ssb_d / 60.0 / 24.0
            b = target.dec.sin * np.sin(ssb_dec)
            c = target.dec.cos * np.cos(ssb_dec) * np.cos(target.ra.rad - ssb_ra)

            return jd - a * (b + c) + ssb_dt / 60.0 / 60.0 / 24.0

        def func(x, jd):
            return _jd_utc_to_bjd_tdb(jd) - bjd_tdb

        JD.__init__(self, curve_fit(func, [0], [0], p0=[bjd_tdb])[0][0], history=self.history)


class BJDUTC(JD):

    def __init__(self, date, target):

        try:
            bjd_utc = float(date)
        except:
            raise PyLCInputError('Not valid date format')

        if not isinstance(target, Target):
            raise PyLCInputError('A plc.Target object needs to be passed.')

        self.history = []
        self.history.append('Converted to UTC from BJD_UTC: {0}, with target at RA: {1} deg, DEC: {2} deg'.format(
            bjd_utc, target.ra.deg, target.dec.deg))

        def _jd_utc_to_bjd_utc(jd):

            ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.barycentre(jd)

            a = ssb_d / 60.0 / 24.0
            b = target.dec.sin * np.sin(ssb_dec)
            c = target.dec.cos * np.cos(ssb_dec) * np.cos(target.ra.rad - ssb_ra)

            return jd - a * (b + c)

        def func(x, jd):
            return _jd_utc_to_bjd_utc(jd) - bjd_utc

        JD.__init__(self, curve_fit(func, [0], [0], p0=[bjd_utc])[0][0], history=self.history)


class HJDTDB(JD):

    def __init__(self, date, target):

        try:
            hjd_tdb = float(date)
        except:
            raise PyLCInputError('Not valid date format')

        if not isinstance(target, Target):
            raise PyLCInputError('A plc.Target object needs to be passed.')

        self.history = []
        self.history.append('Converted to UTC from HJD_TDB: {0}, with target at RA: {1} deg, DEC: {2} deg'.format(
            hjd_tdb, target.ra.deg, target.dec.deg))

        def _jd_utc_to_hjd_tdb(jd):

            ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.heliocentre(jd)

            a = ssb_d / 60.0 / 24.0
            b = target.dec.sin * np.sin(ssb_dec)
            c = target.dec.cos * np.cos(ssb_dec) * np.cos(target.ra.rad - ssb_ra)

            return jd - a * (b + c) + ssb_dt / 60.0 / 60.0 / 24.0

        def func(x, jd):
            return _jd_utc_to_hjd_tdb(jd) - hjd_tdb

        JD.__init__(self, curve_fit(func, [0], [0], p0=[hjd_tdb])[0][0], history=self.history)


class HJDUTC(JD):

    def __init__(self, date, target):

        try:
            hjd_utc = float(date)
        except:
            raise PyLCInputError('Not valid date format')

        if not isinstance(target, Target):
            raise PyLCInputError('A plc.Target object needs to be passed.')

        self.history = []
        self.history.append('Converted to UTC from HJD_UTC: {0}, with target at RA: {1} deg, DEC: {2} deg'.format(
            hjd_utc, target.ra.deg, target.dec.deg))

        def _jd_utc_to_hjd_utc(jd):

            ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.heliocentre(jd)

            a = ssb_d / 60.0 / 24.0
            b = target.dec.sin * np.sin(ssb_dec)
            c = target.dec.cos * np.cos(ssb_dec) * np.cos(target.ra.rad - ssb_ra)

            return jd - a * (b + c)

        def func(x, jd):
            return _jd_utc_to_hjd_utc(jd) - hjd_utc

        JD.__init__(self, curve_fit(func, [0], [0], p0=[hjd_utc])[0][0], history=self.history)


class LT(_Time):

    def __init__(self, date, observatory):

        try:
            lt = UTC(date).utc
        except:
            raise PyLCInputError('Not valid date format')

        if not isinstance(observatory, Observatory):
            raise PyLCInputError('A plc.Observatory object needs to be passed.')

        self.history = []
        self.history.append('Converted to JD from LT: {0}, with observatory time zone {1} hours'.format(
            lt.isoformat(), observatory.time_zone))

        _Time.__init__(self, lt - observatory.time_zone_shift, history=self.history)


def is_time(obsject):
    if (isinstance(obsject, NOW) or isinstance(obsject, UTC) or isinstance(obsject, JD) or isinstance(obsject, MJD)
       or isinstance(obsject, BJDTDB) or isinstance(obsject, BJDUTC) or isinstance(obsject, HJDTDB)
       or isinstance(obsject, HJDUTC) or isinstance(obsject, LT)):
        return True
    else:
        return False

# Observation

class Observation:

    def __init__(self, target, observatory):

        if isinstance(target, Target):
            self.target = target
        else:
            raise PyLCInputError('A plc.Target object needs to be passed.')

        if isinstance(observatory, Observatory):
            self.observatory = observatory
        else:
            raise PyLCInputError('A plc.Observatory object needs to be passed.')

        self.max_altitude = Degrees(90 - self.observatory.latitude.deg + self.target.dec.deg)
        self.min_altitude = Degrees(-(90 - self.observatory.latitude.deg - self.target.dec.deg))

    def azimuth_altitude(self, utc):
        ha = Degrees(utc.lst(self.observatory).deg - self.target.ra.deg)

        altitude = Rad(np.arcsin(np.clip(self.target.dec.sin * self.observatory.latitude.sin
                                         + self.target.dec.cos * self.observatory.latitude.cos * ha.cos, -1,
                                         1)))

        if ha.hours < 12:
            azimuth = Rad(np.pi - np.arccos(np.clip((self.target.dec.sin -
                                                     altitude.sin * self.observatory.latitude.sin) /
                                                    (altitude.cos * self.observatory.latitude.cos), -1, 1)))
        else:

            azimuth = Rad(np.pi + np.arccos(np.clip((self.target.dec.sin -
                                                     altitude.sin * self.observatory.latitude.sin) /
                                                    (altitude.cos * self.observatory.latitude.cos), -1, 1)))

        return azimuth, altitude

    def is_target_visible(self, utc):

        sidereal_time_horizon_function, sidereal_time_altitude_function = self._get_target_trajectory()

        target_horizon_diff = sidereal_time_horizon_function(utc.lst(self.observatory).hours)

        if target_horizon_diff > 0:
            return True
        else:
            return False

    def rise_set_events(self, utc0, utc1):

        st_0 = utc0.lst(self.observatory).hours
        jd_0 = utc0.jd
        max_dt = (utc1.jd - utc0.jd) * 24 * (24 / 23.9344696)

        all_target_events = []

        target_events = self._get_target_rise_set_st()

        for target_event in target_events:

            if target_event[0] > st_0:
                dt = target_event[0] - st_0
            else:
                dt = 24 + target_event[0] - st_0

            while dt <= max_dt:
                all_target_events.append([JD(jd_0 + dt * (23.9344696 / 24) / 24), target_event[1]])
                dt += 24

        all_target_events = [[ff[0].jd, ff[0], ff[1]] for ff in all_target_events]
        all_target_events.sort()
        all_target_events = [[ff[1], ff[2]] for ff in all_target_events]

        return all_target_events

    def sun_azimuth_altitude(self, utc):

        sun = utc.get_sun()

        ha = Degrees(utc.lst(self.observatory).deg - sun.ra.deg)

        altitude = Rad(np.arcsin(np.clip(sun.dec.sin * self.observatory.latitude.sin
                                         + sun.dec.cos * self.observatory.latitude.cos * ha.cos, -1,
                                         1)))

        if ha.hours < 12:
            azimuth = Rad(np.pi - np.arccos(np.clip((sun.dec.sin -
                                                     altitude.sin * self.observatory.latitude.sin) /
                                                    (altitude.cos * self.observatory.latitude.cos), -1, 1)))
        else:

            azimuth = Rad(np.pi + np.arccos(np.clip((sun.dec.sin -
                                                     altitude.sin * self.observatory.latitude.sin) /
                                                    (altitude.cos * self.observatory.latitude.cos), -1, 1)))

        return azimuth, altitude

    @lru_cache()
    def _get_target_trajectory(self):

        sidereal_time_list = list(np.arange(0, 24.001, 0.25))
        altitude_list = []
        horizon_list = []

        for sidereal_time in sidereal_time_list:

            ha = Hours(sidereal_time - self.target.ra.hours)

            altitude = Rad(np.arcsin(np.clip(self.target.dec.sin * self.observatory.latitude.sin
                                             + self.target.dec.cos * self.observatory.latitude.cos * ha.cos, -1,
                                             1)))
            if ha.hours < 12:
                azimuth = Rad(np.pi - np.arccos(np.clip((self.target.dec.sin -
                                                         altitude.sin * self.observatory.latitude.sin) /
                                                        (altitude.cos * self.observatory.latitude.cos), -1, 1)))
            else:

                azimuth = Rad(np.pi + np.arccos(np.clip((self.target.dec.sin -
                                                         altitude.sin * self.observatory.latitude.sin) /
                                                        (altitude.cos * self.observatory.latitude.cos), -1, 1)))

            altitude_list.append(altitude.deg_pm)
            horizon_list.append(self.observatory.horizon(azimuth).deg_pm)

        sidereal_time_horizon_function = interp1d(np.array(sidereal_time_list),
                                                  np.array(altitude_list) - np.array(horizon_list))
        sidereal_time_altitude_function = interp1d(np.array(sidereal_time_list), np.array(altitude_list))

        return sidereal_time_horizon_function, sidereal_time_altitude_function

    @lru_cache()
    def _get_target_rise_set_st(self):

        events = []

        sidereal_time_horizon_function, sidereal_time_altitude_function = self._get_target_trajectory()

        sidereal_time_list = list(np.arange(0, 24.001, 0.1))

        def target_horizon_diference(x, st):
            xx = Hours(st)
            return sidereal_time_horizon_function(xx.hours)

        test_alt = sidereal_time_horizon_function(sidereal_time_list[0])

        for sidereal_time in sidereal_time_list[1:]:
            new_test_alt = sidereal_time_horizon_function(sidereal_time)
            if test_alt * new_test_alt < 0:
                popt, pcov = curve_fit(target_horizon_diference, [0], [0], p0=[sidereal_time - 0.1])
                if test_alt < new_test_alt:
                    events.append([popt[0], 'target_rise'])
                else:
                    events.append([popt[0], 'target_set'])
            else:
                pass

            test_alt = new_test_alt

        return events

    def __str__(self):
        return 'plc.Observation(Target: RA = {0} deg, DEC = {1} deg,  Observatory: latitude = {2} deg, ' \
               'longitude = {3} deg, time zone = {4} hours)'.format(self.target.ra.deg, self.target.dec.deg_pm,
                                                                    self.observatory.latitude.deg,
                                                                    self.observatory.longitude.deg,
                                                                    self.observatory.time_zone)

    def __repr__(self):
        return self.__str__()

