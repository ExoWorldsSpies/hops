__all__ = ['DTime', 'UTC', 'JD', 'MJD', 'now']

import datetime
from functools import lru_cache

from pylightcurve.errors import *
from pylightcurve.__databases__ import plc_data
from pylightcurve.spacetime.angles import *


class _Time:

    def __init__(self, date):

        self.utc = date

        self._definition = '{0}'.format(date)

    @lru_cache()
    def iso(self):
        return self.utc.isoformat()

    @lru_cache()
    def jd(self):
        time_dif = self.utc - datetime.datetime(1600, 1, 1, 0, 0, 0, 0)
        return time_dif.days + (time_dif.seconds + time_dif.microseconds / 1000000.0) / 60.0 / 60.0 / 24.0 + 2305447.5

    @lru_cache()
    def mjd(self):
        return self.jd() - 2400000.5

    @lru_cache()
    def gmst(self):

        ut1 = self.jd() + plc_data.earth_rotation(self.jd()) / 60.0 / 60.0 / 24.0

        tu = ut1 - 2451545.0

        tt = self.jd() + (32.184 + plc_data.leap_seconds(self.utc)) / 60.0 / 60.0 / 24.0

        t = (tt - 2451545.0) / 36525.0

        gmst = (1296000 * (0.7790572732640 + 1.00273781191135448 * tu)
                + 0.014506
                + 4612.15739966 * t
                + 1.39667721 * (t ** 2)
                - 0.00009344 * (t ** 3)
                + 0.00001882 * (t ** 4)
                )

        return Degrees(0, 0, gmst)

    def __add__(self, dt):
        if isinstance(dt, DTime):
            return UTC(self.utc + dt.dtime)
        else:
            raise PyLCInputError('Only a plc.DTime object can be added to a {0} object.'.format(self._get_class()))

    def __sub__(self, dt):
        if isinstance(dt, DTime):
            return UTC(self.utc - dt.dtime)
        else:
            raise PyLCInputError('Only a plc.DTime object can be subtracted from a {0} object.'.format(
                self._get_class()))

    @lru_cache()
    def _get_class(self):
        return 'plc.{0}'.format(str(self.__class__).split('.')[-1][:-2])

    def __repr__(self):
        return '{0}({1} UTC, defined as {2})'.format(self._get_class(), self.iso(), self._definition)


class DTime:

    def __init__(self, days=0.0, hours=0.0, minutes=0.0, seconds=0.0):
        try:
            days, hours, minutes, seconds = float(days), float(hours), float(minutes), float(seconds)
            if days >= 0 and hours >= 0 and minutes >= 0 and seconds >= 0:
                self.dtime = datetime.timedelta(seconds=seconds + minutes * 60.0 + hours * 3600.0 + days * 86400.0)
            else:
                raise PyLCInputError('Only positive time intervals are allowed.')
        except:
            raise PyLCInputError('Only positive time intervals are allowed.')

    def __repr__(self):
        return 'plc.DTime({0} seconds)'.format(self.dtime.total_seconds())


class UTC(_Time):

    def __init__(self, date):

        if isinstance(date, datetime.datetime):

            utc = date

        elif isinstance(date, str):

            try:

                date = date.replace('T', ' ')

                year, month, day = list(map(float, date.split(' ')[0].split('-')))
                hour_minutes_seconds = date.split(' ')[1]

                if '+' in hour_minutes_seconds or '-' in hour_minutes_seconds:
                    if '+' in hour_minutes_seconds:
                        time_zone_sign = 1.0
                        hour_minutes_seconds, time_zone = hour_minutes_seconds.split('+')
                    else:
                        time_zone_sign = -1.0
                        hour_minutes_seconds, time_zone = hour_minutes_seconds.split('-')

                    if ':' in time_zone:
                        time_zone_hours = float(time_zone.split(':')[0])
                        time_zone_minutes = float(time_zone.split(':')[1])
                    else:
                        time_zone_hours = float(time_zone[:2])
                        time_zone_minutes = float(time_zone[2:])

                    time_zone = time_zone_sign * (time_zone_hours + time_zone_minutes / 60)

                else:
                    time_zone = 0

                hour, minutes, seconds = hour_minutes_seconds.split(':')
                hour, minutes, seconds = [float(ff) for ff in [hour, minutes, seconds]]

                if (int(year) != year or int(month) != month or int(day) != day or int(hour) != hour or
                        int(minutes) != minutes):
                    raise PyLCInputError('Not valid date format')

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

        _Time.__init__(self, utc)

        self._definition = '{0} UTC'.format(utc)


class JD(_Time):

    def __init__(self, date):

        try:
            jd = float(date)
        except:
            raise PyLCInputError('Not valid date format')

        if jd < 2305447.5:
            raise PyLCInputError('Not valid date format')

        _Time.__init__(self, datetime.datetime(1600, 1, 1, 0, 0, 0, 0) + datetime.timedelta(
                        days=(jd - 2305447.5)))

        self._definition = '{0} JD'.format(jd)


class MJD(_Time):

    def __init__(self, date):

        try:
            mjd = float(date)
        except:
            raise PyLCInputError('Not valid date format')

        if mjd < -94553.0:
            raise PyLCInputError('Not valid date format')

        _Time.__init__(self, datetime.datetime(1600, 1, 1, 0, 0, 0, 0) + datetime.timedelta(
                    days=(mjd - 2305447.5 + 2400000.5)))

        self._definition = '{0} MJD'.format(mjd)


def _is_time(item):
    return isinstance(item, UTC) or isinstance(item, JD) or isinstance(item, MJD)


def _request_time(item):
    if _is_time(item):
        pass
    else:
        raise PyLCInputError('A time object is required (plc.UTC, plc.JD or plc.MJD)')


def now():
    return UTC(datetime.datetime.now(datetime.timezone.utc).isoformat())
