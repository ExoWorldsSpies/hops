
__all__ = ['FixedTarget', 'Moon', 'Sun']

import numpy as np
from functools import lru_cache

from ..errors import *
from ..__databases__ import plc_data
from ..analysis.curve_fit import curve_fit
from .angles import *
from .angles import _request_angle
from ..spacetime.times import _request_time


class _Target:

    def __init__(self, ra, dec, name=None, all_names=[], epoch=None):

        self._notes = None
        self.ra, self.dec, self.name, self.all_names, self.epoch = ra, dec, name, all_names, epoch

        _request_angle(ra)
        _request_angle(dec)

        if 270 > dec.deg() > 90:
            raise PyLCInputError('Declination must be within -90, 90 degrees')

        if self.all_names != []:

            test = list(set(self.all_names).intersection(list(plc_data.ecc()['hosts'])))

            if len(test) > 0:
                self._notes = '/ Host of: ' + ', '.join(plc_data.ecc()['hosts'][test[0]])

    @lru_cache()
    def _get_class(self):
        return 'plc.{0}'.format(str(self.__class__).split('.')[-1][:-2])

    def reset(self, epoch):
        pass

    def __repr__(self):
        string = '{0}({1}'.format(self._get_class(), self.coord())
        if self.name:
            string += ', {0}'.format(self.name)
        if self._notes:
            string += ', {0}'.format(self._notes)
        return string + ')'

    @lru_cache()
    def coord(self):
        return '{0} {1}'.format(self.ra.hms(), self.dec.dms_coord())

    def distance_on_sphere(self, other):

        _request_target(other)

        return arccos(self.dec.sin() * other.dec.sin() +
                      self.dec.cos() * other.dec.cos() * (self.ra - other.ra).cos())


class FixedTarget(_Target):

    def __init__(self, *args, **kwargs):

        _Target.__init__(self, *args, **kwargs)

    def convert_to_bjd_tdb(self, time, time_format):

        if isinstance(time, float):

            if time_format in ['BJD_TDB', 'BJD_TT']:
                return time
            elif time_format == 'JD_UTC':
                return self._bjd_tdb(None, time)
            elif time_format == 'MJD_UTC':
                return self._bjd_tdb(None, time + 2400000.5)
            elif time_format == 'BJD_UTC':
                return self._bjd_tdb(None, curve_fit(self._bjd_utc, [0], [time], p0=[time])[0][0])
            elif time_format in ['HJD_TDB', 'HJD_TT']:
                return self._bjd_tdb(None, curve_fit(self._hjd_tdb, [0], [time], p0=[time])[0][0])
            elif time_format == 'HJD_UTC':
                return self._bjd_tdb(None, curve_fit(self._hjd_utc, [0], [time], p0=[time])[0][0])
            else:
                raise PyLCInputError(
                    'Not valid time format. Available formats: JD_UTC, MJD_UTC, BJD_UTC, BJD_TDB, BJD_TT, '
                    'HJD_UTC, HJD_BJD, HJD_TT')

        else:
            try:
                time = np.array(time, dtype=float)

                if time_format in ['BJD_TDB', 'BJD_TT']:
                    return time
                elif time_format == 'JD_UTC':
                    return np.array([self._bjd_tdb(None, ff) for ff in time])
                elif time_format == 'MJD_UTC':
                    return np.array([self._bjd_tdb(None, ff + 2400000.5) for ff in time])
                elif time_format == 'BJD_UTC':
                    return np.array([self._bjd_tdb(None, curve_fit(self._bjd_utc, [0], [ff], p0=[ff])[0][0])
                                     for ff in time])
                elif time_format in ['HJD_TDB', 'HJD_TT']:
                    return np.array([self._bjd_tdb(None, curve_fit(self._hjd_tdb, [0], [ff], p0=[ff])[0][0])
                                     for ff in time])
                elif time_format == 'HJD_UTC':
                    return np.array([self._bjd_tdb(None, curve_fit(self._hjd_utc, [0], [ff], p0=[ff])[0][0])
                                     for ff in time])
                else:
                    raise PyLCInputError(
                        'Not valid time format. Available formats: JD_UTC, MJD_UTC, BJD_UTC, BJD_TDB, BJD_TT, '
                        'HJD_UTC, HJD_BJD, HJD_TT')

            except:
                raise PyLCInputError('Not valid input for time')

    def convert_to_jd_utc(self, time, time_format):

        if isinstance(time, float):

            if time_format in ['BJD_TDB', 'BJD_TT']:
                return curve_fit(self._bjd_tdb, [0], [time], p0=[time])[0][0]
            elif time_format == 'JD_UTC':
                return time
            elif time_format == 'MJD_UTC':
                return time + 2400000.5
            elif time_format == 'BJD_UTC':
                return curve_fit(self._bjd_utc, [0], [time], p0=[time])[0][0]
            elif time_format in ['HJD_TDB', 'HJD_TT']:
                return curve_fit(self._hjd_tdb, [0], [time], p0=[time])[0][0]
            elif time_format == 'HJD_UTC':
                return curve_fit(self._hjd_utc, [0], [time], p0=[time])[0][0]
            else:
                raise PyLCInputError(
                    'Not valid time format. Available formats: JD_UTC, MJD_UTC, BJD_UTC, BJD_TDB, BJD_TT, '
                    'HJD_UTC, HJD_BJD, HJD_TT')

        else:
            try:
                time = np.array(time, dtype=float)

                if time_format in ['BJD_TDB', 'BJD_TT']:
                    return np.array([curve_fit(self._bjd_tdb, [0], [ff], p0=[time])[0][0] for ff in time])
                elif time_format == 'JD_UTC':
                    return time
                elif time_format == 'MJD_UTC':
                    return time + 2400000.5
                elif time_format == 'BJD_UTC':
                    return np.array([curve_fit(self._bjd_utc, [0], [ff], p0=[time])[0][0] for ff in time])
                elif time_format in ['HJD_TDB', 'HJD_TT']:
                    return np.array([curve_fit(self._hjd_tdb, [0], [ff], p0=[time])[0][0] for ff in time])
                elif time_format == 'HJD_UTC':
                    return np.array([curve_fit(self._hjd_utc, [0], [ff], p0=[time])[0][0] for ff in time])
                else:
                    raise PyLCInputError(
                        'Not valid time format. Available formats: JD_UTC, MJD_UTC, BJD_UTC, BJD_TDB, BJD_TT, '
                        'HJD_UTC, HJD_BJD, HJD_TT')

            except:
                raise PyLCInputError('Not valid input for time')

    def _hjd_utc(self, x, jd):

        ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.heliocentre(jd)

        a = ssb_d / 60.0 / 24.0
        b = self.dec.sin() * np.sin(ssb_dec)
        c = self.dec.cos() * np.cos(ssb_dec) * np.cos(self.ra.rad() - ssb_ra)

        return jd - a * (b + c)

    def _hjd_tdb(self, x, jd):

        ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.heliocentre(jd)

        a = ssb_d / 60.0 / 24.0
        b = self.dec.sin() * np.sin(ssb_dec)
        c = self.dec.cos() * np.cos(ssb_dec) * np.cos(self.ra.rad() - ssb_ra)

        return jd - a * (b + c) + ssb_dt / 60.0 / 60.0 / 24.0

    def _bjd_utc(self, x, jd):

        ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.barycentre(jd)

        a = ssb_d / 60.0 / 24.0
        b = self.dec.sin() * np.sin(ssb_dec)
        c = self.dec.cos() * np.cos(ssb_dec) * np.cos(self.ra.rad() - ssb_ra)

        return jd - a * (b + c)

    def _bjd_tdb(self, x, jd):

        ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.barycentre(jd)

        a = ssb_d / 60.0 / 24.0
        b = self.dec.sin() * np.sin(ssb_dec)
        c = self.dec.cos() * np.cos(ssb_dec) * np.cos(self.ra.rad() - ssb_ra)

        return jd - a * (b + c) + ssb_dt / 60.0 / 60.0 / 24.0


class Sun(_Target):

    def __init__(self, epoch):

        _request_time(epoch)

        ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.heliocentre(epoch.jd())

        _Target.__init__(self, Rad(ssb_ra), Rad(ssb_dec), epoch=epoch)

    def reset(self, epoch):

        _request_time(epoch)

        if epoch.jd() != self.epoch.jd():

            ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.heliocentre(epoch.jd())

            _Target.__init__(self, Rad(ssb_ra), Rad(ssb_dec), epoch=epoch)

    def __repr__(self):
        return 'plc.Sun(RA(hms)/DEC(dms): {0} at {1} UTC)'.format(self.coord(), self.epoch.utc)


class Moon(_Target):

    def __init__(self, epoch):

        _request_time(epoch)

        ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.mooncentre(epoch.jd())

        _Target.__init__(self, Rad(ssb_ra), Rad(ssb_dec), epoch=epoch)

    def reset(self, epoch):

        _request_time(epoch)

        if epoch.jd() != self.epoch.jd():

            ssb_ra, ssb_dec, ssb_d, ssb_dt = plc_data.mooncentre(epoch.jd())

            _Target.__init__(self, Rad(ssb_ra), Rad(ssb_dec), epoch=epoch)

    def illumination(self):

        sun_ra, sun_dec, sun_d, sun_dt = plc_data.heliocentre(self.epoch.jd())
        moon_ra, moon_dec, moon_d, moon_dt = plc_data.mooncentre(self.epoch.jd())

        theta = self.distance_on_sphere(FixedTarget(Rad(sun_ra), Rad(sun_dec)))

        sun_moon_d = np.sqrt(sun_d ** 2 + moon_d ** 2 - 2 * sun_d * moon_d * theta.cos())
        f = arcsin(theta.sin() * (sun_d / sun_moon_d))

        if theta.deg() < 90:
            return 0.5 * (1.0 - f.cos())
        else:
            return 0.5 * (1.0 + f.cos())

    def __repr__(self):
        return 'plc.Moon(RA(hms)/DEC(dms): {0} at {1} UTC)'.format(self.coord(), self.epoch.utc)


def _is_target(item):
    return isinstance(item, FixedTarget) or isinstance(item, Sun) or isinstance(item, Moon)


def _request_target(item):
    if _is_target(item):
        pass
    else:
        raise PyLCInputError('A plc.Target object is required')
