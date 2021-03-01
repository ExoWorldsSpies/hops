
import os
import glob
import time
import datetime
import numpy as np
import matplotlib
from matplotlib.cm import Greys, Greys_r
import matplotlib.patches as mpatches
import hops.pylightcurve3 as plc
from matplotlib.backend_bases import MouseEvent as mpl_MouseEvent
from astroquery.simbad import Simbad

from hops.application_windows import MainWindow

from scipy.optimize import curve_fit as scipy_curve_fit
import warnings

def curve_fit(*args, **kwargs):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message='Covariance of the parameters could not be estimated')
        return scipy_curve_fit(*args, **kwargs)


def test_coordinates2(ra_string, dec_string):

    try:
        coord = plc.Target(plc.Hours(ra_string), plc.Degrees(dec_string))
        return [True, 'Coordinates\naccepted']
    except:
        return [False, 'Wrong\ncoordinates']



def test_date(year_month_string):

    pass_test = True

    if len(year_month_string.split()) != 2:
        pass_test = False

    elif len(year_month_string.split()[0]) != 4:
        pass_test = False

    elif len(year_month_string.split()[1]) != 2:
        pass_test = False
    else:
        try:
            year = int(year_month_string.split()[0])
            month = int(year_month_string.split()[1])
            if int(month) < 1 or int(month) > 12:
                pass_test = False
            if int(year) < 0:
                pass_test = False
        except:
            pass_test = False

    if pass_test:
        return[True, 'Date\naccepted']
    else:
        return [False, 'Wrong\ndate']


class ObservingPlannerWindow(MainWindow):

    def __init__(self, log):

        MainWindow.__init__(self, log, name='HOPS - Observing Planner', position=2)

        # get variables from log and set as tk variables those to be modified

        self.target = self.Label(text=' ', instance=str)
        self.target_search = self.Entry(value='M42', instance=str)

        self.target_ra = self.Entry(value=' ', instance=str, command=self.update_window)
        self.target_dec = self.Entry(value=' ', instance=str, command=self.update_window)
        self.target_ra_dec_test = self.Label(text=' ')

        now = datetime.datetime.now()
        self.obs_year_month = self.Entry(value='{0} {1}'.format(now.year, str(now.month).zfill(2)), instance=str,
                                         command=self.update_window)
        self.obs_year_month_test = self.Label(text=' ')

        self.latitude = self.Entry(value=self.log.get_param('observatory_lat'), instance=str)
        self.longitude = self.Entry(value=self.log.get_param('observatory_long'), instance=str)
        self.horizon_s = self.Entry(value=self.log.get_param('observatory_horizon_s'), instance=float)
        self.horizon_sw = self.Entry(value=self.log.get_param('observatory_horizon_sw'), instance=float)
        self.horizon_w = self.Entry(value=self.log.get_param('observatory_horizon_w'), instance=float)
        self.horizon_nw = self.Entry(value=self.log.get_param('observatory_horizon_nw'), instance=float)
        self.horizon_n = self.Entry(value=self.log.get_param('observatory_horizon_n'), instance=float)
        self.horizon_ne = self.Entry(value=self.log.get_param('observatory_horizon_ne'), instance=float)
        self.horizon_e = self.Entry(value=self.log.get_param('observatory_horizon_e'), instance=float)
        self.horizon_se = self.Entry(value=self.log.get_param('observatory_horizon_se'), instance=float)
        self.observatory = self.Entry(value=self.log.get_param('observatory'), instance=str)
        self.timezone = self.Entry(value=self.log.get_param('observatory_time_zone'), instance=float)

        # create widgets

        self.figure = self.FigureWindow(figsize=(6, 6.5), show_nav=True)
        self.ax1 = self.figure.figure.add_subplot(111)
        self.figure.figure.subplots_adjust(left=0.1, right=1 - 0.05, bottom=0.1, top=0.85)

        self.plot_button = self.Button(text='RESET PLOT', command=self.plot)

        # setup window

        self.setup_window([
            [[self.figure, 4, 1, 30]],
            [[self.Label(text='Observatory'), 2]],
            [[self.observatory, 2]],
            [],
            [[self.Label(text='Latitude'), 1],
             [self.Label(text='Horizon altitude S (deg)'), 2],
             [self.horizon_s, 3]],
            [[self.latitude, 1],
             [self.Label(text='Horizon altitude SW (deg)'), 2],
             [self.horizon_sw, 3]],
            [[self.Label(text='Horizon altitude W (deg)'), 2],
             [self.horizon_w, 3]],
            [[self.Label(text='Longitude'), 1],
             [self.Label(text='Horizon altitude NW (deg)'), 2],
             [self.horizon_nw, 3]],
            [[self.longitude, 1],
             [self.Label(text='Horizon altitude N (deg)'), 2],
             [self.horizon_n, 3]],
            [[self.Label(text='Horizon altitude NE (deg)'), 2],
             [self.horizon_ne,3]],
            [[self.Label(text='Time Zone'), 1],
             [self.Label(text='Horizon altitude E (deg)'), 2],
             [self.horizon_e, 3]],
            [[self.timezone, 1],
             [self.Label(text='Horizon altitude SE (deg)'), 2],
             [self.horizon_se, 3]],
            [],
            [[self.Label(text='Object'), 1]],
            [[self.target_search, 1],
             [self.Button(text='SEARCH RA/DEC', command=self.search_object), 2], [self.target, 3]],
            [],
            [[self.Label(text='Manual target RA DEC\n(hh:mm:ss +/-dd:mm:ss)'), 1, 1, 2], [self.target_ra, 2],
             [self.target_ra_dec_test, 3, 1, 2]],
            [[self.target_dec, 2]],
            [[self.Label(text='Observation year and month\n(yyyy mm)'), 1], [self.obs_year_month, 2],
             [self.obs_year_month_test, 3]],
            [],
            [[self.plot_button, 2]],
            [],
        ])

        self.search_object()

    def update_window(self, *event):

        if not event:
            pass

        check_ra_dec = test_coordinates2(self.target_ra.get(), self.target_dec.get())
        self.target_ra_dec_test.set(check_ra_dec[1])
        check_year_month = test_date(self.obs_year_month.get())
        self.obs_year_month_test.set(check_year_month[1])

        if check_ra_dec[0] and check_year_month[0]:
            self.plot_button.activate()
        else:
            self.plot_button.disable()

    def search_object(self):

        try:
            result_table = Simbad.query_object(self.target_search.get())[0]

            if result_table:
                try:
                    xx = result_table['MAIN_ID'].decode("utf-8")
                except:
                    xx = result_table['MAIN_ID']
                self.target.set(xx)
                self.target_ra.set(str(result_table['RA']))
                self.target_dec.set(str(result_table['DEC']))
            else:
                self.target_ra.set('')
                self.target_dec.set('')
                self.target.set('No results')

        except:
            self.target_ra.set('')
            self.target_dec.set('')
            self.target.set('No results')

        self.plot()
        self.update_window()

    def plot(self):
        self.ax1.cla()

        try:

            name = self.target.get()
            if name == 'No results':
                name = self.target_search.get()

            self.avc_plot(self.latitude.get(), self.longitude.get(), self.timezone.get(),
                          [
                              [0, self.horizon_s.get()],
                              [45, self.horizon_sw.get()],
                              [90, self.horizon_w.get()],
                              [135, self.horizon_nw.get()],
                              [180, self.horizon_n.get()],
                              [225, self.horizon_ne.get()],
                              [270, self.horizon_e.get()],
                              [315, self.horizon_se.get()],

                          ],
                          self.target_ra.get(), self.target_dec.get(), self.obs_year_month.get(), self.ax1,
                          name, self.observatory.get())

        except:
            pass

        self.figure.draw()

    def target_azimuth_altitude(self, target_ra, target_dec, observatory_latitude, sidereal_time):

        observatory_latitude *= np.pi / 180
        target_dec *= np.pi / 180

        ha = sidereal_time - target_ra / 15
        if ha < 0:
            ha += 24

        altitude = np.arcsin(np.clip(np.sin(target_dec) * np.sin(observatory_latitude)
                                     + np.cos(target_dec) * np.cos(observatory_latitude) * np.cos(
            ha * 15 * np.pi / 180), -1, 1))

        azimuth = np.pi - np.arccos(np.clip((np.sin(target_dec) - np.sin(altitude) * np.sin(observatory_latitude)) /
                                            (np.cos(altitude) * np.cos(observatory_latitude)), -1, 1))

        if ha >= 12:
            azimuth = 2 * np.pi - azimuth

        return azimuth * 180 / np.pi, altitude * 180 / np.pi

    def get_target_events(self, target_ra, target_dec, observatory_latitude, horizon):

        # horizon

        horizon_list = []
        for horizon_line in horizon.split('\n'):
            if horizon_line != '':
                horizon_list.append(horizon_line.split())
        if len(horizon_list) == 1:
            def horizon(azimuth):
                return float(horizon_list[0][0])
        else:
            horizon_list.append([360.0, horizon_list[0][0]])
            horizon_data = np.swapaxes(np.array(horizon_list, dtype=float), 0, 1)
            horizon = plc.interp1d(horizon_data[0], horizon_data[1])

        # horizon

        sidereal_time_list = list(np.arange(0, 24, 0.1)) + [24]
        sidereal_time_altitude_list = []

        for sidereal_time in sidereal_time_list:
            azimuth, altitude = self.target_azimuth_altitude(target_ra, target_dec, observatory_latitude, sidereal_time)
            sidereal_time_altitude_list.append([sidereal_time, altitude - horizon(azimuth)])

        sidereal_time_altitude_list.sort()
        sidereal_time_altitude_list = np.swapaxes(sidereal_time_altitude_list, 0, 1)
        sidereal_time_altitude_function = plc.interp1d(np.array(sidereal_time_altitude_list[0]),
                                                       np.array(sidereal_time_altitude_list[1]))

        def target_horizon_diference(x, st):
            return sidereal_time_altitude_function(st)

        # find events

        events = []

        test_alt = sidereal_time_altitude_function(sidereal_time_list[0])
        for sidereal_time in sidereal_time_list:
            new_test_alt = sidereal_time_altitude_function(sidereal_time)
            if test_alt * new_test_alt < 0:
                popt, pcov = curve_fit(target_horizon_diference, [0], [0], p0=[sidereal_time - 0.1])
                if test_alt < new_test_alt:
                    events.append([popt[0], 'target_rise'])
                else:
                    events.append([popt[0], 'target_set'])
            else:
                pass

            test_alt = new_test_alt

        return events, sidereal_time_altitude_function

    def avc_plot(self, latitude, longitude, tmzn, horizon, target_ra, target_dec, year_mont_string, ax, name,
                 observatory_name):
        ax.cla()

        target = plc.Target(plc.Hours(target_ra), plc.Degrees(target_dec))
        observatory = plc.Observatory(plc.Degrees(latitude), plc.Degrees(longitude), tmzn, horizon)
        observation = plc.Observation(target, observatory)

        year = int(year_mont_string.split()[0])
        month = int(year_mont_string.split()[1])

        months = ['xx', 'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September',
                  'October',
                  'November', 'December']
        if (year - 2000) / 4.0 - int((year - 2000) / 4.0) == 0:
            days = [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        else:
            days = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        time_0 = plc.LT('{0}-{1}-1 12:00:00'.format(year, month), observatory)
        jd_0 = time_0.jd

        time_1 = plc.JD(jd_0 + days[month])

        events = []

        # mid-day splits

        for jj in range(days[month] + 1):
            events.append([plc.JD(time_0.jd + jj), 'mid-day'])

        # target rise/set events

        events += observation.rise_set_events(time_0, time_1)

        # sun rise/set events

        for jj in range(days[month]):

            check_time = plc.JD(time_0.jd + jj)
            check_st = check_time.lst(observatory).hours

            sun = check_time.get_sun()

            # sun rise/set

            if - (90 - observatory.latitude.deg_pm) < sun.dec.deg_pm < 90 - observatory.latitude.deg_pm:

                rise_ha = np.arccos(- sun.dec.tan * observatory.latitude.tan) * 12 / np.pi

                if rise_ha < 12:
                    set_ha = rise_ha
                    rise_ha = 24 - rise_ha
                else:
                    set_ha = 24 - rise_ha

                rise_st = rise_ha + sun.ra.hours
                if rise_st > 24:
                    rise_st -= 24

                set_st = set_ha + sun.ra.hours
                if set_st > 24:
                    set_st -= 24

                if rise_st < check_st:
                    next_rise_in_st_hours = 24 + rise_st - check_st
                else:
                    next_rise_in_st_hours = rise_st - check_st

                if set_st < check_st:
                    next_set_in_st_hours = 24 + set_st - check_st
                else:
                    next_set_in_st_hours = set_st - check_st

                dt = next_rise_in_st_hours * (23.9344696 / 24)
                if dt < 24:
                    events.append([plc.JD(check_time.jd + dt / 24), 'sun_rise'])

                dt = next_set_in_st_hours * (23.9344696 / 24)
                if dt < 24:
                    events.append([plc.JD(check_time.jd + dt / 24), 'sun_set'])

            # sun -18 rise/set

            if - (90 - observatory.latitude.deg_pm + 18.0) < sun.dec.deg_pm < 90 - (observatory.latitude.deg_pm + 18):

                rise_ha = np.arccos(np.sin((-18.0) * np.pi / 180) / sun.dec.cos / observatory.latitude.cos
                                    - sun.dec.tan * observatory.latitude.tan) * 12 / np.pi
                if rise_ha < 12:
                    set_ha = rise_ha
                    rise_ha = 24 - rise_ha
                else:
                    set_ha = 24 - rise_ha

                rise_st = rise_ha + sun.ra.hours
                if rise_st > 24:
                    rise_st -= 24

                set_st = set_ha + sun.ra.hours
                if set_st > 24:
                    set_st -= 24

                if rise_st < check_st:
                    next_rise_in_st_hours = 24 + rise_st - check_st
                else:
                    next_rise_in_st_hours = rise_st - check_st

                if set_st < check_st:
                    next_set_in_st_hours = 24 + set_st - check_st
                else:
                    next_set_in_st_hours = set_st - check_st

                dt = next_rise_in_st_hours * (23.9344696 / 24)
                if dt < 24:
                    events.append([plc.JD(check_time.jd + dt / 24), 'sun_rise_18'])

                dt = next_set_in_st_hours * (23.9344696 / 24)
                if dt < 24:
                    events.append([plc.JD(check_time.jd + dt / 24), 'sun_set_18'])

        events2 = [[ff[0].jd, ff[0], ff[1]] for ff in events]
        events2.sort(key=lambda ff: ff[0])

        #
        maxalt = str(round(observation.max_altitude.deg_pm, 1))

        ax.xaxis.tick_top()

        ax.set_title(observatory_name + '\n' + name + '   ' + months[month] + ' ' + str(year) +
                     '    max. alt. = ' + maxalt + ' degrees')
        ax.set_xlim((0, 1))
        ax.set_xlabel('HOUR (UTC{0:+.1f})'.format(tmzn))
        ax.set_xticks(np.arange(0, 24.5, 1))
        ax.set_xticklabels(('12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23',
                            '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'))
        ax.set_ylim((0, days[month] + 1))
        ax.set_ylabel('DAY')
        ax.set_yticks(np.arange(1, days[month] + 0.5, 1))
        ax.tick_params(bottom=True, top=True, left=True, right=True, labelbottom=True, labeltop=True,
                       labelright=False, labelleft=True)
        ax.grid(True, axis='y', linestyle='--')

        check_full_moon = plc.UTC('2000-1-21 04:41:00')

        for jj, ii in enumerate(events2[:-1]):

            moonphase = (np.sin((float(ii[0] + 0.5) - float(check_full_moon.jd)) * np.pi / 29.530589)) ** 2

            test_jd = 0.5 * (ii[0] + events2[jj + 1][0])
            dt_jd = 0.5 * (events2[jj + 1][0] - ii[0])

            day = 1 + int(test_jd - jd_0)

            time_range = [(ii[0] - jd_0 - int(ii[0] - jd_0)) * 24,
                          (events2[jj + 1][0] - jd_0 - int(events2[jj + 1][0] - jd_0)) * 24]
            if time_range[1] == 0:
                time_range[1] = 24

            alpha = 1
            if not observation.is_target_visible(plc.JD(ii[0] + dt_jd)):
                color = 'w'
                alpha = 0
            else:
                sun_az, sun_alt = observation.sun_azimuth_altitude(plc.JD(ii[0] + dt_jd))
                if sun_alt.deg_pm > 0:
                    color = 'y'
                elif sun_alt.deg_pm > -18:
                    color = 'r'
                else:
                    color = str(0.8 * (1 - moonphase))

            ax.plot(time_range, [day, day], linewidth=2.5, color=color, alpha=alpha)

            shift = {'left': +0.3, 'right': -0.3}

            if ii[2] == 'target_set':
                ax.plot(time_range[0], day, 'k*', mec='k', markersize=8)
                if time_range[0] > 20.5:
                    align = 'right'
                else:
                    align = 'left'
                ax.text(time_range[0] + shift[align], day + 0.4,
                        'set: ' + (ii[1].utc + datetime.timedelta(days=tmzn / 24)).isoformat().split('T')[1][:5],
                        va='center', ha=align, fontsize=9)

            if ii[2] == 'target_rise':
                ax.plot(time_range[0], day, 'w*', mec='k', markersize=8, markeredgewidth=0.5)
                if time_range[0] < 3.5:
                    align = 'left'
                else:
                    align = 'right'
                ax.text(time_range[0] + shift[align], day + 0.4,
                        'rise: ' + (ii[1].utc + datetime.timedelta(days=tmzn / 24)).isoformat().split('T')[1][:5],
                        va='center', ha=align, fontsize=9)
