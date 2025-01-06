
import numpy as np
import datetime
import exoclock

from hops.application_windows import MainWindow


def test_coordinates2(ra_string, dec_string):

    try:
        exoclock.FixedTarget(exoclock.Hours(ra_string), exoclock.Degrees(dec_string))
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
            target = exoclock.simbad_search_by_name(self.target_search.get())

            if target:
                self.target.set(target.name)
                self.target_ra.set(target.ra.hms())
                self.target_dec.set(target.dec.dms_coord())
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

    def avc_plot(self, latitude, longitude, tmzn, horizon, target_ra, target_dec, year_mont_string, ax, name,
                 observatory_name):
        ax.cla()

        target = exoclock.FixedTarget(exoclock.Hours(target_ra), exoclock.Degrees(target_dec))
        observatory = exoclock.Observatory(exoclock.Degrees(latitude), exoclock.Degrees(longitude), tmzn, horizon)

        year = int(year_mont_string.split()[0])
        month = int(year_mont_string.split()[1])

        months = ['xx', 'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September',
                  'October',
                  'November', 'December']
        if (year - 2000) / 4.0 - int((year - 2000) / 4.0) == 0:
            days = [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        else:
            days = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        jd_0 = exoclock.Moment('{0}-{1}-01 12:00:00'.format(year, month)).jd_utc() - tmzn / 24
        time_0 = exoclock.Moment(jd_utc=jd_0)

        events = []

        # mid-day splits

        for jj in range(days[month] + 1):
            events.append([exoclock.Moment(jd_utc=time_0.jd_utc() + jj), 'mid-day'])

        # target rise/set events

        events += observatory.target_horizon_crossings(target, time_0, days[month])

        # sun rise/set events

        for jj in range(days[month]):

            check_time = exoclock.Moment(jd_utc=time_0.jd_utc() + jj)

            sun = exoclock.Sun(check_time)
            fixed_sun = exoclock.FixedTarget(sun.ra, sun.dec)

            # # sun rise/set
            sun_events = observatory.target_horizon_crossings(fixed_sun, check_time, 1)
            for i in range(len(sun_events)):
                if sun_events[i][1] == 'rise':
                    sun_events[i][1] = 'sun_rise'
                elif sun_events[i][1] == 'set':
                    sun_events[i][1] = 'sun_set'

            events += sun_events

            # sun -18 rise/set

            twilight_events = observatory.target_altitude_crossings(fixed_sun, check_time, 1, exoclock.Degrees(-18))
            for i in range(len(twilight_events)):
                if twilight_events[i][1] == 'rise':
                    twilight_events[i][1] = 'sun_rise_18'
                elif twilight_events[i][1] == 'set':
                    twilight_events[i][1] = 'sun_set_18'

            events += twilight_events

        events2 = [[ff[0].jd_utc(), ff[0], ff[1]] for ff in events]
        events2.sort(key=lambda ff: ff[0])

        # for i in events2:
        #     print(i)

        ax.xaxis.tick_top()

        ax.set_title(observatory_name + '\n' + name + '   ' + months[month] + ' ' + str(year))
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

        for jj, ii in enumerate(events2[:-1]):

            test_jd = 0.5 * (ii[0] + events2[jj + 1][0])
            check_time = exoclock.Moment(jd_utc=test_jd)

            day = 1 + int(test_jd - jd_0)

            time_range = [(ii[0] - jd_0 - int(ii[0] - jd_0)) * 24,
                          (events2[jj + 1][0] - jd_0 - int(events2[jj + 1][0] - jd_0)) * 24]
            if time_range[1] == 0:
                time_range[1] = 24

            alpha = 1
            if not observatory.is_target_visible(target, check_time):
                color = 'w'
                alpha = 0
            else:
                sun_az, sun_alt = observatory.target_azimuth_altitude(exoclock.Sun(check_time), check_time)
                if sun_alt.deg_coord() > 0:
                    color = 'y'
                elif sun_alt.deg_coord() > -18:
                    color = 'r'
                else:
                    color = str(0.8 * exoclock.Moon(check_time).illumination())

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
