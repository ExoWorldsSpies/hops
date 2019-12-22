from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .hops_basics import *


def initialise_window(window, window_name, windows_to_hide, windows_to_close, exit_python):

    def exit_command():

        for i in windows_to_close:
            i.destroy()

        for i in windows_to_hide:
            i.withdraw()

        if exit_python:
            os._exit(-1)

    window.wm_title(window_name)
    window.protocol('WM_DELETE_WINDOW', exit_command)

    window.withdraw()


def setup_window(window, objects, main_font=None, button_font=None, entries_bd=3):

    if button_font is None:
        button_font = ['times', 15, 'bold']

    if main_font is None:
        main_font = ['times', 15]

    for row in range(len(objects)):
        if len(objects[row]) == 0:
            label_empty = Label(window, text='')
            label_empty.grid(row=row, column=100)
        else:
            for obj in objects[row]:

                if obj[0].winfo_class() == 'Button':
                    obj[0].configure(font=button_font)
                elif obj[0].winfo_class() == 'Entry':
                    obj[0].configure(bd=entries_bd, font=main_font)
                elif obj[0].winfo_class() in ['Label', 'Radiobutton']:
                    obj[0].configure(font=main_font)

                if len(obj) == 4:
                    obj[0].grid(row=row, column=obj[1], columnspan=obj[2], rowspan=obj[3])
                elif len(obj) == 3:
                    obj[0].grid(row=row, column=obj[1], columnspan=obj[2])
                else:
                    obj[0].grid(row=row, column=obj[1])


def finalise_window(window, position=5, topmost=False):

    window.update_idletasks()

    if position == 1:
        x = 0
        y = 0

    elif position == 2:
        x = (window.winfo_screenwidth() - window.winfo_reqwidth()) / 2
        y = 0

    elif position == 3:
        x = window.winfo_screenwidth() - window.winfo_reqwidth()
        y = 0

    elif position == 4:
        x = 0
        y = (window.winfo_screenheight() - window.winfo_reqheight()) / 2

    elif position == 5:
        x = (window.winfo_screenwidth() - window.winfo_reqwidth()) / 2
        y = (window.winfo_screenheight() - window.winfo_reqheight()) / 2

    elif position == 6:
        x = window.winfo_screenwidth() - window.winfo_reqwidth()
        y = (window.winfo_screenheight() - window.winfo_reqheight()) / 2

    elif position == 7:
        x = 0
        y = window.winfo_screenheight() - window.winfo_reqheight()

    elif position == 8:
        x = (window.winfo_screenwidth() - window.winfo_reqwidth()) / 2
        y = window.winfo_screenheight() - window.winfo_reqheight()

    elif position == 9:
        x = window.winfo_screenwidth() - window.winfo_reqwidth()
        y = window.winfo_screenheight() - window.winfo_reqheight()

    else:
        x = 0
        y = 0

    window.geometry('+%d+%d' % (x, y))

    window.update_idletasks()

    window.lift()
    window.wm_attributes("-topmost", 1)
    if not topmost:
        window.after_idle(window.attributes, '-topmost', 0)

    window.deiconify()


def test_float_positive_input(input_str, typing):

    if typing == '1':
        try:
            if float(input_str) >= 0:
                return True
            else:
                return False
        except ValueError:
            return False

    else:
        return True


def test_coordinates(ra_string, dec_string):

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


def target_azimuth_altitude(target_ra, target_dec, observatory_latitude, sidereal_time):

    observatory_latitude *= np.pi / 180
    target_dec *= np.pi / 180

    ha = sidereal_time - target_ra / 15
    if ha < 0:
        ha += 24

    altitude = np.arcsin(np.clip(np.sin(target_dec) * np.sin(observatory_latitude)
                         + np.cos(target_dec) * np.cos(observatory_latitude) * np.cos(ha * 15 * np.pi / 180), -1, 1))

    azimuth = np.pi - np.arccos(np.clip((np.sin(target_dec) - np.sin(altitude) * np.sin(observatory_latitude))/
                                (np.cos(altitude) * np.cos(observatory_latitude)), -1, 1))

    if ha >= 12:
        azimuth = 2 * np.pi - azimuth

    return azimuth * 180 / np.pi, altitude * 180 / np.pi


def get_target_events(target_ra, target_dec, observatory_latitude, horizon):

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
        azimuth, altitude = target_azimuth_altitude(target_ra, target_dec, observatory_latitude, sidereal_time)
        sidereal_time_altitude_list.append([sidereal_time, altitude - horizon(azimuth)])

    sidereal_time_altitude_list.sort()
    sidereal_time_altitude_list = np.swapaxes(sidereal_time_altitude_list, 0, 1)
    sidereal_time_altitude_function = plc.interp1d(np.array(sidereal_time_altitude_list[0]), np.array(sidereal_time_altitude_list[1]))

    def target_horizon_diference(x, st):
        return sidereal_time_altitude_function(st)

    # plt.plot(sidereal_time_altitude_list[0], sidereal_time_altitude_list[1], 'o')
    # plt.show()

    # find events

    events = []

    test_alt = sidereal_time_altitude_function(sidereal_time_list[0])
    for sidereal_time in sidereal_time_list:
        new_test_alt = sidereal_time_altitude_function(sidereal_time)
        if test_alt * new_test_alt < 0:
            popt, pcov = plc.curve_fit(target_horizon_diference, [0], [0], p0=[sidereal_time - 0.1])
            if test_alt < new_test_alt:
                events.append([popt[0], 'target_rise'])
            else:
                events.append([popt[0], 'target_set'])
        else:
            pass

        test_alt = new_test_alt

    return events, sidereal_time_altitude_function


def avc_plot(latitude, longitude, tmzn, horizon, target_ra, target_dec, year_mont_string, ax, name, observatory_name):
    ax.cla()

    target = plc.Target(plc.Hours(target_ra), plc.Degrees(target_dec))
    observatory = plc.Observatory(plc.Degrees(latitude), plc.Degrees(longitude), tmzn, horizon)
    observation = plc.Observation(target, observatory)

    year = int(year_mont_string.split()[0])
    month = int(year_mont_string.split()[1])

    months = ['xx', 'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October',
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
    ax.set_xlabel('HOUR / UTC {0:+.1f}'.format(tmzn))
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

        test_jd = 0.5*(ii[0] + events2[jj+1][0])
        dt_jd = 0.5*(events2[jj+1][0] - ii[0])

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
            elif sun_alt .deg_pm > -18:
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
            ax.text(time_range[0] + shift[align], day + 0.4, 'set: ' + (ii[1].utc + datetime.timedelta(days=tmzn/24)).isoformat().split('T')[1][:5],
                    va='center', ha=align, fontsize=9)

        if ii[2] == 'target_rise':
            ax.plot(time_range[0], day, 'w*', mec='k', markersize=8, markeredgewidth=0.5)
            if time_range[0] < 3.5:
                align = 'left'
            else:
                align = 'right'
            ax.text(time_range[0] + shift[align], day + 0.4, 'rise: ' + (ii[1].utc + datetime.timedelta(days=tmzn/24)).isoformat().split('T')[1][:5],
                    va='center', ha=align, fontsize=9)


def run_observing_planner():

    # #########
    # create and initialise the window
    # #########

    root = Tk()
    root2 = Tk()

    initialise_window(root, 'Observation planner', [], [root, root2], False)
    initialise_window(root2, 'Observation planner', [root2], [], False)

    # get variables from log and set as tk variables those to be modified

    object_search = StringVar(root, value='Orion nebula')
    skyobject = StringVar(root, value=' ')
    target_ra = StringVar(root, value=' ')
    target_dec = StringVar(root, value=' ')
    now = datetime.datetime.now()
    obs_year_month = StringVar(root, value='{0} {1}'.format(now.year, now.month))
    latitude = StringVar(root, value=read_local_log_profile('observatory_lat'))
    longitude = StringVar(root, value=read_local_log_profile('observatory_long'))
    horizon_s = StringVar(root, value=read_local_log_profile('observatory_horizon_s'))
    horizon_sw = StringVar(root, value=read_local_log_profile('observatory_horizon_sw'))
    horizon_w = StringVar(root, value=read_local_log_profile('observatory_horizon_w'))
    horizon_nw = StringVar(root, value=read_local_log_profile('observatory_horizon_nw'))
    horizon_n = StringVar(root, value=read_local_log_profile('observatory_horizon_n'))
    horizon_ne = StringVar(root, value=read_local_log_profile('observatory_horizon_ne'))
    horizon_e = StringVar(root, value=read_local_log_profile('observatory_horizon_e'))
    horizon_se = StringVar(root, value=read_local_log_profile('observatory_horizon_se'))
    observatory = StringVar(root, value=read_local_log_profile('observatory'))
    timezone = IntVar(root, value=read_local_log_profile('observatory_time_zone'))

    # set progress variables, useful for updating the window

    update_object = BooleanVar(root, value=True)
    update_object_list = BooleanVar(root, value=True)
    open_root2 = BooleanVar(root, value=False)

    # create the plot in the additional window

    figure = matplotlib.figure.Figure(figsize=(7, 7))
    figure.patch.set_facecolor('white')
    ax1 = figure.add_subplot(111)
    figure.subplots_adjust(top=0.75)
    canvas = FigureCanvasTkAgg(figure, root2)
    canvas.get_tk_widget().pack()
    NavigationToolbar2TkAgg(canvas, root2)

    # create widgets

    observatory_label = Label(root, text='Observatory')
    observatory_entry = Entry(root, textvariable=observatory, width=30)

    latitude_label = Label(root, text='Latitude')
    latitude_entry = Entry(root, textvariable=latitude)

    longitude_label = Label(root, text='Longitude')
    longitude_entry = Entry(root, textvariable=longitude)

    horizon_s_label = Label(root, text='Horizon altitude S (deg)')
    horizon_s_entry = Entry(root, textvariable=horizon_s)

    horizon_sw_label = Label(root, text='Horizon altitude SW (deg)')
    horizon_sw_entry = Entry(root, textvariable=horizon_sw)

    horizon_w_label = Label(root, text='Horizon altitude W (deg)')
    horizon_w_entry = Entry(root, textvariable=horizon_w)

    horizon_nw_label = Label(root, text='Horizon altitude NW (deg)')
    horizon_nw_entry = Entry(root, textvariable=horizon_nw)

    horizon_n_label = Label(root, text='Horizon altitude N (deg)')
    horizon_n_entry = Entry(root, textvariable=horizon_n)

    horizon_ne_label = Label(root, text='Horizon altitude NE (deg)')
    horizon_ne_entry = Entry(root, textvariable=horizon_ne)

    horizon_e_label = Label(root, text='Horizon altitude E (deg)')
    horizon_e_entry = Entry(root, textvariable=horizon_e)

    horizon_se_label = Label(root, text='Horizon altitude SE (deg)')
    horizon_se_entry = Entry(root, textvariable=horizon_se)

    timezone_label = Label(root, text='Time Zone')
    timezone_entry = Entry(root, textvariable=timezone)

    target_ra_dec_label = Label(root, text='Manual target RA DEC\n(hh:mm:ss +/-dd:mm:ss)')
    target_ra_entry = Entry(root, textvariable=target_ra, width=30)
    target_dec_entry = Entry(root, textvariable=target_dec, width=30)
    target_ra_dec_test = Label(root, text=' ')

    obs_year_month_label = Label(root, text='Observation year and month\n(yyyy mm)')
    obs_year_month_entry = Entry(root, textvariable=obs_year_month, width=30)
    obs_year_month_test = Label(root, text=' ')

    object_label = Label(root, text='     Object     ')
    object_search_entry = Entry(root, textvariable=object_search)
    combostyle = ttk.Style()
    try:
        combostyle.theme_create('combostyle', parent='alt',
                            settings={'TCombobox': {'configure':
                                                    {'selectbackground': 'white',
                                                     'fieldbackground': 'white',
                                                     'background': 'white'}}})
    except:
        pass

    combostyle.theme_use('combostyle')
    object_entry = ttk.Combobox(root, textvariable=skyobject, state='readonly')

    search_object_button = Button(root, text='SEARCH')

    plot_button = Button(root, text='RESET PLOT')

    exit_avc_button = Button(root, text='EXIT')

    # define the function that updates the window

    def update_window(*event):

        if not event:
            pass

        if update_object_list.get():

            try:
                result_table = Simbad.query_object(object_search.get(), wildcard=True)

                if result_table:
                    object_entry['values'] = tuple([str(ff['MAIN_ID'])[2:-1] for ff in result_table])
                else:
                    object_entry['values'] = tuple([])

                if len(object_entry['values']) == 1:
                    skyobject.set(object_entry['values'][0])
                    update_object.set(True)
                else:
                    skyobject.set('Choose Object')

            except requests.exceptions.ConnectionError:
                skyobject.set('No connection')

            update_object_list.set(False)

        if update_object.get():

            try:
                result_table = Simbad.query_object(skyobject.get())[0]
                target_ra.set(str(result_table['RA']))
                target_dec.set(str(result_table['DEC']))

            except requests.exceptions.ConnectionError:
                target_ra.set('Import manually')
                target_dec.set('Import manually')
                skyobject.set('No connection')

            update_object.set(False)

        object_entry.selection_clear()

        check_ra_dec = test_coordinates(target_ra_entry.get(), target_dec_entry.get(),)
        target_ra_dec_test.configure(text=check_ra_dec[1])
        check_year_month = test_date(obs_year_month_entry.get())
        obs_year_month_test.configure(text=check_year_month[1])

        if check_ra_dec[0] and check_year_month[0]:
            plot_window = True
        else:
            plot_window = False

        if plot_window:

            # try:
                avc_plot(latitude.get(), longitude.get(), timezone.get(),
                         [
                             [0, horizon_s.get()],
                             [45, horizon_sw.get()],
                             [90, horizon_w.get()],
                             [135, horizon_nw.get()],
                             [180, horizon_n.get()],
                             [225, horizon_ne.get()],
                             [270, horizon_e.get()],
                             [315, horizon_se.get()],

                         ],
                         target_ra.get(), target_dec.get(), obs_year_month.get(), ax1, skyobject.get(),
                         observatory.get())
            #
            # except IndexError:
            #     pass

        canvas.draw()

    update_window()

    # define actions for the different buttons, including calls to the function that updates the window

    def choose_object(entry):

        if not entry:
            return 0

        update_object.set(True)
        update_window()

    def search_object():

        update_object_list.set(True)
        update_window()

    def plot():

        open_root2.set(True)
        update_window()

        root2.deiconify()

    def exit_acv():
        root.destroy()
        root2.destroy()

    # connect actions to widgets

    object_entry.bind('<<ComboboxSelected>>', choose_object)
    search_object_button['command'] = search_object
    plot_button['command'] = plot
    exit_avc_button['command'] = exit_acv

    # setup window

    setup_window(root, [
        [],
        [[observatory_label, 2]],
        [[observatory_entry, 2]],
        [],
        [[latitude_label, 1], [horizon_s_label, 2], [horizon_s_entry, 3]],
        [[latitude_entry, 1], [horizon_sw_label, 2], [horizon_sw_entry, 3]],
        [[horizon_w_label, 2], [horizon_w_entry, 3]],
        [[longitude_label, 1],[horizon_nw_label, 2], [horizon_nw_entry, 3]],
        [[longitude_entry, 1], [horizon_n_label, 2], [horizon_n_entry, 3]],
        [[horizon_ne_label, 2], [horizon_ne_entry, 3]],
        [[timezone_label, 1], [horizon_e_label, 2], [horizon_e_entry, 3]],
        [[timezone_entry, 1], [horizon_se_label, 2], [horizon_se_entry, 3]],
        [],
        [[object_label, 1]],
        [[object_search_entry, 1], [object_entry, 2]],
        [[search_object_button, 1]],
        [],
        [[target_ra_dec_label, 1], [target_ra_entry, 2], [target_ra_dec_test, 3]],
        [[target_dec_entry, 2]],
        [[obs_year_month_label, 1], [obs_year_month_entry, 2], [obs_year_month_test, 3]],
        [],
        [[plot_button, 2]],
        [],
        [[exit_avc_button, 2]],
        [],
    ])

    # finalise and show  window

    finalise_window(root, 1)
    finalise_window(root2, 3)
    root.mainloop()






