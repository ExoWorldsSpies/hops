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


def test_coordinates(ra_dec_string):

    try:
        coord = SkyCoord(ra_dec_string, unit=(u.hourangle, u.deg))
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


def avc_plot(lat, long, tmzn, horizon, elevation, target_ra_dec_string, year_mont_string, ax, name, observatory_name):

    ax.cla()

    coord = SkyCoord(target_ra_dec_string, unit=(u.hourangle, u.deg))
    targetra = coord.ra.deg
    targetdec = coord.dec.deg
    year = int(year_mont_string.split()[0])
    month = int(year_mont_string.split()[1])
    horizon = float(horizon)

    months = ['xx', 'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October',
              'November', 'December']
    if (year - 2000) / 4.0 - int((year - 2000) / 4.0) == 0:
        days = [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    else:
        days = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    bear_mountain = EarthLocation(lat=lat, lon=long)
    utcoffset = tmzn * u.hour

    time_0 = Time('{0}-{1}-1 12:00:00'.format(year, month), location=bear_mountain) - utcoffset
    jd_0 = time_0.jd
    st_0 = time_0.sidereal_time('apparent').hour

    events = []

    # mid-day splits

    for jj in range(days[month] + 1):
        events.append([time_0 + jj * u.day, 'mid-day'])

    # target rise/set events

    if targetdec > - (90 - bear_mountain.lat.deg - horizon) and targetdec < 90 - (bear_mountain.lat.deg - horizon):

        alt_rise = horizon
        rise_ha = np.arccos(np.sin(alt_rise*np.pi/180)/np.cos(targetdec*np.pi/180)/np.cos(bear_mountain.lat.deg*np.pi/180)
                            - np.tan(targetdec*np.pi/180)*np.tan(bear_mountain.lat.deg*np.pi/180))*12/np.pi
        if rise_ha < 12:
            rise_ha = 24 - rise_ha

        rise_st = rise_ha + targetra / 15
        if rise_st > 24:
            rise_st -= 24

        alt_set = horizon
        set_ha = np.arccos(np.sin(alt_set*np.pi/180)/np.cos(targetdec*np.pi/180)/np.cos(bear_mountain.lat.deg*np.pi/180)
                           - np.tan(targetdec*np.pi/180)*np.tan(bear_mountain.lat.deg*np.pi/180))*12/np.pi
        if set_ha > 12:
            set_ha = 24 - set_ha

        set_st = set_ha + targetra / 15
        if set_st > 24:
            set_st -= 24

        if rise_st < st_0:
            next_rise_in_st_hours = 24 + rise_st - st_0
        else:
            next_rise_in_st_hours = rise_st - st_0

        if set_st < st_0:
            next_set_in_st_hours = 24 + set_st - st_0
        else:
            next_set_in_st_hours = set_st - st_0

        for jj in range(days[month] + 1):
            dt = (jj*24 + next_rise_in_st_hours) * (365.25/366.25)
            if dt < days[month] * 24:
                xx = Time('{0}-{1}-1 12:00:00'.format(year, month, jj), location=bear_mountain) - utcoffset + dt * u.hour
                events.append([xx, 'target_rise'])

        for jj in range(days[month] + 1):
            dt = (jj*24 + next_set_in_st_hours) * (365.25/366.25)
            if dt < days[month] * 24:
                xx = Time('{0}-{1}-1 12:00:00'.format(year, month, jj), location=bear_mountain) - utcoffset + dt * u.hour
                events.append([xx, 'target_set'])

    # sun rise/set events

    def get_sun_altitude(time_object):

        check_jd = int(time_object.jd)
        sun_t = plc.hjd_dict[check_jd]['t']
        sun_ra = plc.hjd_dict[check_jd]['ra']
        sun_dec = plc.hjd_dict[check_jd]['dec']

        sun_ra = plc.interp1d(sun_t, sun_ra, kind='cubic')(check_jd)
        sun_ra = sun_ra - int(sun_ra / (2 * np.pi)) * 2 * np.pi
        sun_ra *= 180 / np.pi
        sun_dec = plc.interp1d(sun_t, sun_dec, kind='cubic')(check_jd)

        st = time_object.sidereal_time('apparent').hour
        ha = st - sun_ra / 15
        if ha < 0:
            ha += 24

        alt = np.arcsin(np.sin(sun_dec)*np.sin(bear_mountain.lat.deg*np.pi/180)
                        + np.cos(sun_dec)*np.cos(bear_mountain.lat.deg*np.pi/180)*np.cos(ha*15*np.pi/180))

        return alt * 180 / np.pi

    def get_target_altitude(time_object):

        st = time_object.sidereal_time('apparent').hour
        ha = st - targetra / 15
        if ha < 0:
            ha += 24

        alt = np.arcsin(np.sin(targetdec*np.pi/180)*np.sin(bear_mountain.lat.deg*np.pi/180)
                        + np.cos(targetdec*np.pi/180)*np.cos(bear_mountain.lat.deg*np.pi/180)*np.cos(ha*15*np.pi/180))

        return alt * 180 / np.pi

    for jj in range(days[month]):

        check_time = Time('{0}-{1}-{2} 12:00:00'.format(year, month, int(jj+1)), location=bear_mountain) - utcoffset
        check_st = check_time.sidereal_time('apparent').hour

        check_jd = int((check_time + 12 * u.hour).jd)
        sun_t = plc.hjd_dict[check_jd]['t']
        sun_ra = plc.hjd_dict[check_jd]['ra']
        sun_dec = plc.hjd_dict[check_jd]['dec']

        sun_ra = plc.interp1d(sun_t, sun_ra, kind='cubic')(check_jd)
        sun_ra = sun_ra - int(sun_ra / (2 * np.pi)) * 2 * np.pi
        sun_ra *= 180 / np.pi
        sun_dec = plc.interp1d(sun_t, sun_dec, kind='cubic')(check_jd)

        # sun rise/set

        if sun_dec * 100 / np.pi > - (90 - bear_mountain.lat.deg) and sun_dec * 100 / np.pi < 90 - (bear_mountain.lat.deg):

            rise_ha = np.arccos(- np.tan(sun_dec)*np.tan(bear_mountain.lat.deg*np.pi/180))*12/np.pi

            if rise_ha < 12:
                set_ha = rise_ha
                rise_ha = 24 - rise_ha
            else:
                set_ha = 24 - rise_ha

            rise_st = rise_ha + sun_ra / 15
            if rise_st > 24:
                rise_st -= 24

            set_st = set_ha + sun_ra / 15
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

            dt = next_rise_in_st_hours * (365.25/366.25)
            if dt < 24:
                events.append([check_time + dt * u.hour, 'sun_rise'])

            dt = next_set_in_st_hours * (365.25 / 366.25)
            if dt < 24:
                events.append([check_time + dt * u.hour, 'sun_set'])

        # sun -18 rise/set

        if sun_dec * 100 / np.pi > - (90 - bear_mountain.lat.deg + 18) and sun_dec * 100 / np.pi < 90 - (bear_mountain.lat.deg + 18):

            rise_ha = np.arccos(np.sin((-18.0)*np.pi/180)/np.cos(sun_dec)/np.cos(bear_mountain.lat.deg*np.pi/180)
                                - np.tan(sun_dec)*np.tan(bear_mountain.lat.deg*np.pi/180))*12/np.pi
            if rise_ha < 12:
                set_ha = rise_ha
                rise_ha = 24 - rise_ha
            else:
                set_ha = 24 - rise_ha

            rise_st = rise_ha + sun_ra / 15
            if rise_st > 24:
                rise_st -= 24

            set_st = set_ha + sun_ra / 15
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

            dt = next_rise_in_st_hours * (365.25/366.25)
            if dt < 24:
                events.append([check_time + dt * u.hour, 'sun_rise 18'])

            dt = next_set_in_st_hours * (365.25 / 366.25)
            if dt < 24:
                events.append([check_time + dt * u.hour, 'sun_set 18'])

    events2 = [[ff[0].jd, ff[0], ff[1]] for ff in events]
    events2.sort()

    #
    maxalt = int(round(90 - bear_mountain.lat.deg + targetdec, 0))
    if maxalt > 90:
        maxalt = 180 - maxalt
    maxalt = str(maxalt)

    ax.xaxis.tick_top()

    ax.set_title(observatory_name + '\n' + name + '   ' + months[month] + ' ' + str(year) + '    max. alt. = ' + maxalt + ' degrees')
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

    check_full_moon = Time('1992-5-16 16:02:33')

    for jj, ii in enumerate(events2[:-1]):

        moonphase = (np.sin((float(ii[0] + 0.5) - float(check_full_moon.jd)) * np.pi / 29.530589)) ** 2

        test_jd = 0.5*(ii[0] + events2[jj+1][0])
        dt_jd = 0.5*(events2[jj+1][0] - ii[0])

        day = 1 + int(test_jd - jd_0)

        time_range = [(ii[0] - jd_0 - int(ii[0] - jd_0)) * 24,
                      (events2[jj + 1][0] - jd_0 - int(events2[jj + 1][0] - jd_0)) * 24]
        if time_range[1] == 0:
            time_range[1] = 24

        test_target = get_target_altitude(ii[1] + dt_jd*u.day)

        alpha=1
        if test_target < horizon:
            color = 'w'
            alpha = 0
        else:
            sun_alt = get_sun_altitude(ii[1] + dt_jd*u.day)
            if sun_alt > 0:
                color = 'y'
            elif sun_alt > -18:
                color = 'r'
            else:
                color = str(0.8 * (1 - moonphase))

        # print(day, time_range, color)

        ax.plot(time_range, [day, day], linewidth=2.5, color=color, alpha=alpha)

        shift = {'left': +0.3, 'right': -0.3}

        if ii[2] == 'target_set':
            ax.plot(time_range[0], day, 'k*', mec='k', markersize=8)
            if time_range[0] > 20.5 / 24:
                align = 'right'
            else:
                align = 'left'
            ax.text(time_range[0] + shift[align], day + 0.4, 'set: ' + str(ii[1] + tmzn*u.hour).split()[1][:5], va='center', ha=align, fontsize=9)

        if ii[2] == 'target_rise':
            ax.plot(time_range[0], day, 'w*', mec='k', markersize=8, markeredgewidth=0.5)
            if time_range[0] < 3.5 / 24:
                align = 'left'
            else:
                align = 'right'
            ax.text(time_range[0] + shift[align], day + 0.4, 'rise: ' + str(ii[1] + tmzn*u.hour).split()[1][:5], va='center', ha=align, fontsize=9)


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
    target_ra_dec = StringVar(root, value=' ')
    now = datetime.datetime.now()
    obs_year_month = StringVar(root, value='{0} {1}'.format(now.year, now.month))
    latitude = StringVar(root, value=read_local_log_profile('observatory_lat'))
    longitude = StringVar(root, value=read_local_log_profile('observatory_long'))
    horizon = StringVar(root, value='20.0')
    observatory = StringVar(root, value=read_local_log_profile('observatory'))
    elevation = IntVar(root, value=read_local_log_profile('observatory_elev'))
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

    horizon_label = Label(root, text='Horizon (deg)')
    horizon_entry = Entry(root, textvariable=horizon)

    timezone_label = Label(root, text='Time Zone')
    timezone_entry = Entry(root, textvariable=timezone)

    elevation_label = Label(root, text='Elevation (m)')
    elevation_entry = Entry(root, textvariable=elevation)

    target_ra_dec_label = Label(root, text='Manual target RA DEC\n(hh:mm:ss +/-dd:mm:ss)')
    target_ra_dec_entry = Entry(root, textvariable=target_ra_dec, width=30)
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
                target_ra_dec.set('{0} {1}'.format(result_table['RA'], result_table['DEC']))

            except requests.exceptions.ConnectionError:
                target_ra_dec.set('Import manually')
                skyobject.set('No connection')

            update_object.set(False)

        object_entry.selection_clear()

        check_ra_dec = test_coordinates(target_ra_dec_entry.get())
        target_ra_dec_test.configure(text=check_ra_dec[1])
        check_year_month = test_date(obs_year_month_entry.get())
        obs_year_month_test.configure(text=check_year_month[1])

        if check_ra_dec[0] and check_year_month[0]:
            plot_window = True
        else:
            plot_window = False

        if plot_window:

            try:
                avc_plot(latitude.get(), longitude.get(), timezone.get(), horizon.get(), elevation.get(),
                         target_ra_dec.get(), obs_year_month.get(), ax1, skyobject.get(), observatory.get())

            except IndexError:
                pass

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
        [[latitude_label, 1], [elevation_label, 2]],
        [[latitude_entry, 1], [elevation_entry, 2], [timezone_label, 3]],
        [[longitude_label, 1], [horizon_label, 2], [timezone_entry, 3]],
        [[longitude_entry, 1], [horizon_entry, 2]],
        [],
        [[object_label, 2]],
        [[object_search_entry, 2]],
        [[search_object_button, 2]],
        [],
        [[object_entry, 2]],
        [],
        [[target_ra_dec_label, 1], [target_ra_dec_entry, 2], [target_ra_dec_test, 3]],
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






