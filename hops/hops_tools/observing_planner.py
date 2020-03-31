

import sys


from tkinter import Tk, TclError
from tkinter import Label, Button, Entry, Checkbutton, Scrollbar, Listbox, PhotoImage, Radiobutton, Scale, Frame
from tkinter import StringVar, BooleanVar, DoubleVar, IntVar
from tkinter import DISABLED, NORMAL, END, RIGHT, LEFT, BOTH, Y, HORIZONTAL
from tkinter.ttk import Combobox, Style, Progressbar

import tkinter.ttk as ttk
import tkinter.filedialog as tkFileDialog
from tkinter.messagebox import *
from urllib.request import urlopen

import warnings
warnings.filterwarnings(
    'ignore', message='Matplotlib is building the font cache using fc-list. This may take a moment.')
warnings.filterwarnings(
    'ignore', message='The installed version of numexpr 2.4.4 is not supported in pandas and will be not be used')

import matplotlib
matplotlib.use('TkAgg')

import datetime
import os
import sys
import glob
import time
import yaml
import numpy as np
import shutil
import hops.pylightcurve3 as plc
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.patches as mpatch

from astropy.io import fits as pf
from scipy.optimize import curve_fit
from matplotlib.figure import Figure
from matplotlib.offsetbox import AnchoredText
from matplotlib.backend_bases import key_press_handler, MouseEvent
import matplotlib.gridspec as gridspec
try:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    NavigationToolbar2TkAgg = NavigationToolbar2Tk
except ImportError:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import FigureCanvasBase
import matplotlib.image as mpimg

from astroquery.simbad import Simbad
import requests
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation

import tkinter.scrolledtext as scrolledtext

import webbrowser

from hops.hops_tools.windows import MainWindow, AddOnWindow, SideWindow, setup_window

from hops.hops_tools.logs import log
from hops.hops_tools.tests import *
#
#
# def initialise_window(window, window_name, windows_to_hide, windows_to_close, exit_python):
#
#     def exit_command():
#
#         for i in windows_to_close:
#             i.destroy()
#
#         for i in windows_to_hide:
#             i.withdraw()
#
#         if exit_python:
#             os._exit(-1)
#
#     window.wm_title(window_name)
#     window.protocol('WM_DELETE_WINDOW', exit_command)
#
#     window.withdraw()
#
#
# def setup_window(window, objects, main_font=None, button_font=None, entries_bd=3):
#
#     if button_font is None:
#         button_font = ['times', 15, 'bold']
#
#     if main_font is None:
#         main_font = ['times', 15]
#
#     for row in range(len(objects)):
#         if len(objects[row]) == 0:
#             label_empty = Label(window, text='')
#             label_empty.grid(row=row, column=100)
#         else:
#             for obj in objects[row]:
#
#                 if obj[0].winfo_class() == 'Button':
#                     obj[0].configure(font=button_font)
#                 elif obj[0].winfo_class() == 'Entry':
#                     obj[0].configure(bd=entries_bd, font=main_font)
#                 elif obj[0].winfo_class() in ['Label', 'Radiobutton']:
#                     obj[0].configure(font=main_font)
#
#                 if len(obj) == 4:
#                     obj[0].grid(row=row, column=obj[1], columnspan=obj[2], rowspan=obj[3])
#                 elif len(obj) == 3:
#                     obj[0].grid(row=row, column=obj[1], columnspan=obj[2])
#                 else:
#                     obj[0].grid(row=row, column=obj[1])
#
#
# def finalise_window(window, position=5, topmost=False):
#
#     window.update_idletasks()
#
#     if position == 1:
#         x = 0
#         y = 0
#
#     elif position == 2:
#         x = (window.winfo_screenwidth() - window.winfo_reqwidth()) / 2
#         y = 0
#
#     elif position == 3:
#         x = window.winfo_screenwidth() - window.winfo_reqwidth()
#         y = 0
#
#     elif position == 4:
#         x = 0
#         y = (window.winfo_screenheight() - window.winfo_reqheight()) / 2
#
#     elif position == 5:
#         x = (window.winfo_screenwidth() - window.winfo_reqwidth()) / 2
#         y = (window.winfo_screenheight() - window.winfo_reqheight()) / 2
#
#     elif position == 6:
#         x = window.winfo_screenwidth() - window.winfo_reqwidth()
#         y = (window.winfo_screenheight() - window.winfo_reqheight()) / 2
#
#     elif position == 7:
#         x = 0
#         y = window.winfo_screenheight() - window.winfo_reqheight()
#
#     elif position == 8:
#         x = (window.winfo_screenwidth() - window.winfo_reqwidth()) / 2
#         y = window.winfo_screenheight() - window.winfo_reqheight()
#
#     elif position == 9:
#         x = window.winfo_screenwidth() - window.winfo_reqwidth()
#         y = window.winfo_screenheight() - window.winfo_reqheight()
#
#     else:
#         x = 0
#         y = 0
#
#     window.geometry('+%d+%d' % (x, y))
#
#     window.update_idletasks()
#
#     window.lift()
#     window.wm_attributes("-topmost", 1)
#     if not topmost:
#         window.after_idle(window.attributes, '-topmost', 0)
#
#     window.deiconify()
#
#
# def test_float_positive_input(input_str, typing):
#
#     if typing == '1':
#         try:
#             if float(input_str) >= 0:
#                 return True
#             else:
#                 return False
#         except ValueError:
#             return False
#
#     else:
#         return True
#
#

def test_coordinates2(ra_string, dec_string):

    try:
        coord = plc.Target(plc.Hours(ra_string), plc.Degrees(dec_string))
        return [True, 'Coordinates\naccepted']
    except:
        return [False, 'Wrong\ncoordinates']




class ObservingPlanner(AddOnWindow):
    
    def __init__(self):
    
        AddOnWindow.__init__(self, 'Observing Planner')
    
        self.frame1 = self.Frame()
        self.frame1.pack(side=LEFT)

    
        # get variables from log and set as tk variables those to be modified
    
        self.object_search = StringVar(self.frame1, value='Orion nebula')
        self.skyobject = StringVar(self.frame1, value=' ')
        self.target_ra = StringVar(self.frame1, value=' ')
        self.target_dec = StringVar(self.frame1, value=' ')
        now = datetime.datetime.now()
        self.obs_year_month = StringVar(self.frame1, value='{0} {1}'.format(now.year, str(now.month).zfill(2)))
        self.latitude = StringVar(self.frame1, value=log.read_local_log_profile('observatory_lat'))
        self.longitude = StringVar(self.frame1, value=log.read_local_log_profile('observatory_long'))
        self.horizon_s = StringVar(self.frame1, value=log.read_local_log_profile('observatory_horizon_s'))
        self.horizon_sw = StringVar(self.frame1, value=log.read_local_log_profile('observatory_horizon_sw'))
        self.horizon_w = StringVar(self.frame1, value=log.read_local_log_profile('observatory_horizon_w'))
        self.horizon_nw = StringVar(self.frame1, value=log.read_local_log_profile('observatory_horizon_nw'))
        self.horizon_n = StringVar(self.frame1, value=log.read_local_log_profile('observatory_horizon_n'))
        self.horizon_ne = StringVar(self.frame1, value=log.read_local_log_profile('observatory_horizon_ne'))
        self.horizon_e = StringVar(self.frame1, value=log.read_local_log_profile('observatory_horizon_e'))
        self.horizon_se = StringVar(self.frame1, value=log.read_local_log_profile('observatory_horizon_se'))
        self.observatory = StringVar(self.frame1, value=log.read_local_log_profile('observatory'))
        self.timezone = IntVar(self.frame1, value=log.read_local_log_profile('observatory_time_zone'))
    
        # set progress variables, useful for updating the window
    
        self.update_object = BooleanVar(self.frame1, value=True)
        self.update_object_list = BooleanVar(self.frame1, value=True)
    
        # create the plot in the additional window
    
        # figure = matplotlib.figure.Figure()
        figure = matplotlib.figure.Figure(figsize=(6, 7))
        figure.patch.set_facecolor('white')
        self.ax1 = figure.add_subplot(111)
        figure.subplots_adjust(left=0.1, right=1-0.05, bottom=0.1, top=0.85)
        self.canvas = self.FigureCanvasTkAgg(figure)
        self.canvas.get_tk_widget().pack(side=RIGHT)
        # self.NavigationToolbar2Tk(self.canvas)
    
        # create widgets
    
        self.observatory_label = Label(self.frame1, text='Observatory')
        self.observatory_entry = Entry(self.frame1, textvariable=self.observatory, width=30)
    
        self.latitude_label = Label(self.frame1, text='Latitude')
        self.latitude_entry = Entry(self.frame1, textvariable=self.latitude)
    
        self.longitude_label = Label(self.frame1, text='Longitude')
        self.longitude_entry = Entry(self.frame1, textvariable=self.longitude)
    
        self.horizon_s_label = Label(self.frame1, text='Horizon altitude S (deg)')
        self.horizon_s_entry = Entry(self.frame1, textvariable=self.horizon_s, width=10)
    
        self.horizon_sw_label = Label(self.frame1, text='Horizon altitude SW (deg)')
        self.horizon_sw_entry = Entry(self.frame1, textvariable=self.horizon_sw, width=10)
    
        self.horizon_w_label = Label(self.frame1, text='Horizon altitude W (deg)')
        self.horizon_w_entry = Entry(self.frame1, textvariable=self.horizon_w, width=10)
    
        self.horizon_nw_label = Label(self.frame1, text='Horizon altitude NW (deg)')
        self.horizon_nw_entry = Entry(self.frame1, textvariable=self.horizon_nw, width=10)
    
        self.horizon_n_label = Label(self.frame1, text='Horizon altitude N (deg)')
        self.horizon_n_entry = Entry(self.frame1, textvariable=self.horizon_n, width=10)
    
        self.horizon_ne_label = Label(self.frame1, text='Horizon altitude NE (deg)')
        self.horizon_ne_entry = Entry(self.frame1, textvariable=self.horizon_ne, width=10)
    
        self.horizon_e_label = Label(self.frame1, text='Horizon altitude E (deg)')
        self.horizon_e_entry = Entry(self.frame1, textvariable=self.horizon_e, width=10)
    
        self.horizon_se_label = Label(self.frame1, text='Horizon altitude SE (deg)')
        self.horizon_se_entry = Entry(self.frame1, textvariable=self.horizon_se, width=10)
    
        self.timezone_label = Label(self.frame1, text='Time Zone')
        self.timezone_entry = Entry(self.frame1, textvariable=self.timezone)
    
        self.target_ra_dec_label = Label(self.frame1, text='Manual target RA DEC\n(hh:mm:ss +/-dd:mm:ss)')
        self.target_ra_entry = Entry(self.frame1, textvariable=self.target_ra, width=30)
        self.target_ra_entry.bind(sequence='<KeyRelease>', func=self.update_window)
        self.target_dec_entry = Entry(self.frame1, textvariable=self.target_dec, width=30)
        self.target_dec_entry.bind(sequence='<KeyRelease>', func=self.update_window)
        self.target_ra_dec_test = Label(self.frame1, text=' ')
    
        self.obs_year_month_label = Label(self.frame1, text='Observation year and month\n(yyyy mm)')
        self.obs_year_month_entry = Entry(self.frame1, textvariable=self.obs_year_month, width=30)
        self.obs_year_month_test = Label(self.frame1, text=' ')
    
        self.object_label = Label(self.frame1, text='     Object     ')
        self.object_search_entry = Entry(self.frame1, textvariable=self.object_search)
        combostyle = self.Style()
        try:
            combostyle.theme_create('combostyle', parent='alt',
                                settings={'TCombobox': {'configure':
                                                        {'selectbackground': 'white',
                                                         'fieldbackground': 'white',
                                                         'background': 'white'}}})
        except:
            pass
    
        combostyle.theme_use('combostyle')
        self.object_entry = Combobox(self.frame1, textvariable=self.skyobject, state='readonly')
    
        self.search_object_button = Button(self.frame1, text='SEARCH')
    
        self.plot_button = Button(self.frame1, text='RESET PLOT')
    
        self.exit_avc_button = Button(self.frame1, text='EXIT')

        self.object_entry.bind('<<ComboboxSelected>>', self.choose_object)
        self.search_object_button['command'] = self.search_object
        self.plot_button['command'] = self.plot

        # setup window

        setup_window(self.frame1, [
            [],
            [[self.observatory_label, 2]],
            [[self.observatory_entry, 2]],
            [],
            [[self.latitude_label, 1], [self.horizon_s_label, 2], [self.horizon_s_entry, 3]],
            [[self.latitude_entry, 1], [self.horizon_sw_label, 2], [self.horizon_sw_entry, 3]],
            [[self.horizon_w_label, 2], [self.horizon_w_entry, 3]],
            [[self.longitude_label, 1], [self.horizon_nw_label, 2], [self.horizon_nw_entry, 3]],
            [[self.longitude_entry, 1], [self.horizon_n_label, 2], [self.horizon_n_entry, 3]],
            [[self.horizon_ne_label, 2], [self.horizon_ne_entry, 3]],
            [[self.timezone_label, 1], [self.horizon_e_label, 2], [self.horizon_e_entry, 3]],
            [[self.timezone_entry, 1], [self.horizon_se_label, 2], [self.horizon_se_entry, 3]],
            [],
            [[self.object_label, 1]],
            [[self.object_search_entry, 1], [self.object_entry, 2]],
            [[self.search_object_button, 1]],
            [],
            [[self.target_ra_dec_label, 1], [self.target_ra_entry, 2], [self.target_ra_dec_test, 3]],
            [[self.target_dec_entry, 2]],
            [[self.obs_year_month_label, 1], [self.obs_year_month_entry, 2], [self.obs_year_month_test, 3]],
            [],
            [[self.plot_button, 2]],
            [],
        ])

        self.update_window()
        self.plot()

    # define the function that updates the window

    def update_window(self, *event):

        if not event:
            pass

        if self.update_object_list.get():

            try:
                result_table = Simbad.query_object(self.object_search.get(), wildcard=True)

                if result_table:
                    self.object_entry['values'] = tuple([str(ff['MAIN_ID'])[2:-1] for ff in result_table])
                else:
                    self.object_entry['values'] = tuple([])

                if len(self.object_entry['values']) == 1:
                    self.skyobject.set(self.object_entry['values'][0])
                    self.update_object.set(True)
                else:
                    self.skyobject.set('Choose Object')

            except requests.exceptions.ConnectionError:
                self.skyobject.set('No connection')

            self.update_object_list.set(False)

        if self.update_object.get():

            try:
                result_table = Simbad.query_object(self.skyobject.get())[0]
                self.target_ra.set(str(result_table['RA']))
                self.target_dec.set(str(result_table['DEC']))

            except requests.exceptions.ConnectionError:
                self.target_ra.set('Import manually')
                self.target_dec.set('Import manually')
                self.skyobject.set('No connection')

            self.update_object.set(False)

        self.object_entry.selection_clear()

        check_ra_dec = test_coordinates2(self.target_ra_entry.get(), self.target_dec_entry.get())
        self.target_ra_dec_test.configure(text=check_ra_dec[1])
        check_year_month = test_date(self.obs_year_month_entry.get())
        self.obs_year_month_test.configure(text=check_year_month[1])

        if check_ra_dec[0] and check_year_month[0]:
            self.plot_button['state'] = NORMAL
        else:
            self.plot_button['state'] = DISABLED

    def choose_object(self, entry):

        if not entry:
            return 0

        self.update_object.set(True)
        self.update_window()

    def search_object(self):

        self.update_object_list.set(True)
        self.update_window()

    def plot(self):
        self.ax1.cla()

        try:
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
                     self.target_ra.get(), self.target_dec.get(), self.obs_year_month.get(), self.ax1, self.skyobject.get(),
                     self.observatory.get())

        except:
            pass

        self.canvas.draw()

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
                popt, pcov = plc.curve_fit(target_horizon_diference, [0], [0], p0=[sidereal_time - 0.1])
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








