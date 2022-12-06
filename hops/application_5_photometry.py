
import os
import glob
import numpy as np
import shutil
import warnings
import matplotlib
import webbrowser
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import hops.pylightcurve41 as plc

from astropy.io import fits as pf
from matplotlib.cm import Greys, Greys_r


from hops.application_windows import MainWindow, AddOnWindow
from hops.hops_tools.images_find_stars import *


class PhotometryWindow(MainWindow):

    def __init__(self, log):

        MainWindow.__init__(self, log, name='HOPS - Photometry', position=2)

        # set variables, create and place widgets, main window

        self.advanced_options_window = AddOnWindow(self, name='Advanced settings', position=2)
        self.show_advanced_options_button = self.Button(text='Advanced settings', command=self.advanced_options_window.show)

        self.use_variable_aperture = self.advanced_options_window.CheckButton(
            text='Vary the aperture size proportionally to the variations of the PSF size.',
            initial=self.log.get_param('use_variable_aperture'),
            command=self.update_window
        )

        self.use_geometric_center = self.advanced_options_window.CheckButton(
            text='Align the aperture with the geometric center instead of the PSF peak.',
            initial=self.log.get_param('use_geometric_center'),
            command=self.update_window
        )

        self.sky_inner_aperture = self.advanced_options_window.Entry(
            value=self.log.get_param('sky_inner_aperture'),
            instance=float,
            command=self.update_window
        )

        self.sky_outer_aperture = self.advanced_options_window.Entry(
            value=self.log.get_param('sky_outer_aperture'),
            instance=float,
            command=self.update_window
        )

        self.saturation = self.advanced_options_window.Entry(
            value=self.log.get_param('saturation'),
            instance=float,
            command=self.update_window
        )

        self.advanced_options_window.setup_window([
            [],
            [
                [self.advanced_options_window.Label(text='Inner sky-ring radius, relatively to the aperture'), 0],
                [self.sky_inner_aperture, 1],
                [self.advanced_options_window.Label(text='(default = 1.7)'), 2]
             ],
            [
                [self.advanced_options_window.Label(text='Outer sky-ring radius, relatively to the aperture'), 0],
                [self.sky_outer_aperture, 1],
                [self.advanced_options_window.Label(text='(default = 2.4)'), 2]
            ],
            [[self.advanced_options_window.Label(text='Note: the bright pixels inside the sky-ring\nare NOT taken into account for the sky background estimation.'), 0, 3]],
            [],
            [
                [self.advanced_options_window.Label(text='Saturation warning limit, relatively to the full-well depth'), 0],
                [self.saturation, 1],
                [self.advanced_options_window.Label(text='(default = 0.95)'), 2]],
            [],
            [[self.use_variable_aperture, 0, 3]],
            [[self.use_geometric_center, 0, 3]],
            [],
        ])

        self.bin_fits = self.log.get_param('bin_fits')
        self.burn_limit = self.log.get_param('burn_limit')
        self.max_targets = self.log.get_param('max_comparisons') + 1
        self.target_ra_dec = self.log.get_param('target_ra_dec')
        self.visible_fov_x_min = self.log.get_param('min_x')
        self.visible_fov_y_min = self.log.get_param('min_y')
        self.visible_fov_x_max = self.log.get_param('max_x')
        self.visible_fov_y_max = self.log.get_param('max_y')

        self.all_frames = plc.open_dict(self.log.all_frames)
        self.science_files = []
        self.exposure_time = []

        for science_file in self.all_frames:
            if not self.all_frames[science_file][self.log.skip_key]:
                self.science_files.append([self.all_frames[science_file][self.log.time_key], science_file])
                self.exposure_time.append(self.all_frames[science_file][self.log.get_param('exposure_time_key')])

        self.science_files.sort()
        self.science_files = [ff[1] for ff in self.science_files]
        self.exposure_time = np.median(self.exposure_time)

        self.all_stars = plc.open_dict(self.log.all_stars)
        self.in_fov = self.all_stars['in_fov']
        self.all_stars = self.all_stars['all_stars']

        self.fits = plc.open_fits(os.path.join(self.log.reduction_directory, self.science_files[0]))[1]

        self.show_good_comparisons = self.CheckButton(text='Show stars with flux similar\nto the target (+/- %):',
                                                      initial=self.log.get_param('show_good_comparisons'),
                                                      command=self.update_window)
        self.show_good_comparisons_percent = self.DropDown(initial=float(self.log.get_param('show_good_comparisons_percent')),
                                                           options=[10, 20, 30, 40, 50], instance=float,
                                                           command=self.update_window, width=8)

        self.targets_indication = self.IntVar(0)
        self.targets_indication_entry = [
            self.Radiobutton(text='      Target           ', variable=self.targets_indication, value=0)]

        try:
            self.targets_x_position = [self.Label(text=self.log.get_param('target_x_position'), instance=float)]
        except:
            self.targets_x_position = [self.Label(text=0, instance=float)]
        try:
            self.targets_y_position = [self.Label(text=self.log.get_param('target_y_position'), instance=float)]
        except:
            self.targets_y_position = [self.Label(text=0.0, instance=float)]
        self.targets_peak = [self.Label(text=0.0, instance=float)]
        self.targets_max_hwhm = [self.Label(text=0.0, instance=float)]
        self.targets_total_flux = [self.Label(text=0.0, instance=float)]
        try:
            self.targets_aperture = [self.Entry(value=self.log.get_param('target_aperture'), instance=float, command=self.update_window)]
        except:
            self.targets_aperture = [self.Entry(value=0.0, instance=float, command=self.update_window)]
        self.targets_clear = [self.Label(text='')]
        self.targets_warning = [self.Label()]
        self.targets_flux = [0.0]
        self.targets_wx = [0.0]
        self.targets_wy = [0.0]

        for comparison in range(self.log.get_param('max_comparisons')):

            self.targets_indication_entry.append(self.Radiobutton(text='Comparison {0}     '.format(comparison + 1),
                                                                  variable=self.targets_indication,
                                                                  value=comparison + 1))
            try:
                self.targets_x_position.append(self.Label(text=self.log.get_param('comparison_{0}_x_position'.format(comparison + 1)), instance=float))
            except:
                self.targets_x_position.append(self.Label(text=0.0, instance=float))
            try:
                self.targets_y_position.append(self.Label(text=self.log.get_param('comparison_{0}_y_position'.format(comparison + 1)), instance=float))
            except:
                self.targets_y_position.append(self.Label(text=0.0, instance=float))

            self.targets_peak.append(self.Label(text=0.0, instance=float))
            self.targets_max_hwhm.append(self.Label(text=0.0, instance=float))
            self.targets_total_flux.append(self.Label(text=0.0, instance=float))
            try:
                self.targets_aperture.append(self.Entry(value=self.log.get_param('comparison_{0}_aperture'.format(comparison + 1)),
                                                        instance=float, command=self.update_window))
            except:
                self.targets_aperture.append(self.Entry(value=0.0, instance=float, command=self.update_window))
            self.targets_clear.append(self.Button(text='CLEAR', command=self.clear_target(comparison + 1)))
            self.targets_warning.append(self.Label())
            self.targets_flux.append(0)
            self.targets_wx.append(0)
            self.targets_wy.append(0)

            if self.targets_indication.get() == 0 and self.targets_x_position[-1].get() == 0:
                self.targets_indication.set(comparison + 1)

        if self.targets_x_position[0].get() == 0:
            self.targets_indication.set(0)

        self.photometry_button = self.Button(text='RUN PHOTOMETRY', command=self.run_photometry, bg='green', highlightbackground='green')

        # plot

        self.fits_figure = self.FitsWindow(input=self.fits, input_name=self.science_files[0],
                                           input_options=self.log.get_param('photometry_fov_options'),
                                           show_controls=True, show_axes=True,
                                           subplots_adjust=(0.07, 0.99, 0.05, 0.85))

        self.fits_figure.canvas.callbacks.connect('button_press_event', self.update_window)

        self.circles_radius = 20 * self.fits.header[self.log.psf_key]

        self.fits_figure.ax.add_patch(mpatches.Rectangle((self.visible_fov_x_min + 1, self.visible_fov_y_min + 1),
                                             self.visible_fov_x_max - self.visible_fov_x_min - 2,
                                             self.visible_fov_y_max - self.visible_fov_y_min - 2,
                                             ec='r', fill=False, label='Available FOV'))

        if self.targets_x_position[0].get() ** 2 + self.targets_y_position[0].get() ** 2 == 0:
            self.targets_box = [[mpatches.Circle((-1000, -1000),
                                                self.targets_aperture[0].get(),
                                                ec='r', fill=False) for ff in range(3)]]
            self.targets_text = [self.fits_figure.ax.text(-1000,
                                                          -1000 - self.targets_aperture[
                                                              0].get() - 1, 'T',
                                                          color='r', fontsize=15, va='top')]
        else:
            self.targets_box = [[mpatches.Circle((self.targets_x_position[0].get(), self.targets_y_position[0].get()),
                                                self.targets_aperture[0].get(),
                                                ec='r', fill=False) for ff in range(3)]]
            self.targets_text = [self.fits_figure.ax.text(self.targets_x_position[0].get(),
                                                          self.targets_y_position[0].get() - self.targets_aperture[
                                                              0].get() - 1, 'T',
                                                          color='r', fontsize=15, va='top')]

        for comparison in range(self.log.get_param('max_comparisons')):
            if self.targets_x_position[comparison + 1].get() ** 2 + self.targets_y_position[comparison + 1].get() ** 2 == 0:
                self.targets_box.append([mpatches.Circle((-1000, -1000),
                                                         self.targets_aperture[comparison + 1].get(),
                                                         ec='#07fefc', fill=False) for ff in range(3)])
                self.targets_text.append(self.fits_figure.ax.text(-1000
                                                                  + self.targets_aperture[comparison + 1].get() + 1,
                                                                  -1000
                                                                  - self.targets_aperture[comparison + 1].get() - 1,
                                                                  'C{0}'.format(comparison + 1), color='#07fefc',
                                                                  fontsize=15, va='top'))
            else:
                self.targets_box.append([mpatches.Circle((self.targets_x_position[comparison + 1].get(),
                                                          self.targets_y_position[comparison + 1].get()),
                                                         self.targets_aperture[comparison + 1].get(),
                                                         ec='#07fefc', fill=False) for ff in range(3)])
                self.targets_text.append(self.fits_figure.ax.text(self.targets_x_position[comparison + 1].get()
                                                                  + self.targets_aperture[comparison + 1].get() + 1,
                                                                  self.targets_y_position[comparison + 1].get()
                                                                  - self.targets_aperture[comparison + 1].get() - 1,
                                                                  'C{0}'.format(comparison + 1), color='#07fefc',
                                                                  fontsize=15, va='top'))

        self.good_comps_boxes1 = []
        for comparison in range(100):
            if comparison == 0:
                circle1 = mpatches.Rectangle((-1000, -1000), self.circles_radius, self.circles_radius, ec='y', fill=False,
                                             label='Stars of similar flux to the target')
            else:
                circle1 = mpatches.Rectangle((-1000, -1000), self.circles_radius, self.circles_radius, ec='y', fill=False)

            self.good_comps_boxes1.append(circle1)

        for box in self.targets_box:
            for ff in range(3):
                self.fits_figure.ax.add_patch(box[ff])

        for circle1 in self.good_comps_boxes1:
            self.fits_figure.ax.add_patch(circle1)

        self.fits_figure.ax.legend(loc=(0, 1.01))

        photometry_folders = (glob.glob(os.path.join('{0}*'.format(self.log.photometry_directory_base))))

        def photometry_order(path):
            try:
                return float(path.split('_')[1])
            except:
                return 1

        photometry_folders = sorted(photometry_folders, key=lambda x: photometry_order(x))
        photometry_folders = ['Load options from previous run'] + photometry_folders
        self.photometry_folder_to_load = self.DropDown(initial='Load options from previous run',
                                               options=photometry_folders,
                                               width=40, command=self.load_options)

        self.run_message = self.Label(text='')

        self.proceed_button = self.Button(text='PROCEED TO FITTING', command=self.proceed)
        if self.log.get_param('photometry_complete'):
            self.proceed_button.activate()
        else:
            self.proceed_button.disable()

        self.save_and_return_button = self.Button(text='SAVE OPTIONS & RETURN TO MAIN MENU', command=self.save_and_return)

        setup_list = [
            [[self.show_good_comparisons, 0],
             [self.show_good_comparisons_percent, 1],
             [self.Button(text='Check SIMBAD', command=self.openweb_simbad), 2],
             [self.Label(text="Remember, the best comparison stars need to be:\n"
                              "a) close to your target, b) of similar magnitude to the target,\n"
                              "c) of similar colour to the target, d) photometrically stable, i.e. "
                              "not variables!"), 3, 9]
             ],
            [[self.fits_figure, 0, 3, 24],],
            [[self.Label(text='X'), 5], [self.Label(text='Y'), 6], [self.Label(text='Total\ncounts'), 7],
             [self.Label(text='Max\ncounts'), 8], [self.Label(text='Max\nHWHM'), 9],
             [self.Label(text='Aperture\nradius (>1.5)'), 10]],
        ]

        for target in range(self.max_targets):

            setup_list.append([[self.targets_indication_entry[target], 3],
                               [self.targets_clear[target], 4],
                               [self.targets_x_position[target], 5],
                               [self.targets_y_position[target], 6],
                               [self.targets_total_flux[target], 7],
                               [self.targets_peak[target], 8],
                               [self.targets_max_hwhm[target], 9],
                               [self.targets_aperture[target], 10],
                               [self.targets_warning[target], 11]])

        setup_list += [
            [[self.show_advanced_options_button, 7, 4]],
            [],
            [[self.photometry_folder_to_load, 3, 8]],
            [[self.Button(text='RETURN TO MAIN MENU', command=self.close), 3, 2],
             [self.save_and_return_button, 5, 6]],
            [[self.photometry_button, 3, 2],
             [self.proceed_button, 5, 6]],
            [[self.run_message, 3, 8]],
            []
        ]

        self.setup_window(setup_list, entries_wd=5)

        self.update_window(None)

    def load_options(self):

        try:

            other_log = self.log.open_yaml(os.path.join(self.photometry_folder_to_load.get(), 'log.yaml'))

            self.use_variable_aperture.set(other_log['use_variable_aperture'])
            self.use_geometric_center.set(other_log['use_geometric_center'])
            self.sky_inner_aperture.set(other_log['sky_inner_aperture'])
            self.sky_outer_aperture.set(other_log['sky_outer_aperture'])
            self.saturation.set(other_log['saturation'])
            self.show_good_comparisons.set(other_log['show_good_comparisons'])
            self.show_good_comparisons_percent.set(other_log['show_good_comparisons_percent'])

            self.targets_indication.set(0)
            try:
                self.targets_x_position[0].set(other_log['target_x_position'])
                self.targets_y_position[0].set(other_log['target_y_position'])
                self.targets_aperture[0].set(other_log['target_aperture'])
            except:
                self.targets_x_position[0].set(0)
                self.targets_y_position[0].set(0)
                self.targets_aperture[0].set(0)

            for comparison in range(self.log.get_param('max_comparisons')):

                try:
                    self.targets_x_position[1 + comparison].set(other_log['comparison_{0}_x_position'.format(comparison + 1)])
                    self.targets_y_position[1 + comparison].set(other_log['comparison_{0}_y_position'.format(comparison + 1)])
                    self.targets_aperture[1 + comparison].set(other_log['comparison_{0}_aperture'.format(comparison + 1)])
                except:
                    self.targets_x_position[1 + comparison].set(0)
                    self.targets_y_position[1 + comparison].set(0)
                    self.targets_aperture[1 + comparison].set(0)

            for comparison in range(self.log.get_param('max_comparisons')):
                if self.targets_x_position[comparison].get() == 0:
                    self.targets_indication.set(comparison)
                    break

            self.update_window(None)

        except FileNotFoundError:
            pass

    def get_star(self, x, y):

        star = find_single_star(
            self.fits.data, x, y,
            mean=self.fits.header[self.log.mean_key], std=self.fits.header[self.log.std_key],
            burn_limit=self.burn_limit * self.bin_fits * self.bin_fits,
            star_std=self.fits.header[self.log.psf_key])

        if not star:
            return None

        else:
            app = round(3 * self.fits.header[self.log.psf_key], 1)
            x = round(star[0], 1)
            y = round(star[1], 1)

            x1 = int(star[0] - app)
            x2 = x1 + int(2 * app) + 2
            y1 = int(star[1] - app)
            y2 = y1 + int(2 * app) + 2

            peak = round(np.max(self.fits.data[y1:y2, x1:x2]), 1)

            flux = round(2 * np.pi * star[2] * star[4] * star[5], 1)

            max_hwhm = round(0.5 * 2.355 * max(star[4], star[5]), 1)

            x1 = int(star[0] - 3 * self.fits.header[self.log.psf_key])
            x2 = x1 + int(6 * self.fits.header[self.log.psf_key]) + 1
            y1 = int(star[1] - 3 * self.fits.header[self.log.psf_key])
            y2 = y1 + int(6 * self.fits.header[self.log.psf_key]) + 1

            area = self.fits.data[y1:y2, x1:x2]

            area_x, area_y = np.meshgrid(
                np.arange(x1, x2) + 0.5,
                np.arange(y1, y2) + 0.5)

            wx = (round(np.sum(area_x * area) / np.sum(area), 1))
            wy = (round(np.sum(area_y * area) / np.sum(area), 1))

            return x, y, app, peak, flux, max_hwhm, wx, wy

    def update_window(self, event=None):

        photometry_active = True

        if isinstance(event, matplotlib.backend_bases.MouseEvent):

            if event.inaxes is None:
                return None

            if event.dblclick:

                star = self.get_star(event.xdata, event.ydata)

                if not star:
                    self.showinfo('Star not acceptable.',
                                  'Star could not be located or it is saturated.')
                    return None

                else:

                    x, y, app, peak, flux, max_hwhm, wx, wy = star

                    if (x < self.visible_fov_x_min or x > self.visible_fov_x_max
                            or y < self.visible_fov_y_min or y > self.visible_fov_y_max):

                        self.showinfo('Star not acceptable.', 'Star moves outside the FOV later.')
                        return None

                    else:

                        self.targets_x_position[self.targets_indication.get()].set(x)
                        self.targets_y_position[self.targets_indication.get()].set(y)
                        self.targets_aperture[self.targets_indication.get()].set(app)
                        self.targets_peak[self.targets_indication.get()].set(peak)
                        self.targets_flux[self.targets_indication.get()] = flux
                        self.targets_wx[self.targets_indication.get()] = wx
                        self.targets_wy[self.targets_indication.get()] = wy
                        self.targets_total_flux[self.targets_indication.get()].set(flux)
                        self.targets_max_hwhm[self.targets_indication.get()].set(max_hwhm)

            else:
                return None

        for i_target in range(self.max_targets):

            if 0 in [self.targets_x_position[i_target].get(), self.targets_y_position[i_target].get()]:

                self.targets_aperture[i_target].disable()

                self.targets_box[i_target][0].set_center((-10000, -10000))
                self.targets_box[i_target][1].set_center((-10000, -10000))
                self.targets_box[i_target][2].set_center((-10000, -10000))
                self.targets_text[i_target].set_x(-10000)
                self.targets_text[i_target].set_y(-10000)

                self.targets_x_position[i_target].set(0)
                self.targets_y_position[i_target].set(0)
                self.targets_aperture[i_target].set(0)
                self.targets_peak[i_target].set(0)
                self.targets_flux[i_target] = 0
                self.targets_wx[i_target] = 0
                self.targets_wy[i_target] = 0
                self.targets_total_flux[i_target].set(0)
                self.targets_max_hwhm[i_target].set(0)
                self.targets_warning[i_target].set('')

                if i_target == 0:
                    for num in range(100):
                        self.good_comps_boxes1[num].set_xy((-1000, -1000))

            else:

                self.targets_aperture[i_target].activate()

                if self.targets_flux[i_target] == 0:

                    star = self.get_star(self.targets_x_position[i_target].get(),
                                         self.targets_y_position[i_target].get())

                    if not star:
                        xx = self.clear_target(i_target)
                        xx()

                    else:

                        x, y, app, peak, flux, max_hwhm, wx, wy = star

                        self.targets_x_position[i_target].set(x)
                        self.targets_y_position[i_target].set(y)
                        if self.targets_aperture[i_target].get() == 0:
                            self.targets_aperture[i_target].set(app)
                        self.targets_peak[i_target].set(peak)
                        self.targets_flux[i_target] = flux
                        self.targets_wx[i_target] = wx
                        self.targets_wy[i_target] = wy
                        self.targets_total_flux[i_target].set(flux)
                        self.targets_max_hwhm[i_target].set(max_hwhm)

                if self.targets_aperture[i_target].get() < 1.5:
                    self.targets_warning[i_target].set('Ap. too small')
                else:
                    self.targets_warning[i_target].set('')

                # for compatibility with older versions - end

                fov_options = self.fits_figure.get_fov_options()
                text_y_drift = [-1, 1][int(fov_options[3])]
                text_x_drift = [1, -1][int(fov_options[4])]

                if i_target == 0:

                    good_comps = []

                    for j, comp_star in enumerate(self.all_stars):
                        if np.sqrt((comp_star[0] - self.targets_x_position[i_target].get())**2 +
                                   (comp_star[1] - self.targets_y_position[i_target].get())**2) > self.fits.header[self.log.psf_key]:
                            if self.in_fov[j]:
                                if comp_star[-1] < (1 + self.show_good_comparisons_percent.get() / 100) * self.targets_flux[i_target]:
                                    if comp_star[-1] > (1 - self.show_good_comparisons_percent.get() / 100) * self.targets_flux[i_target]:
                                        good_comps.append(comp_star)

                    good_comps = sorted(good_comps,
                                        key=lambda x: np.sqrt((x[0] -
                                                               self.targets_x_position[i_target].get()) ** 2 +
                                                              (x[1] -
                                                               self.targets_y_position[i_target].get()) ** 2))

                    for num in range(100):

                        if num < len(good_comps):
                            self.good_comps_boxes1[num].set_xy((good_comps[num][0] - self.circles_radius / 2,
                                                                good_comps[num][1] - self.circles_radius / 2))
                        else:
                            self.good_comps_boxes1[num].set_xy((-1000, -1000))

                        if not self.show_good_comparisons.get():
                            self.good_comps_boxes1[num].set_alpha(0)
                        else:
                            self.good_comps_boxes1[num].set_alpha(1)

                else:

                    if 0 not in [self.targets_x_position[0].get(), self.targets_y_position[0].get()]:
                        if self.targets_flux[i_target] > 2 * self.targets_flux[0]:
                            if self.targets_warning[i_target].get() == '':
                                self.targets_warning[i_target].set('Comp. too bright')
                            else:
                                self.targets_warning[i_target].set(
                                    ', '.join(self.targets_warning[i_target].get().split(', ') + ['Comp. too bright']))
                        elif self.targets_flux[i_target] < 0.5 * self.targets_flux[0]:
                            if self.targets_warning[i_target].get() == '':
                                self.targets_warning[i_target].set('Comp. too faint')
                            else:
                                self.targets_warning[i_target].set(
                                    ', '.join(self.targets_warning[i_target].get().split(', ') + ['Comp. too faint']))

                for ff in range(3):
                    if self.use_geometric_center.get():
                        self.targets_box[i_target][ff].set_center((self.targets_wx[i_target],
                                                                   self.targets_wy[i_target]))
                    else:
                        self.targets_box[i_target][ff].set_center((self.targets_x_position[i_target].get(),
                                                                 self.targets_y_position[i_target].get()))

                self.targets_box[i_target][0].set_radius(self.targets_aperture[i_target].get())
                self.targets_box[i_target][1].set_radius(self.targets_aperture[i_target].get() * self.sky_inner_aperture.get())
                self.targets_box[i_target][1].set_linestyle(":")
                self.targets_box[i_target][2].set_radius(self.targets_aperture[i_target].get() * self.sky_outer_aperture.get())
                self.targets_box[i_target][2].set_linestyle(":")

                self.targets_text[i_target].set_x(self.targets_x_position[i_target].get() + text_x_drift * self.targets_aperture[i_target].get())
                self.targets_text[i_target].set_y(self.targets_y_position[i_target].get() + text_y_drift * self.targets_aperture[i_target].get())

                if self.targets_peak[i_target].get() > self.saturation.get() * self.burn_limit * self.bin_fits * self.bin_fits:
                    if self.targets_warning[i_target].get() == '':
                        self.targets_warning[i_target].set('Star close to saturation')
                    else:
                        self.targets_warning[i_target].set(
                            ', '.join(self.targets_warning[i_target].get().split(', ') + ['Star close to saturation']))

        self.run_message.set('')
        for i_target in range(self.max_targets):
            if 0 not in [self.targets_x_position[i_target].get(), self.targets_y_position[i_target].get()]:
                if self.targets_aperture[i_target].get() < 1.5:
                    photometry_active = False
                    self.run_message.set('All aperture radii must be larget than 1.5 pixels')

        if 0 in [self.targets_x_position[0].get(), self.targets_y_position[0].get(),
                 self.targets_x_position[1].get(), self.targets_y_position[1].get()]:
            self.run_message.set('You need to select the target and at least one comparison star to proceed.')
            photometry_active = False

        if self.sky_inner_aperture.get() <= 1.0 or self.sky_outer_aperture.get() <= 1.0:
            self.run_message.set('Inner and outer sky-ring radii must be larger than 1.')
            photometry_active = False

        if self.sky_inner_aperture.get() >= self.sky_outer_aperture.get():
            self.run_message.set('Inner sky-ring radius must be smaller than the outer sky-ring radius.')
            photometry_active = False

        if self.saturation.get() <= 0:
            self.run_message.set('Saturation limit must be larger than 0.')
            photometry_active = False

        if self.saturation.get() >= 1:
            self.run_message.set('Saturation limit must be lower than 1.')
            photometry_active = False

        if photometry_active:
            self.photometry_button.activate()
            self.save_and_return_button.activate()
        else:
            self.photometry_button.disable()
            self.save_and_return_button.disable()

        for i_target in range(self.max_targets):
            self.targets_indication_entry[i_target]['state'] = self.DISABLED
        for i_target in range(self.max_targets):
            if 0 not in [self.targets_x_position[i_target].get(), self.targets_y_position[i_target].get()]:
                self.targets_indication_entry[i_target]['state'] = self.NORMAL
            else:
                self.targets_indication_entry[i_target]['state'] = self.NORMAL
                break

        self.fits_figure.draw()

    # define actions for the different buttons, including calls to the function that updates the window

    def clear_target(self, num):

        def clear_target_num():

            self.targets_x_position[num].set(0.0)
            self.targets_y_position[num].set(0.0)
            self.targets_aperture[num].set(0.0)
            self.targets_peak[num].set(0.0)
            self.targets_max_hwhm[num].set(0.0)
            self.targets_total_flux[num].set(0.0)
            self.targets_flux[num] = 0.0
            self.targets_warning[num].set('')

            for i_target in range(1, self.max_targets - 1):
                if self.targets_x_position[i_target].get() == 0 and self.targets_x_position[i_target + 1].get() != 0:
                    self.targets_x_position[i_target].set(self.targets_x_position[i_target + 1].get())
                    self.targets_y_position[i_target].set(self.targets_y_position[i_target + 1].get())
                    self.targets_aperture[i_target].set(self.targets_aperture[i_target + 1].get())
                    self.targets_wx[i_target] = self.targets_wx[i_target + 1]
                    self.targets_wy[i_target] = self.targets_wy[i_target + 1]
                    self.targets_peak[i_target].set(self.targets_peak[i_target + 1].get())
                    self.targets_max_hwhm[i_target].set(self.targets_max_hwhm[i_target + 1].get())
                    self.targets_total_flux[i_target].set(self.targets_total_flux[i_target + 1].get())
                    self.targets_flux[i_target] = self.targets_flux[i_target + 1]
                    self.targets_warning[i_target].set(self.targets_warning[i_target + 1].get())

                    self.targets_x_position[i_target + 1].set(0.0)
                    self.targets_y_position[i_target + 1].set(0.0)
                    self.targets_aperture[i_target + 1].set(0.0)
                    self.targets_wx[i_target + 1] = 0
                    self.targets_wy[i_target + 1] = 0
                    self.targets_peak[i_target + 1].set(0.0)
                    self.targets_max_hwhm[i_target + 1].set(0.0)
                    self.targets_total_flux[i_target + 1].set(0.0)
                    self.targets_flux[i_target + 1] = 0.0
                    self.targets_warning[i_target + 1].set('')

            if self.targets_indication.get() > num:
                self.targets_indication.set(self.targets_indication.get() - 1)

            self.update_window(None)

        return clear_target_num

    def openweb_simbad(self):

        radec_string = self.log.get_param('target_ra_dec').replace('+', '%2B').replace(' ', '+')
        webbrowser.open("http://simbad.u-strasbg.fr/simbad/sim-coo?Coord={0}&CooFrame=ICRS&"
                        "CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=20&"
                        "Radius.unit=arcmin&submit=submit+query&CoordList=".format(radec_string), new=1)

    def save_settings(self):

        self.log.set_param('show_good_comparisons', self.show_good_comparisons.get())
        self.log.set_param('show_good_comparisons_percent', self.show_good_comparisons_percent.get())
        self.log.set_param('show_good_comparisons_percent', self.show_good_comparisons_percent.get())
        self.log.set_param('sky_inner_aperture', self.sky_inner_aperture.get())
        self.log.set_param('sky_outer_aperture', self.sky_outer_aperture.get())
        self.log.set_param('use_variable_aperture', self.use_variable_aperture.get())
        self.log.set_param('use_geometric_center', self.use_geometric_center.get())
        self.log.set_param('photometry_fov_options', self.fits_figure.get_fov_options())

        # targets

        targets = 1

        self.log.set_param('target_x_position', self.targets_x_position[0].get())
        self.log.set_param('target_y_position', self.targets_y_position[0].get())
        self.log.set_param('target_aperture', self.targets_aperture[0].get())
        target_polar = plc.cartesian_to_polar(self.targets_x_position[0].get(), self.targets_y_position[0].get(),
                                              self.fits.header[self.log.align_x0_key],
                                              self.fits.header[self.log.align_y0_key])
        self.log.set_param('target_r_position', float(target_polar[0]))
        self.log.set_param('target_u_position', float(target_polar[1]))

        for i_comparison in range(self.log.get_param('max_comparisons')):
            self.log.set_param('comparison_{0}_x_position'.format(i_comparison + 1),
                               self.targets_x_position[i_comparison + 1].get())
            self.log.set_param('comparison_{0}_y_position'.format(i_comparison + 1),
                               self.targets_y_position[i_comparison + 1].get())
            self.log.set_param('comparison_{0}_aperture'.format(i_comparison + 1),
                               self.targets_aperture[i_comparison + 1].get())

            if 0 not in [self.targets_x_position[i_comparison + 1].get(),
                         self.targets_y_position[i_comparison + 1].get()]:

                target_polar = plc.cartesian_to_polar(self.targets_x_position[i_comparison + 1].get(),
                                                      self.targets_y_position[i_comparison + 1].get(),
                                                      self.fits.header[self.log.align_x0_key],
                                                      self.fits.header[self.log.align_y0_key])
                targets += 1

            else:

                target_polar = [0, 0]

            self.log.set_param('comparison_{0}_r_position'.format(i_comparison + 1), float(target_polar[0]))
            self.log.set_param('comparison_{0}_u_position'.format(i_comparison + 1), float(target_polar[1]))

        self.targets_indication.set(targets)

        # save log

        self.log.save_local_log()

    def save_and_return(self):

        self.save_settings()

        self.close()

    def run_photometry(self):

        self.save_settings()

        self.log.set_param('proceed', 'run_photometry')
        self.close()

    def proceed(self):
        self.log.set_param('proceed', True)
        self.close()


class PhotometryProgressWindow(MainWindow):

    def __init__(self, log):

        MainWindow.__init__(self, log, name='HOPS - Photometry progress', position=2)

        # load variabvles

        self.bin_fits = self.log.get_param('bin_fits')
        self.burn_limit = self.log.get_param('burn_limit')
        self.saturation = self.log.get_param('saturation')
        self.max_targets = self.log.get_param('max_comparisons') + 1
        self.target_ra_dec = self.log.get_param('target_ra_dec')
        self.visible_fov_x_min = self.log.get_param('min_x')
        self.visible_fov_y_min = self.log.get_param('min_y')
        self.visible_fov_x_max = self.log.get_param('max_x')
        self.visible_fov_y_max = self.log.get_param('max_y')
        self.use_geometric_center = self.log.get_param('use_geometric_center')
        self.use_variable_aperture = self.log.get_param('use_variable_aperture')
        self.sky_inner_aperture = self.log.get_param('sky_inner_aperture')
        self.sky_outer_aperture = self.log.get_param('sky_outer_aperture')

        # load science files and targets

        self.all_frames = plc.open_dict(self.log.all_frames)
        self.science_files = []
        self.exposure_time = []

        for science_file in self.all_frames:
            if not self.all_frames[science_file][self.log.skip_key]:
                self.science_files.append(science_file)
                self.exposure_time.append(self.all_frames[science_file][self.log.get_param('exposure_time_key')])

        self.exposure_time = np.median(self.exposure_time)
        self.science_files.sort()

        self.first_fits = plc.open_fits(os.path.join(self.log.reduction_directory, self.science_files[0]))

        self.jd = []
        self.file_names = []
        self.psf_ratio = []
        for science_file in self.science_files:
            metadata = self.all_frames[science_file]
            self.jd.append(metadata[self.log.time_key])
            self.file_names.append(science_file)
            self.psf_ratio.append(metadata[self.log.psf_key])

        if self.use_variable_aperture:
            self.psf_ratio = np.array(self.psf_ratio) / self.psf_ratio[0]
        else:

            self.psf_ratio = np.ones_like(self.psf_ratio)

        self.targets_aperture = [self.log.get_param('target_aperture')]
        self.targets_x_position = [self.log.get_param('target_x_position')]
        self.targets_y_position = [self.log.get_param('target_y_position')]
        self.targets_r_position = [self.log.get_param('target_r_position')]
        self.targets_u_position = [self.log.get_param('target_u_position')]

        for i_comparison in range(self.log.get_param('max_comparisons')):
            if 0 not in [self.log.get_param('comparison_{0}_x_position'.format(i_comparison + 1)),
                         self.log.get_param('comparison_{0}_y_position'.format(i_comparison + 1))]:
                self.targets_aperture.append(self.log.get_param('comparison_{0}_aperture'.format(i_comparison + 1)))
                self.targets_x_position.append(self.log.get_param('comparison_{0}_x_position'.format(i_comparison + 1)))
                self.targets_y_position.append(self.log.get_param('comparison_{0}_y_position'.format(i_comparison + 1)))
                self.targets_r_position.append(self.log.get_param('comparison_{0}_r_position'.format(i_comparison + 1)))
                self.targets_u_position.append(self.log.get_param('comparison_{0}_u_position'.format(i_comparison + 1)))

        self.targets = len(self.targets_aperture)
        self.max_aperture = np.max(self.targets_aperture)

        # results placeholder

        self.results = {
            'jd': np.array(self.jd),
            'file_names': np.array(self.file_names),
            'gauss_x_position': np.array([[np.nan for fff in self.science_files] for ff in range(self.targets)]),
            'gauss_y_position': np.array([[np.nan for fff in self.science_files] for ff in range(self.targets)]),
            'gauss_x_std': np.array([[np.nan for fff in self.science_files] for ff in range(self.targets)]),
            'gauss_y_std': np.array([[np.nan for fff in self.science_files] for ff in range(self.targets)]),
            'gauss_flux': np.array([[np.nan for fff in self.science_files] for ff in range(self.targets)]),
            'gauss_flux_error': np.array([[np.nan for fff in self.science_files] for ff in range(self.targets)]),
            'gauss_sky': np.array([[np.nan for fff in self.science_files] for ff in range(self.targets)]),
            'gauss_sky_error': np.array([[np.nan for fff in self.science_files] for ff in range(self.targets)]),
            'aperture_x_position': np.array([[np.nan for fff in self.science_files] for ff in range(self.targets)]),
            'aperture_y_position': np.array([[np.nan for fff in self.science_files] for ff in range(self.targets)]),
            'aperture_flux': np.array([[np.nan for fff in self.science_files] for ff in range(self.targets)]),
            'aperture_flux_error': np.array([[np.nan for fff in self.science_files] for ff in range(self.targets)]),
            'aperture_sky': np.array([[np.nan for fff in self.science_files] for ff in range(self.targets)]),
            'aperture_sky_error': np.array([[np.nan for fff in self.science_files] for ff in range(self.targets)]),
        }

        self.results_saved = False

        # progress figure

        self.progress_figure = self.FigureWindow()

        self.axes_fits = []
        self.circles_fits = []
        self.images_fits = []
        self.star_image_lim = int(1.5 * np.max(self.max_aperture) + 5)
        self.star_plot_lim = 1.5 * np.max(self.max_aperture)
        self.axes_lc = []
        self.images_text_lc = []
        self.images_gauss_lc = []
        self.images_aperture_lc = []
        self.axes_lc_ylims = [0.99, 1.01]

        xx = 1 / (self.targets + 2)
        gs = gridspec.GridSpec(self.targets, 4, self.progress_figure.figure, 0.01, xx, 0.99, 1 - xx, 1.0, 0.1)

        for target in range(self.targets):
            self.axes_fits.append(self.progress_figure.figure.add_subplot(gs[target, 0]))

            circle = mpatches.Circle((0, 0), self.targets_aperture[target], ec='r', fill=False)
            self.axes_fits[-1].add_patch(circle)
            self.circles_fits.append(circle)

            self.images_fits.append(self.axes_fits[-1].imshow(
                0 * self.first_fits[1].data[0: 2 * self.star_image_lim, 0: 2 * self.star_image_lim],
                origin='lower',
                extent=(- self.star_image_lim, self.star_image_lim, - self.star_image_lim, self.star_image_lim),
                cmap=Greys_r,
                vmin=self.first_fits[1].header[self.log.mean_key] + self.log.frame_low_std * self.first_fits[1].header[self.log.std_key],
                vmax=self.first_fits[1].header[self.log.mean_key] + self.log.frame_upper_std * self.first_fits[1].header[self.log.std_key]
            ))

            self.axes_fits[-1].axis('off')
            self.axes_fits[-1].set_xlim(- self.star_plot_lim, self.star_plot_lim)
            self.axes_fits[-1].set_ylim(- self.star_plot_lim, self.star_plot_lim)

            self.axes_lc.append(self.progress_figure.figure.add_subplot(gs[target, 1:]))

            self.images_text_lc.append(
                self.axes_lc[-1].text((self.results['jd'][-1] - self.results['jd'][0]) * 24 / 2, 1, ' ', fontsize=10,
                                      va='center', ha='center'))
            self.images_gauss_lc.append(
                self.axes_lc[-1].plot((np.array(self.results['jd']) - self.results['jd'][0]) * 24,
                                      self.results['gauss_flux'][0], 'ro', ms=3, label='PSF'))
            self.images_aperture_lc.append(
                self.axes_lc[-1].plot((np.array(self.results['jd']) - self.results['jd'][0]) * 24,
                                      self.results['aperture_flux'][0], 'ko', ms=3, label='Aperture'))

            self.axes_lc[-1].set_ylim(self.axes_lc_ylims[0], self.axes_lc_ylims[1])
            self.axes_lc[-1].set_xlim((-0.01) * 24, (self.results['jd'][-1] - self.results['jd'][0] + 0.01) * 24)
            if target == 0:
                self.axes_lc[-1].set_ylabel('T', fontsize=12)
            else:
                self.axes_lc[-1].set_ylabel('C{0}'.format(target), fontsize=12)

            if target != self.targets - 1:
                self.axes_lc[-1].tick_params(labelbottom=False, labelsize=8)
            else:
                self.axes_lc[-1].set_xlabel('Time (hours in observation)', fontsize=8)
                self.axes_lc[-1].tick_params(labelsize=8)

        self.axes_lc[0].legend(bbox_to_anchor=(1, 2.2), fontsize='x-small')

        if self.targets == 2:
            self.images_text_lc[1].set_text('Can\'t plot relative curve with only 1 active comparison')

        # progress bar and diagnostics

        self.progress_photometry = self.Progressbar(task="Photometry")

        self.target_label = self.Label(text='Target')

        self.progress_diagnostics = [self.target_label]
        self.progress_active = [self.CheckButton(text='',initial=1)]
        self.progress_active[0].disable()
        for target in range(1, self.targets):
            self.progress_diagnostics.append(self.Label(text='Comparison {0}'.format(target)))
            self.progress_active.append(self.CheckButton(text='', initial=1, command=self.replot_results))

        for i in self.progress_active:
            i.disable()

        # buttons

        self.stop_and_return_button = self.Button(
            text='STOP PHOTOMETRY', command=self.trigger_exit)
        self.save_button = self.Button(
            text='SAVE RESULTS', command=self.save)
        self.return_to_photomertry_button = self.Button(
            text='RETURN TO PHOTOMETRY MENU', command=self.close)
        self.proceed_button = self.Button(
            text='PROCEED TO FITTING MENU', command=self.proceed)

        # structure

        structure = [
            [[self.Label(text=' '), 0], [self.Label(text='Active'), 1], [self.progress_figure, 2, 2, self.targets + 2]],
            [[self.progress_diagnostics[0], 0], [self.progress_active[0], 1]],
        ]

        for i in range(1, self.targets):
            structure.append([[self.progress_diagnostics[i], 0], [self.progress_active[i], 1]])

        structure.append([[self.Label(text='...'), 0]])
        structure.append([])

        structure.append([[self.progress_photometry, 2], [self.stop_and_return_button, 3]])
        structure.append([[self.save_button, 0, 2], [self.return_to_photomertry_button, 2], [self.proceed_button, 3]])
        structure.append([])

        self.set_close_button_function(self.trigger_exit)

        self.setup_window(structure, main_font='Courier')

        self.stop_and_return_button.activate()
        self.save_button.disable()
        self.return_to_photomertry_button.disable()
        self.proceed_button.disable()

    def run_photometry(self):

        self.log.set_param('proceed', 'return_to_photometry')

        self.science_counter = 0
        self.target_fails = np.zeros(self.targets)

        self.after(self.measure)

    def measure(self):

        if self.exit:
            self.after(self.options)

        else:

            if self.science_counter == 0:
                self.progress_photometry.initiate(len(self.science_files))

            counter = self.science_counter
            science_file = self.science_files[counter]

            fits = pf.open(os.path.join(self.log.reduction_directory, science_file), memmap=False, mode='update')
            metadata = self.all_frames[science_file]

            ref_x_position = metadata[self.log.align_x0_key]
            ref_y_position = metadata[self.log.align_y0_key]
            ref_u_position = metadata[self.log.align_u0_key]

            for target in range(self.targets):

                if self.exit:
                    break

                star = find_single_star(
                    fits[1].data,
                    (ref_x_position + self.targets_r_position[target] *
                     np.cos(ref_u_position + self.targets_u_position[target])),
                    (ref_y_position + self.targets_r_position[target] *
                     np.sin(ref_u_position + self.targets_u_position[target])),
                    mean=fits[1].header[self.log.mean_key],
                    std=fits[1].header[self.log.std_key],
                    burn_limit=self.burn_limit * self.bin_fits * self.bin_fits,
                    star_std=fits[1].header[self.log.psf_key]
                )

                if star:

                    x_mean, y_mean, norm, floor, x_std, y_std, centroid_x, centroid_y = star

                    # save psf photometry

                    self.results['gauss_x_position'][target][counter] = x_mean
                    self.results['gauss_y_position'][target][counter] = y_mean
                    self.results['gauss_x_std'][target][counter] = x_std
                    self.results['gauss_y_std'][target][counter] = y_std
                    self.results['gauss_flux'][target][counter] = 2 * np.pi * norm * x_std * y_std
                    self.results['gauss_flux_error'][target][counter] = np.sqrt(
                        np.abs(2 * np.pi * norm * x_std * y_std) + 2 * np.abs(9 * floor * x_std * y_std))
                    self.results['gauss_sky'][target][counter] = floor
                    self.results['gauss_sky_error'][target][counter] = np.sqrt(floor)

                    try:

                        if self.use_geometric_center:
                            x1 = int(star[0] - 3 * fits[1].header[self.log.psf_key])
                            x2 = x1 + int(6 * fits[1].header[self.log.psf_key]) + 1
                            y1 = int(star[1] - 3 * fits[1].header[self.log.psf_key])
                            y2 = y1 + int(6 * fits[1].header[self.log.psf_key]) + 1

                            area = fits[1].data[y1:y2, x1:x2]

                            area_x, area_y = np.meshgrid(
                                np.arange(x1, x2) + 0.5,
                                np.arange(y1, y2) + 0.5)

                            x_mean = np.sum(area * area_x)/np.sum(area)
                            y_mean = np.sum(area * area_y)/np.sum(area)

                        # progress plot

                        plot_x_min = int(x_mean - self.star_image_lim)
                        plot_y_min = int(y_mean - self.star_image_lim)

                        self.images_fits[target].set_data(
                            fits[1].data[plot_y_min: plot_y_min + 2 * self.star_image_lim,
                            plot_x_min: plot_x_min + 2 * self.star_image_lim])
                        self.images_fits[target].set_clim(
                            fits[1].header[self.log.mean_key] + self.log.frame_low_std * fits[1].header[
                                self.log.std_key],
                            fits[1].header[self.log.mean_key] + self.log.frame_upper_std * fits[1].header[
                                self.log.std_key])
                        self.images_fits[target].set_extent((plot_x_min, plot_x_min + 2 * self.star_image_lim,
                                                             plot_y_min, plot_y_min + 2 * self.star_image_lim))

                        self.axes_fits[target].set_xlim(x_mean - self.star_plot_lim,
                                                        x_mean + self.star_plot_lim)
                        self.axes_fits[target].set_ylim(y_mean - self.star_plot_lim,
                                                        y_mean + self.star_plot_lim)

                        self.circles_fits[target].set_center((x_mean, y_mean))

                        self.circles_fits[target].set_radius(self.targets_aperture[target] * self.psf_ratio[counter])

                        # aperture photometry

                        sky_area_1 = int(round(
                            self.sky_inner_aperture * self.targets_aperture[target] * self.psf_ratio[
                                counter]))
                        sky_area_2 = int(round(
                            self.sky_outer_aperture * self.targets_aperture[target] * self.psf_ratio[
                                counter]))
                        sky_area_2 = max(sky_area_2, sky_area_1 + 3)

                        sky_area = fits[1].data[int(y_mean) - sky_area_2:int(y_mean) + sky_area_2 + 1,
                                   int(x_mean) - sky_area_2:int(x_mean) + sky_area_2 + 1]

                        sky_area = np.ones_like(sky_area) * sky_area

                        sky_center = int(len(sky_area) / 2)

                        sky_area[sky_center - sky_area_1:sky_center + sky_area_1 + 1,
                        sky_center - sky_area_1:sky_center + sky_area_1 + 1] = -100000

                        sky_area = sky_area[np.where((sky_area != -100000) &
                                                     (sky_area < fits[1].header[self.log.mean_key] + 3 * fits[1].header[
                                                         self.log.std_key]))]

                        sky_total = np.sum(sky_area)
                        sky = np.sum(sky_area) / sky_area.size
                        sky_error = np.sqrt(np.abs(sky_total)) / sky_area.size

                        flux_area = fits[1].data[
                                    int(y_mean - self.targets_aperture[target] * self.psf_ratio[counter]) - 2:
                                    int(y_mean + self.targets_aperture[target] * self.psf_ratio[counter]) + 3,
                                    int(x_mean - self.targets_aperture[target] * self.psf_ratio[counter]) - 2:
                                    int(x_mean + self.targets_aperture[target] * self.psf_ratio[counter]) + 3]

                        flux_area_x, flux_area_y = np.meshgrid(
                            np.arange(
                                max(0, int(x_mean - self.targets_aperture[target] * self.psf_ratio[counter]) - 2),
                                min(len(fits[1].data[0]),
                                    int(x_mean + self.targets_aperture[target] * self.psf_ratio[counter]) + 3),
                                1) + 0.5,
                            np.arange(
                                max(0, int(y_mean - self.targets_aperture[target] * self.psf_ratio[counter]) - 2),
                                min(len(fits[1].data),
                                    int(y_mean + self.targets_aperture[target] * self.psf_ratio[counter]) + 3),
                                1) + 0.5)

                        flux_pixels = np.concatenate(np.swapaxes([flux_area, flux_area_x, flux_area_y], 0, 2))

                        aperture_flux = 0
                        aperture_flux_error = 0
                        for pixel in flux_pixels:
                            overlap = pixel_to_aperture_overlap(pixel[1], pixel[2], x_mean, y_mean,
                                                                    self.targets_aperture[target] * self.psf_ratio[
                                                                        counter])
                            aperture_flux += (pixel[0] - sky) * overlap
                            aperture_flux_error = np.sqrt(
                                aperture_flux_error ** 2 + (np.abs(pixel[0]) + sky_error ** 2) * (overlap ** 2))

                        # save aperture photometry

                        self.results['aperture_x_position'][target][counter] = x_mean
                        self.results['aperture_y_position'][target][counter] = y_mean
                        self.results['aperture_flux'][target][counter] = aperture_flux
                        self.results['aperture_flux_error'][target][counter] = aperture_flux_error
                        self.results['aperture_sky'][target][counter] = sky
                        self.results['aperture_sky_error'][target][counter] = sky_error

                    except Exception as e:
                        print(e)
                        self.target_fails[target] += 1

                else:
                    self.target_fails[target] += 1

                if self.target_fails[target] > 0:
                    if target == 0:
                        name = 'Target'
                    else:
                        name = 'Comparison {0}'.format(target)
                    self.progress_diagnostics[target].set(
                        '{0}\nStar not found\nin {1} frame(s)'.format(name, int(self.target_fails[target])))

            for target in range(self.targets):

                if self.exit:
                    break

                if target == 0:
                    gauss_lc = self.results['gauss_flux'][0] / np.sum(self.results['gauss_flux'][1:], 0)
                    aperture_lc = self.results['aperture_flux'][0] / np.sum(self.results['aperture_flux'][1:], 0)

                else:
                    gauss_lc = self.results['gauss_flux'][target] / (
                            np.sum(self.results['gauss_flux'][1:target], 0) +
                            np.sum(self.results['gauss_flux'][target:], 0))
                    aperture_lc = self.results['aperture_flux'][target] / (
                            np.sum(self.results['aperture_flux'][1:target], 0) +
                            np.sum(self.results['aperture_flux'][target:], 0))

                gauss_lc = gauss_lc / np.nanmedian(gauss_lc)
                aperture_lc = aperture_lc / np.nanmedian(aperture_lc)

                min_data_point = np.nanmin([gauss_lc, aperture_lc])
                max_data_point = np.nanmax([gauss_lc, aperture_lc])

                if min_data_point < self.axes_lc_ylims[0]:
                    self.axes_lc_ylims[0] = self.axes_lc_ylims[0] - (
                                int((self.axes_lc_ylims[0] - min_data_point) / 0.01) + 1) * 0.01
                if max_data_point > self.axes_lc_ylims[1]:
                    self.axes_lc_ylims[1] = self.axes_lc_ylims[1] + (
                                int((max_data_point - self.axes_lc_ylims[1]) / 0.01) + 1) * 0.01
                if target == 0 or self.targets > 2:
                    self.images_gauss_lc[target][0].set_ydata(gauss_lc)
                    self.images_aperture_lc[target][0].set_ydata(aperture_lc)

            for ax_lc in self.axes_lc:
                ax_lc.set_ylim(self.axes_lc_ylims[0], self.axes_lc_ylims[1])

            self.progress_figure.draw()
            self.progress_photometry.update()
            self.science_counter += 1

            if self.science_counter >= len(self.science_files):
                self.after(self.options)
            else:
                self.after(self.measure)

    def options(self):

        if self.exit:
            self.def_close()
        else:
            self.save()

            if self.targets > 2:
                for i in self.progress_active[1:]:
                    i.activate()

    def replot_results(self):

        active = [True]

        for target in range(1, self.targets):

            if self.progress_active[target].get():
                active.append(True)
                self.images_text_lc[target].set_text(' ')
            else:
                active.append(False)
                self.images_text_lc[target].set_text('Inactive'.format(target))

        active = np.where(active)

        if len(active[0]) == 2:
            for target in range(1, self.targets):
                if self.progress_active[target].get():
                    self.progress_active[target].disable()
                    self.images_text_lc[target].set_text('Can\'t plot relative curve with only 1 active comparison')
        else:
            for target in range(1, self.targets):
                self.progress_active[target].activate()

        self.axes_lc_ylims = [0.99, 1.01]

        for target in range(self.targets):

            if target == 0:
                gauss_lc = self.results['gauss_flux'][0] / np.sum(self.results['gauss_flux'][active][1:], 0)
                aperture_lc = self.results['aperture_flux'][0] / np.sum(self.results['aperture_flux'][active][1:], 0)

            else:

                if target in active[0]:

                    gauss_comps = self.results['gauss_flux'][target] * 0
                    aperture_comps = self.results['gauss_flux'][target] * 0

                    for comp in range(1, self.targets):
                        if comp in active[0] and comp != target:
                            gauss_comps += self.results['gauss_flux'][comp]
                            aperture_comps += self.results['aperture_flux'][comp]

                    gauss_lc = self.results['gauss_flux'][target] / gauss_comps
                    aperture_lc = self.results['aperture_flux'][target] / aperture_comps

                else:
                    gauss_lc = self.results['gauss_flux'][target] * np.nan
                    aperture_lc = self.results['aperture_flux'][target] * np.nan

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                gauss_lc = gauss_lc / np.nanmedian(gauss_lc)
                aperture_lc = aperture_lc / np.nanmedian(aperture_lc)

                min_data_point = np.nanmin([gauss_lc, aperture_lc])
                max_data_point = np.nanmax([gauss_lc, aperture_lc])

                if min_data_point < self.axes_lc_ylims[0]:
                    self.axes_lc_ylims[0] = self.axes_lc_ylims[0] - (
                            int((self.axes_lc_ylims[0] - min_data_point) / 0.01) + 1) * 0.01
                if max_data_point > self.axes_lc_ylims[1]:
                    self.axes_lc_ylims[1] = self.axes_lc_ylims[1] + (
                            int((max_data_point - self.axes_lc_ylims[1]) / 0.01) + 1) * 0.01

            self.images_gauss_lc[target][0].set_ydata(gauss_lc)
            self.images_aperture_lc[target][0].set_ydata(aperture_lc)

        for ax_lc in self.axes_lc:
            ax_lc.set_ylim(self.axes_lc_ylims[0], self.axes_lc_ylims[1])

        self.progress_figure.draw()

    def save(self):

        self.stop_and_return_button.disable()
        self.set_close_button_function(self.no_action)
        self.save_button.disable()
        self.return_to_photomertry_button.disable()
        self.proceed_button.disable()

        photometry_directory = self.log.photometry_directory_base + '_1'

        if not os.path.isdir(photometry_directory):
            os.mkdir(photometry_directory)
        else:
            fi = 2
            while os.path.isdir('{0}_{1}'.format(self.log.photometry_directory_base, str(fi))):
                fi += 1
            photometry_directory = '{0}_{1}'.format(self.log.photometry_directory_base, str(fi))
            os.mkdir(photometry_directory)

        # save fov

        visible_fov_x_min = self.log.get_param('min_x')
        visible_fov_y_min = self.log.get_param('min_y')
        visible_fov_x_max = self.log.get_param('max_x')
        visible_fov_y_max = self.log.get_param('max_y')

        fov = self.FitsWindow(input=self.first_fits[1], input_options=self.log.get_param('photometry_fov_options'),
                              subplots_adjust=(0.01, 0.99, 0.05, 0.89))
        options = self.log.get_param('photometry_fov_options')

        fov.ax.add_patch(mpatches.Rectangle((visible_fov_x_min + 1, visible_fov_y_min + 1),
                                             visible_fov_x_max - visible_fov_x_min - 2,
                                             visible_fov_y_max - visible_fov_y_min - 2,
                                             ec='r', fill=False))

        text_y_drift = [-1,  1][int(options[3])]
        text_x_drift = [1,  -1][int(options[4])]
        for target in range(self.targets):

            if target == 0:
                fov.ax.text(self.targets_x_position[target] + text_x_drift * self.targets_aperture[target],
                            self.targets_y_position[target] + text_y_drift * self.targets_aperture[target],
                            '$T$({0})'.format(self.targets_aperture[target]), color='r', fontsize=8, va='top', ha='left')
                circle = mpatches.Circle((self.targets_x_position[target], self.targets_y_position[target]),
                                         self.targets_aperture[target],
                                         ec='r', fill=False, label='Target (ap. radius)')
                fov.ax.add_patch(circle)

            else:

                if target == 2:
                    label = 'Comparisons (ap. radius)'
                else:
                    label = None

                if self.progress_active[target].get():
                    fov.ax.text(self.targets_x_position[target] + text_x_drift * self.targets_aperture[target],
                                self.targets_y_position[target] + text_y_drift * self.targets_aperture[target],
                                '$C_{0}$({1})'.format(target, self.targets_aperture[target]), color='#07fefc', fontsize=8, va='top', ha='left')
                    circle = mpatches.Circle((self.targets_x_position[target], self.targets_y_position[target]),
                                             self.targets_aperture[target],
                                             ec='#07fefc', fill=False, label=label)
                    fov.ax.add_patch(circle)

                else:
                    fov.ax.text(self.targets_x_position[target] + text_x_drift * self.targets_aperture[target],
                                self.targets_y_position[target] + text_y_drift * self.targets_aperture[target],
                                '$C_{0}$({1}) - Inactive'.format(target, self.targets_aperture[target]), color='#07fefc', fontsize=8, va='top', ha='left')
                    circle = mpatches.Circle((self.targets_x_position[target], self.targets_y_position[target]),
                                             self.targets_aperture[target],
                                             ec='#07fefc', fill=False)
                    fov.ax.add_patch(circle)

        fov.ax.legend(loc=(0, 1.01))

        fov.figure.savefig(os.path.join(photometry_directory, self.log.fov_figure), dpi=200)

        self.progress_figure.figure.savefig(os.path.join(photometry_directory, self.log.lc_figure), dpi=200)

        # save log

        for target in range(1, self.targets):
            self.log.set_param('comparison_{0}_active'.format(target), self.progress_active[target].get())

        # save lightcurves

        active = []

        for target in range(self.targets):
            if self.progress_active[target].get():
                active.append(True)
            else:
                active.append(False)

        active = np.where(active)

        gauss_lc = self.results['gauss_flux'][0] / np.sum(self.results['gauss_flux'][active][1:], 0)
        gauss_lc_error = np.sqrt(
            (self.results['gauss_flux_error'][0] / self.results['gauss_flux'][0]) ** 2 +
            (np.sqrt(np.sum(self.results['gauss_flux_error'][active][1:] ** 2, 0)) / np.sum(
                self.results['gauss_flux'][active][1:], 0)) ** 2
        ) * gauss_lc
        aperture_lc = self.results['aperture_flux'][0] / np.sum(self.results['aperture_flux'][active][1:], 0)
        aperture_lc_error = np.sqrt(
            (self.results['aperture_flux_error'][0] / self.results['aperture_flux'][0]) ** 2 +
            (np.sqrt(np.sum(self.results['aperture_flux_error'][active][1:] ** 2, 0)) / np.sum(
                self.results['aperture_flux'][active][1:], 0)) ** 2
        ) * aperture_lc

        valid_gauss = np.where(~np.isnan(gauss_lc))
        valid_aperture = np.where(~np.isnan(aperture_lc))

        targets_results = [self.results['file_names'][valid_gauss], self.results['jd'][valid_gauss]]

        for array in ['gauss_x_position', 'gauss_y_position', 'gauss_x_std', 'gauss_y_std', 'gauss_flux',
                      'gauss_flux_error', 'gauss_sky', 'gauss_sky_error']:
            for target in range(self.targets):
                targets_results.append(self.results[array][target][valid_gauss])

        np.savetxt(os.path.join(photometry_directory, self.log.photometry_file.replace('.txt', '_g.txt')),
                   np.swapaxes(targets_results, 0, 1), fmt="%s")

        np.savetxt(os.path.join(photometry_directory, self.log.light_curve_gauss_file), np.swapaxes([
            self.results['jd'][valid_gauss], gauss_lc[valid_gauss], gauss_lc_error[valid_gauss]
        ], 0, 1))

        targets_results = [self.results['file_names'][valid_gauss], self.results['jd'][valid_aperture]]

        for array in ['aperture_x_position', 'aperture_y_position', 'aperture_flux',
                      'aperture_flux_error', 'aperture_sky', 'aperture_sky_error']:
            for target in range(self.targets):
                targets_results.append(self.results[array][target][valid_aperture])

        np.savetxt(os.path.join(photometry_directory, self.log.photometry_file.replace('.txt', '_a.txt')),
                   np.swapaxes(targets_results, 0, 1), fmt="%s")

        np.savetxt(os.path.join(photometry_directory, self.log.light_curve_aperture_file), np.swapaxes([
            self.results['jd'][valid_aperture], aperture_lc[valid_aperture], aperture_lc_error[valid_aperture]
        ], 0, 1))

        # save output description

        shutil.copy(self.log.files['photometry_output_description'], photometry_directory)

        # save exoclock info

        ra_dec_string = self.log.get_param('target_ra_dec')
        ra_dec_string = ra_dec_string.split(' ')
        try:
            planet = plc.locate_planet(plc.Hours(ra_dec_string[0]), plc.Degrees(ra_dec_string[1])).name
        except:
            planet = 'Could not find a planet in the ExoClock catalogue at this location'

        phot_filter = self.log.get_param('filter')

        if np.std(np.abs(gauss_lc[1:] - gauss_lc[:-1])) < np.std(np.abs(aperture_lc[1:] - aperture_lc[:-1])):
            files_to_upload = ['PHOTOMETRY_GAUSS.txt', 'PHOTOMETRY_APERTURE.txt']
        else:
            files_to_upload = ['PHOTOMETRY_APERTURE.txt', 'PHOTOMETRY_GAUSS.txt']

        w = open(os.path.join(photometry_directory, 'ExoClock_info.txt'), 'w')
        w.write('\n'.join([
            'The ExoClock Project is an effort to keep the ephemerides of exoplanets as precise as \n'
            'possible for planning future observations. If you have observed an exoplanet you can\n'
            'contribute your observation at: \n\nhttps://www.exoclock.space\n\n'
            'File to upload: {0} \n(this is a suggestion based on the scatter \nof your light curves, '
            'you can also try \nuploading {1})'.format(*files_to_upload),
            '',
            'Planet: {0}'.format(planet),
            '',
            'Time format: JD_UTC \n(UTC-based Julian date)',
            '',
            'Time stamp: Exposure start \n(the time produced refers to the beginning of each exposure)',
            '',
            'Flux format: Flux \n(flux of target over summed flux of comparisons)',
            '',
            'Filter: {0}'.format(phot_filter),
            '',
            'Exposure time in seconds: {0}'.format(self.exposure_time),
        ]))
        w.close()

        # saving ok

        self.log.set_param('photometry_complete', True)
        self.log.set_param('photometry_version', self.log.version)

        self.log.save_local_log()

        self.log.export_local_log(photometry_directory)

        self.showinfo('Results saved successfully',
                      'Your results have been saved successfully in {0}'.format(photometry_directory))
        self.show()

        self.set_close_button_function(self.close)
        self.save_button.activate()
        self.return_to_photomertry_button.activate()
        self.proceed_button.activate()

    def proceed(self):

        self.log.set_param('proceed', True)
        self.close()
