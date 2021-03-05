
import os
import time
import numpy as np
import matplotlib.patches as mpatches
import hops.pylightcurve3 as plc
from astropy.io import fits as pf

from hops.application_windows import MainWindow


class AlignmentWindow(MainWindow):

    def __init__(self, log):

        MainWindow.__init__(self, log, name='HOPS - Alignment', position=2)

        # set variables, create and place widgets

        self.bin_fits = self.log.get_param('bin_fits')
        self.burn_limit = int(1.1 * self.log.get_param('burn_limit')) * self.bin_fits * self.bin_fits
        self.shift_tolerance_p = self.log.get_param('shift_tolerance_p')
        self.rotation_tolerance = self.log.get_param('rotation_tolerance')
        self.min_calibration_stars_number = int(self.log.get_param('min_calibration_stars_number'))

        self.all_frames = plc.open_dict(self.log.all_frames)
        self.science_files = []
        for science_file in self.all_frames:
            if not self.all_frames[science_file][self.log.skip_key]:
                self.science_files.append([self.all_frames[science_file][self.log.time_key], science_file])
            else:
                self.all_frames[science_file][self.log.align_x0_key] = False
                self.all_frames[science_file][self.log.align_y0_key] = False
                self.all_frames[science_file][self.log.align_u0_key] = False

                fits = pf.open(os.path.join(self.log.reduction_directory, science_file), mode='update')
                fits[1].header.set(self.log.align_x0_key, False)
                fits[1].header.set(self.log.align_y0_key, False)
                fits[1].header.set(self.log.align_u0_key, False)
                fits.flush()
                fits.close()

        self.science_files.sort()
        self.science_files = [ff[1] for ff in self.science_files]

        self.skip_time = 0
        self.science_counter = 0
        self.test_level = None
        self.redraw = None
        self.stars = None
        self.science_file = None
        self.fits = None
        self.std = None
        self.mean = None
        self.star_std = None
        self.int_psf = None
        self.stars_detected = None
        self.rotation_detected = None
        self.check_num = None
        self.check_num_snr = None
        self.x0 = None
        self.y0 = None
        self.u0 = None
        self.f0 = None
        self.comparisons = None
        self.comparisons_snr = None
        self.small_angles = None
        self.large_angles = None
        self.circle = None
        self.settings_to_check = None
        self.comparisons_to_check = None

        # common definitions for all images

        fits = plc.open_fits(os.path.join(self.log.reduction_directory, self.science_files[0]))

        self.shift_tolerance = int(max(len(fits[1].data), len(fits[1].data[0])) * (self.shift_tolerance_p / 100.0))
        self.y_length, self.x_length = fits[1].data.shape
        self.circles_diameter = 0.02 * max(self.y_length, self.x_length)

        # progress window

        y_scale = (self.root.winfo_screenheight() - 500) / self.root.winfo_screenheight()

        self.progress_figure = self.FitsWindow(figsize=(0.5, y_scale, 10, 10, len(fits[1].data[0]) / len(fits[1].data)))
        self.progress_figure.load_fits(fits[1], input_name=self.science_files[0])
        self.progress_all_stars = self.Label(text='')
        self.progress_alignment = self.Progressbar(task="Aligning frames")
        # self.progress_all_frames = self.Progressbar(task="Aligning all stars in all frames")

        self.setup_window([
            [[self.progress_figure, 0, 2]],
            [[self.progress_all_stars, 0, 2]],
            [[self.progress_alignment, 0, 2]],
            [[self.Button(text='STOP ALIGNMENT & RETURN TO MAIN MENU', command=self.trigger_exit), 0, 2]],
            []
        ])

        self.set_close_button_function(self.trigger_exit)

    def run_alignment(self):

        self.close = self.trigger_exit

        if self.log.get_param('alignment_complete'):
            if self.askyesno('Overwrite files', 'Alignment has been completed, do you want to run again?'):
                self.log.set_param('alignment_complete', False)
                self.log.save_local_log()
            else:
                self.log.set_param('proceed', True)
            self.show()

        if not self.log.get_param('alignment_complete'):

            self.progress_all_stars.set('Analysing first frame...')

            self.after(self.find_all_stars)

        else:
            self.def_close()

    def alignment_log(self, *text):
        # print(*text)
        pass

    def find_all_stars(self):

        if self.exit:
            self.after(self.align)

        else:

            fits = plc.open_fits(os.path.join(self.log.reduction_directory, self.science_files[0]))
            metadata = self.all_frames[self.science_files[0]]

            self.progress_figure.load_fits(fits[1], self.science_files[0])

            stars, psf = plc.find_all_stars(fits[1].data,
                                            mean=metadata[self.log.mean_key], std=metadata[self.log.std_key],
                                            std_limit=3, burn_limit=self.burn_limit, star_std=metadata[self.log.psf_key],
                                            progressbar=self.progress_all_stars, progress_window=self,
                                            verbose=True
                                            )

            if self.exit:
                self.after(self.choose_calibartion_stars)

            stars = np.array(stars)

            self.log.save_local_log()

            all_stars_dict = {'all_stars': stars}
            plc.save_dict(all_stars_dict, 'all_stars.pickle')

            self.progress_all_stars.set('Choosing calibrating stars...')

            self.after(self.choose_calibartion_stars)

    def choose_calibartion_stars(self):

        if self.exit:
            self.after(self.align)

        else:

            all_stars_dict = plc.open_dict('all_stars.pickle')
            stars = np.array(all_stars_dict['all_stars'])

            fits = pf.open(os.path.join(self.log.reduction_directory, self.science_files[0]), memmap=False)
            metadata = self.all_frames[self.science_files[0]]

            frame_mean = metadata[self.log.mean_key]
            frame_std = metadata[self.log.std_key]
            frame_star_psf = metadata[self.log.psf_key]

            bright_stars = []
            std_limit = 30
            while len(bright_stars) < 100 and std_limit >= 5.0:
                bright_stars = []
                for star in stars:
                    if star[2] + star[3] < 2.0 * self.burn_limit / 3.0:
                        if star[-1] > (2 * np.pi * (std_limit * frame_std) * frame_star_psf * frame_star_psf):
                            self.alignment_log(star[0], star[1], star[-1],
                                               2 * np.pi * (frame_mean + std_limit * frame_std) * (frame_star_psf ** 2))
                            bright_stars.append(star)
                std_limit -= 5

            if len(bright_stars) < self.min_calibration_stars_number:
                bright_stars = []
                std_limit = 30
                while len(bright_stars) < 100 and std_limit >= 5.0:
                    bright_stars = []
                    for star in stars:
                        if star[-1] > (2 * np.pi * (std_limit * frame_std) * frame_star_psf * frame_star_psf):
                            self.alignment_log(star[0], star[1], star[-1],
                                               2 * np.pi * (frame_mean + std_limit * frame_std) * (frame_star_psf ** 2))
                            bright_stars.append(star)
                    std_limit -= 5

            stars = sorted(bright_stars, key=lambda x: -x[-1] / (x[-2] ** 3))

            x_ref_position = stars[0][0]
            y_ref_position = stars[0][1]
            f_ref = stars[0][-1]

            del stars[0]

            # take the rest as calibration stars and calculate their polar coordinates relatively to the first
            calibration_stars_polar = []
            for star in stars:
                r_position, u_position = plc.cartesian_to_polar(star[0], star[1], x_ref_position, y_ref_position)
                if r_position > 5 * frame_star_psf:
                    calibration_stars_polar.append([r_position, u_position])

            stars = sorted(stars, key=lambda x: -x[-1])

            calibration_stars_polar_snr = []
            for star in stars:
                r_position, u_position = plc.cartesian_to_polar(star[0], star[1], x_ref_position, y_ref_position)
                if r_position > 5 * frame_star_psf:
                    calibration_stars_polar_snr.append([r_position, u_position])

            if len(calibration_stars_polar) <= self.min_calibration_stars_number:
                self.check_num = len(calibration_stars_polar) - 0.5
                self.check_num_snr = len(calibration_stars_polar) - 0.5
            else:
                self.check_num = max(self.min_calibration_stars_number - 0.5, len(calibration_stars_polar) / 10.0 - 0.5)
                self.check_num_snr = max(self.min_calibration_stars_number - 0.5, len(calibration_stars_polar) / 20.0 - 0.5)

            self.x0 = x_ref_position
            self.y0 = y_ref_position
            self.u0 = 0
            self.f0 = f_ref
            self.comparisons = calibration_stars_polar
            self.comparisons_snr = calibration_stars_polar_snr
            fits.close()

            ustep = np.arcsin(float(frame_star_psf) / self.comparisons[int(len(self.comparisons) / 2)][0])
            self.small_angles = np.append(np.arange(-self.rotation_tolerance, self.rotation_tolerance, ustep),
                                          np.arange(-self.rotation_tolerance, self.rotation_tolerance, ustep) + np.pi)

            self.large_angles = np.array([np.pi, 0])
            for ff in range(1, int(np.pi / ustep) + 1):
                self.large_angles = np.append(self.large_angles, np.pi - ff * ustep)
                self.large_angles = np.append(self.large_angles, np.pi + ff * ustep)
                self.large_angles = np.append(self.large_angles, 0 - ff * ustep)
                self.large_angles = np.append(self.large_angles, 0 + ff * ustep)

            # set the looking window and angular step

            self.progress_all_stars.set(' ')

            self.after(self.align)

    def align(self):

        if self.exit:
            self.after(self.plot_current)

        else:

            if self.science_counter == 0:
                self.progress_alignment.initiate(len(self.science_files))

            self.stars = None
            self.science_file = self.science_files[self.science_counter]
            self.fits = pf.open(os.path.join(self.log.reduction_directory, self.science_file), memmap=False, mode='update')
            self.std = self.fits[1].header[self.log.std_key]
            self.mean = self.fits[1].header[self.log.mean_key]
            self.star_std = self.fits[1].header[self.log.psf_key]
            self.int_psf = int(max(1, round(self.fits[1].header[self.log.psf_key])))
            self.stars_detected = False
            self.rotation_detected = False

            self.progress_alignment.show_message(' ')

            self.progress_figure.load_fits(self.fits[1], self.science_file, draw=False)
            self.circle = mpatches.Circle((self.x0, self.y0), self.circles_diameter, ec='r', fill=False)
            self.progress_figure.ax.add_patch(self.circle)

            self.test_level = 1
            self.redraw = 0
            self.skip_time = 0

            self.after(self.detect_stars)

    def detect_stars(self):

        if self.exit:
            self.after(self.plot_current)

        else:

            if self.test_level == 1:
                self.stars = plc.find_single_star(self.fits[1].data, self.x0, self.y0, mean=self.mean, std=self.std,
                                                  burn_limit=self.burn_limit, star_std=self.star_std)
                if self.stars:
                    self.stars.append(2 * np.pi * self.stars[2] * self.stars[4] * self.stars[5])
                    self.stars = [self.stars]

            elif self.test_level == 2:
                self.skip_time = time.time()
                self.stars = plc.find_all_stars(self.fits[1].data, x_low=self.x0 - self.shift_tolerance,
                                                x_upper=self.x0 + self.shift_tolerance,
                                                y_low=self.y0 - self.shift_tolerance,
                                                y_upper=self.y0 + self.shift_tolerance, x_centre=self.x0,
                                                y_centre=self.y0, mean=self.mean,
                                                std=self.std, burn_limit=self.burn_limit,
                                                star_std=self.star_std, verbose=True,
                                                order_by_distance_and_flux=self.f0)[0]
                self.progress_all_stars.set(' ')

            elif self.test_level == 4:
                if self.askyesno('HOPS - Alignment', 'Stars not found close to their previous positions.\n'
                                                     'Do you want to skip this frame?'):
                    self.after(self.plot_current)
                self.stars = plc.find_all_stars(self.fits[1].data, mean=self.mean, std=self.std, burn_limit=self.burn_limit,
                                                star_std=self.star_std, order_by_flux=self.f0, verbose=True)[0]
                self.progress_all_stars.set(' ')

            if self.stars:

                self.settings_to_check = []

                if self.test_level == 1:
                    self.settings_to_check.append([self.stars[0][0], self.stars[0][1], self.u0, self.stars[0]])
                    self.setting_checking = 0
                    self.comparisons_to_check = self.comparisons
                elif self.test_level == 2:
                    for star in self.stars:
                        self.settings_to_check.append([star[0], star[1], self.u0, star])
                    self.comparisons_to_check = self.comparisons_snr
                elif self.test_level == 3:
                    for star in self.stars:
                        for rotation in self.small_angles:
                            self.settings_to_check.append([star[0], star[1], rotation, star])
                elif self.test_level == 4:
                    for star in self.stars:
                        self.settings_to_check.append([star[0], star[1], self.u0, star])
                elif self.test_level == 5:
                    for star in self.stars:
                        for rotation in self.large_angles:
                            self.settings_to_check.append([star[0], star[1], rotation, star])

                self.setting_checking = 0
                self.after(self.check_star)
            else:
                if self.test_level == 1:
                    self.test_level = 2
                    self.progress_figure.draw()
                    self.progress_all_stars.set('Analysing frame...')
                    self.progress_alignment.show_message('Testing small shift...')
                    self.after(self.detect_stars)
                elif self.test_level == 2:
                    self.test_level = 3
                    self.progress_alignment.show_message('Testing small shift & rotation...')
                    self.after(self.detect_stars)
                elif self.test_level == 3:
                    self.test_level = 4
                    self.progress_all_stars.set('Analysing frame...')
                    self.progress_alignment.show_message('Testing large shift...')
                    self.after(self.detect_stars)
                elif self.test_level == 4:
                    self.test_level = 5
                    self.progress_alignment.show_message('Testing large shift & rotation...')
                    self.after(self.detect_stars)
                else:
                    self.after(self.plot_current)

    def check_star(self):

        if self.exit:
            self.after(self.plot_current)

        else:

            x, y, u, star = self.settings_to_check[self.setting_checking]
            self.alignment_log('Checking star at: ', x, y, ', with rotation:', u)

            if self.redraw >= 1:
                self.circle.set_center((x, y))
                self.progress_figure.draw()
                self.redraw = 0
            else:
                self.redraw += 0.01

            test = 0

            for comp in self.comparisons_to_check:

                check_x = int(x + comp[0] * np.cos(u + comp[1]))
                check_y = int(y + comp[0] * np.sin(u + comp[1]))
                if 0 < check_x < self.x_length and 0 < check_y < self.y_length:
                    check_sum = np.sum(self.fits[1].data[check_y - self.int_psf:check_y + self.int_psf + 1,
                                       check_x - self.int_psf:check_x + self.int_psf + 1])
                    check_lim = (self.fits[1].header[self.log.mean_key] +
                                 3 * self.fits[1].header[self.log.std_key]) * ((2 * self.int_psf + 1) ** 2)
                    if check_sum > check_lim:
                        test += 1
                    else:
                        test -= 1
                self.alignment_log('Check ref. star at: ', check_x, check_y, ', Test: ', test)

                if abs(test) > self.check_num:
                    break

            if test >= self.check_num:
                self.stars_detected = True
                if self.test_level > 1:
                    self.rotation_detected = True
                self.x0 = x
                self.y0 = y
                self.u0 = u
                self.f0 = star[-1]
                self.after(self.plot_current)
            else:
                self.setting_checking += 1
                if self.setting_checking < len(self.settings_to_check):
                    self.after(self.check_star)
                else:
                    if self.test_level == 1:
                        self.test_level = 2
                        self.progress_figure.draw()
                        self.progress_all_stars.set('Analysing frame...')
                        self.progress_alignment.show_message('Testing small shift...')
                        self.after(self.detect_stars)
                    elif self.test_level == 2:
                        self.test_level = 3
                        self.progress_alignment.show_message('Testing small shift & rotation...')
                        self.after(self.detect_stars)
                    elif self.test_level == 3:
                        self.test_level = 4
                        self.progress_all_stars.set('Analysing frame...')
                        self.progress_alignment.show_message('Testing large shift...')
                        self.after(self.detect_stars)
                    elif self.test_level == 4:
                        self.test_level = 5
                        self.progress_alignment.show_message('Testing large shift & rotation...')
                        self.after(self.detect_stars)
                    else:
                        self.after(self.plot_current)

    def plot_current(self):

        if self.exit:
            self.after(self.save)

        else:

            if self.stars_detected:

                if self.rotation_detected:

                    test_u0 = []
                    test_cos = []
                    test_sin = []

                    for ii in self.comparisons[:int(self.check_num + 0.5)]:
                        check_x = self.x0 + ii[0] * np.cos(self.u0 + ii[1])
                        check_y = self.y0 + ii[0] * np.sin(self.u0 + ii[1])
                        star = plc.find_single_star(self.fits[1].data, check_x, check_y, mean=self.mean,  std=self.std,
                                                    burn_limit=self.burn_limit, star_std=self.star_std)

                        if star:
                            diff = plc.cartesian_to_polar(star[0], star[1], self.x0, self.y0)[1] - ii[1]
                            if diff < 0:
                                diff += 2 * np.pi
                            test_u0.append(diff)
                            test_cos.append(np.cos(diff))
                            test_sin.append(np.sin(diff))

                    if len(test_u0) > 0:
                        test_cos = np.median(test_cos)
                        test_sin = np.median(test_sin)

                        self.u0 = np.arccos(test_cos)
                        if test_sin < 0:
                            self.u0 = np.pi + (np.pi - self.u0)

                self.fits[1].header.set(self.log.align_x0_key, self.x0)
                self.fits[1].header.set(self.log.align_y0_key, self.y0)
                self.fits[1].header.set(self.log.align_u0_key, self.u0)

                self.circle.set_center((self.x0, self.y0))
                for ii in self.comparisons[:int(self.check_num + 0.5)]:
                    circle = mpatches.Circle((self.x0 + ii[0] * np.cos(self.u0 + ii[1]),
                                              self.y0 + ii[0] * np.sin(self.u0 + ii[1])),
                                             self.circles_diameter, ec='w', fill=False)
                    self.progress_figure.ax.add_patch(circle)

            else:

                self.fits[1].header.set(self.log.align_x0_key, False)
                self.fits[1].header.set(self.log.align_y0_key, False)
                self.fits[1].header.set(self.log.align_u0_key, False)

            self.all_frames[self.science_file][self.log.align_x0_key] = self.fits[1].header[self.log.align_x0_key]
            self.all_frames[self.science_file][self.log.align_y0_key] = self.fits[1].header[self.log.align_y0_key]
            self.all_frames[self.science_file][self.log.align_u0_key] = self.fits[1].header[self.log.align_u0_key]

            if not self.fits[1].header[self.log.align_x0_key]:
                self.all_frames[self.science_file][self.log.skip_key] = True

            self.fits.flush()
            self.fits.close()

            self.progress_figure.draw()
            if self.skip_time == 0:
                self.progress_alignment.update()
            else:
                self.progress_alignment.update(skip=time.time() - self.skip_time)
                self.skip_time = 0
            self.science_counter += 1

            if self.science_counter >= len(self.science_files):
                self.progress_all_stars.set('Aligning all stars in all frames...')
                self.after(self.save)
            else:
                self.after(self.align)

    def save(self):

        if self.exit:
            self.after(self.check_visibility)
        else:
            plc.save_dict(self.all_frames, self.log.all_frames)
            self.after(self.check_visibility)

    def check_visibility(self):

        if self.exit:
            self.def_close()

        else:

            all_stars_dict = plc.open_dict('all_stars.pickle')
            stars = np.array(all_stars_dict['all_stars'])
            fits = plc.open_fits(os.path.join(self.log.reduction_directory, self.science_files[0]))

            polar_coords = []
            for star in all_stars_dict['all_stars']:
                polar_coords.append(plc.cartesian_to_polar(star[0], star[1],
                                                           self.all_frames[self.science_files[0]][
                                                               self.log.align_x0_key],
                                                           self.all_frames[self.science_files[0]][
                                                               self.log.align_y0_key]))

            in_fov = np.ones(len(polar_coords))

            for science_file in self.science_files:

                metadata = self.all_frames[science_file]

                if self.exit:
                    self.def_close()

                ref_x_position = metadata[self.log.align_x0_key]
                ref_y_position = metadata[self.log.align_y0_key]
                ref_u_position = metadata[self.log.align_u0_key]
                star_std = metadata[self.log.psf_key]

                if ref_x_position:

                    in_fov_single = []

                    for star in polar_coords:

                        cartesian_x = ref_x_position + star[0] * np.cos(ref_u_position + star[1])
                        cartesian_y = ref_y_position + star[0] * np.sin(ref_u_position + star[1])

                        if (3 * star_std < cartesian_x < len(fits[1].data[0]) - 3 * star_std
                                and 3 * star_std < cartesian_y < len(fits[1].data) - 3 * star_std):

                            in_fov_single.append(1)
                        else:
                            in_fov_single.append(0)

                    in_fov *= in_fov_single

            all_stars_dict['in_fov'] = np.array(in_fov)

            visible_fov_x_min = np.min(stars[np.where(in_fov), 0]) - 3 * star_std
            visible_fov_x_max = np.max(stars[np.where(in_fov), 0]) + 3 * star_std
            visible_fov_y_min = np.min(stars[np.where(in_fov), 1]) - 3 * star_std
            visible_fov_y_max = np.max(stars[np.where(in_fov), 1]) + 3 * star_std

            self.log.set_param('min_x', float(visible_fov_x_min))
            self.log.set_param('min_y', float(visible_fov_y_min))
            self.log.set_param('max_x', float(visible_fov_x_max))
            self.log.set_param('max_y', float(visible_fov_y_max))

            plc.save_dict(all_stars_dict, 'all_stars.pickle')

            self.log.set_param('alignment_complete', True)
            self.log.set_param('alignment_version', self.log.version)
            self.log.save_local_log()
            self.log.set_param('proceed', True)

            self.def_close()
