
from tkinter import Tk, TclError
from tkinter import Label, Button, Entry, Checkbutton, Scrollbar, Listbox, PhotoImage, Radiobutton, Scale, Frame
from tkinter import StringVar, BooleanVar, DoubleVar, IntVar
from tkinter import DISABLED, NORMAL, END, RIGHT, LEFT, BOTH, Y, HORIZONTAL

import tkinter.ttk as ttk
from tkinter.messagebox import *

import warnings
warnings.filterwarnings(
    'ignore', message='Matplotlib is building the font cache using fc-list. This may take a moment.')
warnings.filterwarnings(
    'ignore', message='The installed version of numexpr 2.4.4 is not supported in pandas and will be not be used')

import matplotlib
matplotlib.use('TkAgg')

import time
import numpy as np
import shutil
import hops.pylightcurve3 as plc
import matplotlib.cm as cm
import matplotlib.patches as mpatches

from astropy.io import fits as pf
from matplotlib.figure import Figure
try:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    NavigationToolbar2TkAgg = NavigationToolbar2Tk
except ImportError:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg


from hops.hops_tools.windows import *
from hops.hops_tools.logs import log

import glob

def find_fits_files(fits_file):

    fits_list = glob.glob('*{0}*.f*t*'.format(fits_file)) + glob.glob('*{0}*.F*T*'.format(fits_file))
    fits_list = list(np.unique(fits_list))
    fits_list.sort()
    return fits_list



def alignment_log(*message):
    # print(*message)
    pass


def alignment():

    print('Alignment...')

    if log.read_local_log('pipeline', 'alignment_complete'):
        if askyesno('Overwrite alignment', 'Alignment has been completed, do you want to run again?'):
            log.write_local_log('pipeline', False, 'alignment_complete')

    test_all_stars = True
    try:
        all_stars = plc.open_dict('all_stars.pickle')
        in_fov = all_stars['in_fov']
        all_stars = all_stars['all_stars']
    except:
        test_all_stars = False

    if log.read_local_log('alignment', 'star_psf_x') == 0 or not log.read_local_log('pipeline', 'alignment_complete'):
        test_all_stars = False

    # if all_stars exist and repeat not requested: exit

    if test_all_stars and log.read_local_log('pipeline', 'alignment_complete'):
        return 0

    # continue to alignment

    # load parameters

    reduction_directory = log.read_local_log('pipeline', 'reduction_directory')

    mean_key = log.read_local_log('pipeline_keywords', 'mean_key')
    std_key = log.read_local_log('pipeline_keywords', 'std_key')
    align_x0_key = log.read_local_log('pipeline_keywords', 'align_x0_key')
    align_y0_key = log.read_local_log('pipeline_keywords', 'align_y0_key')
    align_u0_key = log.read_local_log('pipeline_keywords', 'align_u0_key')

    frame_low_std = log.read_local_log('windows', 'frame_low_std')
    frame_upper_std = log.read_local_log('windows', 'frame_upper_std')

    bin_fits = int(log.read_local_log('reduction', 'bin_fits'))

    burn_limit = int(log.read_local_log('alignment', 'burn_limit')) * bin_fits * bin_fits
    shift_tolerance_p = log.read_local_log('alignment', 'shift_tolerance_p')
    rotation_tolerance = log.read_local_log('alignment', 'rotation_tolerance')
    min_calibration_stars_number = int(log.read_local_log('alignment', 'min_calibration_stars_number'))
    star_std = log.read_local_log('alignment', 'star_std')

    science = find_fits_files(os.path.join(reduction_directory, '*'))

    def find_all_stars():

        stars, psf = plc.find_all_stars(fits[1].data,
                                        mean=fits[1].header[mean_key], std=fits[1].header[std_key],
                                        std_limit=3, burn_limit=burn_limit, star_std=star_std,
                                        order_by_flux=False,
                                        progress_pack=(progress_bar, percent_label))

        stars = np.array(stars)

        log.write_local_log('alignment', int(max(1, round(max(psf[0], psf[1])))), 'star_std')
        log.write_local_log('alignment', float(max(psf[0], psf[1])), 'star_psf')
        log.write_local_log('alignment', float(psf[0]), 'star_psf_x')
        log.write_local_log('alignment', float(psf[1]), 'star_psf_y')

        all_stars_dict = {'all_stars': stars}
        plc.save_dict(all_stars_dict, 'all_stars.pickle')

    def align():

        star_std = log.read_local_log('alignment', 'star_std')
        psf = (log.read_local_log('alignment', 'star_psf_x'), log.read_local_log('alignment', 'star_psf_y'))

        all_stars_dict = plc.open_dict('all_stars.pickle')
        stars = np.array(all_stars_dict['all_stars'])

        fits = pf.open(science[0], memmap=False)
        frame_mean = fits[1].header[mean_key]
        frame_std = fits[1].header[std_key]
        shift_tolerance = int(max(len(fits[1].data), len(fits[1].data[0])) * (shift_tolerance_p / 100.0))
        y_length, x_length = fits[1].data.shape
        circles_diameter = 0.02 * max(y_length, x_length)

        bright_stars = []
        std_limit = 30
        while len(bright_stars) < 100 and std_limit >= 5.0:
            bright_stars = []
            for star in stars:
                if star[2] + star[3] < 2.0 * burn_limit / 3.0:
                    if star[-1] > (2 * np.pi * (std_limit * frame_std) * psf[0] * psf[1]):
                        alignment_log(star[-1], 2 * np.pi * (frame_mean + std_limit * frame_std) * psf[0] * psf[1])
                        bright_stars.append(star)
            std_limit -= 5

        stars = sorted(bright_stars, key=lambda x: -x[-1] / (x[-2]**3))

        x_ref_position = stars[0][0]
        y_ref_position = stars[0][1]

        del stars[0]

        # take the rest as calibration stars and calculate their polar coordinates relatively to the first
        calibration_stars_polar = []
        for star in stars:
            r_position, u_position = plc.cartesian_to_polar(star[0], star[1], x_ref_position, y_ref_position)
            if r_position > 5 * star_std:
                calibration_stars_polar.append([r_position, u_position])

        stars = sorted(stars, key=lambda x: -x[-1])

        calibration_stars_polar_snr = []
        for star in stars:
            r_position, u_position = plc.cartesian_to_polar(star[0], star[1], x_ref_position, y_ref_position)
            if r_position > 5 * star_std:
                calibration_stars_polar_snr.append([r_position, u_position])

        if len(calibration_stars_polar) <= min_calibration_stars_number:
            check_num = len(calibration_stars_polar) - 0.5
            check_num_snr = len(calibration_stars_polar) - 0.5
        else:
            check_num = max(min_calibration_stars_number - 0.5, int(len(calibration_stars_polar)) / 10 - 0.5)
            check_num_snr = max(min_calibration_stars_number - 0.5, int(len(calibration_stars_polar)) / 20 - 0.5)

        x0, y0, u0, comparisons = x_ref_position, y_ref_position, 0, calibration_stars_polar
        x0, y0, u0, comparisons_snr = x_ref_position, y_ref_position, 0, calibration_stars_polar_snr
        fits.close()

        # set the looking window and angular step

        circle = mpatches.Circle((x0, y0), circles_diameter, ec='r', fill=False)
        ax.add_patch(circle)
        for ii in comparisons[:10]:
            circle = mpatches.Circle((x0 + ii[0] * np.cos(u0 + ii[1]), y0 + ii[0] * np.sin(u0 + ii[1])),
                                     circles_diameter, ec='w', fill=False)
            ax.add_patch(circle)

        # for each science_file
        percent = 0
        skip_time = 0
        lt0 = time.time()
        for counter, science_file in enumerate(science):

            if show_progress.exit:
                return None

            alignment_log('\n', science_file)
            label_2.configure(text='Running Alignment: {0}'.format(science_file.split(os.sep)[-1]))
            label_2.update()
            label_3.configure(text='')
            label_3.update()

            fits = pf.open(science_file, mode='update')

            ax.cla()
            ax.imshow(fits[1].data, origin='lower', cmap=cm.Greys_r,
                      vmin=fits[1].header[mean_key] + frame_low_std * fits[1].header[std_key],
                      vmax=fits[1].header[mean_key] + frame_upper_std * fits[1].header[std_key])
            ax.axis('off')
            circle = mpatches.Circle((-100, -100), circles_diameter, ec='r', fill=False)
            ax.add_patch(circle)

            rotation_detected = False

            # super fast detection test
            alignment_log('Test no shift', x0, y0, u0)
            star = plc.find_single_star(fits[1].data, x0, y0, mean=fits[1].header[mean_key], std=fits[1].header[std_key],
                                        burn_limit=burn_limit, star_std=star_std)

            if star:

                tests = []

                max_x = star[0]
                max_y = star[1]
                alignment_log(science_file)
                alignment_log('Testing star at: ', max_x, max_y, ', with rotation:', u0)

                test = 0

                for comp in comparisons:

                    check_x = int(max_x + comp[0] * np.cos(u0 + comp[1]))
                    check_y = int(max_y + comp[0] * np.sin(u0 + comp[1]))
                    if 0 < check_x < x_length and 0 < check_y < y_length:
                        check_sum = np.sum(fits[1].data[check_y - star_std:check_y + star_std + 1,
                                           check_x - star_std:check_x + star_std + 1])
                        check_lim = (fits[1].header[mean_key] + 3 * fits[1].header[std_key]) * ((2 * star_std + 1) ** 2)
                        if check_sum > check_lim:
                            test += 1
                        else:
                            test -= 1
                    alignment_log('Check ref. star at: ', check_x, check_y, ', Test: ', test)

                    if abs(test) > check_num:
                        break

                tests.append([test, max_x, max_y])
                alignment_log([test, max_x, max_y])

                tests.sort()
                test, max_x, max_y = tests[-1]

                if test < check_num:

                    stars_detected = False

                else:
                    stars_detected = True
                    x0 = max_x
                    y0 = max_y
                    u0 = u0

            else:
                stars_detected = False

            alignment_log(stars_detected, x0, y0, u0)
            # super fast detection test

            if show_progress.exit:
                return None

            delta_skip_time = time.time()
            # look for reasonable field shift
            if not stars_detected:

                alignment_log('Testing small shift')
                label_3.configure(text='Testing small shift and rotation')
                label_3.update()

                stars = plc.find_all_stars(fits[1].data, x_low=x0 - shift_tolerance, x_upper=x0 + shift_tolerance,
                                           y_low=y0 - shift_tolerance, y_upper=y0 + shift_tolerance, x_centre=x0,
                                           y_centre=y0, mean=fits[1].header[mean_key], std=fits[1].header[std_key],
                                           burn_limit=burn_limit, star_std=star_std)[0]

                if show_progress.exit:
                    return None

                if stars:

                    tests = []

                    for star in stars:

                        if show_progress.exit:
                            return None

                        max_x = star[0]
                        max_y = star[1]
                        circle.set_center((max_x, max_y))
                        canvas.draw()
                        show_progress.update()

                        alignment_log(science_file)
                        alignment_log('Testing star at: ', max_x, max_y, ',with rotation:', u0)

                        test = 0

                        for comp in comparisons:

                            check_x = int(max_x + comp[0] * np.cos(u0 + comp[1]))
                            check_y = int(max_y + comp[0] * np.sin(u0 + comp[1]))
                            if 0 < check_x < x_length and 0 < check_y < y_length:
                                check_sum = np.sum(fits[1].data[check_y - star_std:check_y + star_std + 1,
                                                   check_x - star_std:check_x + star_std + 1])
                                check_lim = (fits[1].header[mean_key] + 3 * fits[1].header[std_key]) * (
                                            (2 * star_std + 1) ** 2)
                                if check_sum > check_lim:
                                    test += 1
                                else:
                                    test -= 1
                            alignment_log('Check ref. star at: ', check_x, check_y, ', Test: ', test)

                            if abs(test) > check_num:
                                break

                        tests.append([test, max_x, max_y])
                        alignment_log([test, max_x, max_y])

                        if test > check_num:
                            break

                    tests.sort()
                    test, max_x, max_y = tests[-1]

                    if test < check_num:

                        tests = []

                        for star in stars:

                            if show_progress.exit:
                                return None

                            max_x = star[0]
                            max_y = star[1]
                            circle.set_center((max_x, max_y))
                            canvas.draw()
                            show_progress.update()
                            alignment_log(science_file)
                            alignment_log('LOW-SNR Testing star at: ', max_x, max_y, ', with rotation:', u0)

                            test = 0

                            for comp in comparisons_snr:

                                check_x = int(max_x + comp[0] * np.cos(u0 + comp[1]))
                                check_y = int(max_y + comp[0] * np.sin(u0 + comp[1]))
                                if 0 < check_x < x_length and 0 < check_y < y_length:
                                    check_sum = np.sum(fits[1].data[check_y - star_std:check_y + star_std + 1,
                                                       check_x - star_std:check_x + star_std + 1])
                                    check_lim = (fits[1].header[mean_key] + 3 * fits[1].header[std_key]) * (
                                                (2 * star_std + 1) ** 2)
                                    if check_sum > check_lim:
                                        test += 1
                                    else:
                                        test -= 1
                                alignment_log('Check ref. star at: ', check_x, check_y, ', Test: ', test)

                                if abs(test) > check_num_snr:
                                    break

                            tests.append([test, max_x, max_y])
                            alignment_log([test, max_x, max_y])

                            if test > check_num_snr:
                                break

                            tests.sort()
                            test, max_x, max_y = tests[-1]

                            if test < min_calibration_stars_number - 0.5:
                                stars_detected = False

                            else:
                                stars_detected = True
                                x0 = max_x
                                y0 = max_y
                                u0 = u0
                                rotation_detected = True

                    else:
                        stars_detected = True
                        x0 = max_x
                        y0 = max_y
                        u0 = u0
                        rotation_detected = True

                else:
                    stars_detected = False

                alignment_log(stars_detected, x0, y0, u0)
            # look for reasonable field shift

            if show_progress.exit:
                return None

            # look for reasonable field rotation
            if not stars_detected:

                alignment_log('Testing small rotation')

                ustep = np.arcsin(float(star_std) / comparisons[int(len(comparisons) / 2)][0])
                angles = np.append(np.arange(-rotation_tolerance, rotation_tolerance, ustep),
                                   np.arange(-rotation_tolerance, rotation_tolerance, ustep) + np.pi)

                if stars:

                    tests = []

                    for star in stars:

                        if show_progress.exit:
                            return None

                        test = 0

                        max_x = star[0]
                        max_y = star[1]
                        circle.set_center((max_x, max_y))
                        canvas.draw()
                        show_progress.update()

                        for rotation in angles:

                            alignment_log(science_file)
                            alignment_log('Testing star at: ', max_x, max_y, ', with rotation: ', rotation)

                            test = 0

                            for comp in comparisons:

                                check_x = int(max_x + comp[0] * np.cos(rotation + comp[1]))
                                check_y = int(max_y + comp[0] * np.sin(rotation + comp[1]))
                                if 0 < check_x < x_length and 0 < check_y < y_length:
                                    check_sum = np.sum(fits[1].data[check_y - star_std:check_y + star_std + 1,
                                                       check_x - star_std:check_x + star_std + 1])
                                    check_lim = (fits[1].header[mean_key] + 3 * fits[1].header[std_key]) * (
                                            (2 * star_std + 1) ** 2)
                                    if check_sum > check_lim:
                                        test += 1
                                    else:
                                        test -= 1
                                alignment_log('Check ref. star at: ', check_x, check_y, ', Test: ', test)

                                if abs(test) > check_num:
                                    break

                            tests.append([test, max_x, max_y, rotation])
                            alignment_log([test, max_x, max_y, rotation])

                            if test > check_num:
                                break

                        if test > check_num:
                            break

                    tests.sort()
                    test, max_x, max_y, rotation = tests[-1]
                    if test < check_num:

                        for star in stars:

                            if show_progress.exit:
                                return None

                            test = 0

                            max_x = star[0]
                            max_y = star[1]
                            circle.set_center((max_x, max_y))
                            canvas.draw()
                            show_progress.update()

                            for rotation in angles:

                                alignment_log(science_file)
                                alignment_log('LOW SNR Testing star at: ', max_x, max_y, ', with rotation: ', rotation)

                                test = 0

                                for comp in comparisons_snr:

                                    check_x = int(max_x + comp[0] * np.cos(rotation + comp[1]))
                                    check_y = int(max_y + comp[0] * np.sin(rotation + comp[1]))
                                    if 0 < check_x < x_length and 0 < check_y < y_length:
                                        check_sum = np.sum(fits[1].data[check_y - star_std:check_y + star_std + 1,
                                                           check_x - star_std:check_x + star_std + 1])
                                        check_lim = (fits[1].header[mean_key] + 3 * fits[1].header[std_key]) * (
                                                (2 * star_std + 1) ** 2)
                                        if check_sum > check_lim:
                                            test += 1
                                        else:
                                            test -= 1
                                    alignment_log('Check ref. star at: ', check_x, check_y, ', Test: ', test)

                                    if abs(test) > check_num_snr:
                                        break

                                tests.append([test, max_x, max_y, rotation])
                                alignment_log([test, max_x, max_y, rotation])

                                if test > check_num_snr:
                                    break

                            if test > check_num_snr:
                                break

                        tests.sort()
                        test, max_x, max_y, rotation = tests[-1]
                        if test < check_num_snr:
                            stars_detected = False

                        else:
                            stars_detected = True
                            x0 = max_x
                            y0 = max_y
                            u0 = rotation
                            rotation_detected = True

                    else:
                        stars_detected = True
                        x0 = max_x
                        y0 = max_y
                        u0 = rotation
                        rotation_detected = True

                else:
                    stars_detected = False

                alignment_log(stars_detected, x0, y0, u0)
            # look for reasonable field rotation

            if not stars_detected:
                circle.set_center((-100, -100))
                canvas.draw()
                show_progress.update()
                skip_frame = askyesno('Alignment',
                                      'Stars not found close to their previous positions.\n'
                                      'Do you want to skip this frame?',
                                      parent=show_progress.root)
            else:
                skip_frame = False

            if show_progress.exit:
                return None

            # look for large field shift
            if not stars_detected and not skip_frame:

                alignment_log('Testing large shift and rotation')
                label_3.configure(text='Testing large shift and rotation')
                label_3.update()
                canvas.draw()
                show_progress.update()

                stars = plc.find_all_stars(fits[1].data, mean=fits[1].header[mean_key], std=fits[1].header[std_key],
                                           burn_limit=burn_limit, star_std=star_std, order_by_flux=True)[0]

                if show_progress.exit:
                    return None

                if stars:

                    tests = []

                    for star in stars:

                        if show_progress.exit:
                            return None

                        max_x = star[0]
                        max_y = star[1]
                        circle.set_center((max_x, max_y))
                        canvas.draw()
                        show_progress.update()
                        alignment_log(science_file)
                        alignment_log('LOW SNR Testing star at: ', max_x, max_y, ', with rotation: ', u0)

                        test = 0

                        for comp in comparisons_snr:

                            check_x = int(max_x + comp[0] * np.cos(u0 + comp[1]))
                            check_y = int(max_y + comp[0] * np.sin(u0 + comp[1]))
                            if 0 < check_x < x_length and 0 < check_y < y_length:
                                check_sum = np.sum(fits[1].data[check_y - star_std:check_y + star_std + 1,
                                                   check_x - star_std:check_x + star_std + 1])
                                check_lim = (fits[1].header[mean_key] + 3 * fits[1].header[std_key]) * (
                                        (2 * star_std + 1) ** 2)
                                if check_sum > check_lim:
                                    test += 1
                                else:
                                    test -= 1
                                alignment_log('Check ref. star at: ', check_x, check_y, ', Test: ', test)

                            if abs(test) > check_num_snr:
                                break

                        tests.append([test, max_x, max_y])
                        alignment_log([test, max_x, max_y])

                        if test > check_num_snr:
                            break

                    tests.sort()
                    test, max_x, max_y = tests[-1]
                    if test < check_num_snr:
                        stars_detected = False

                    else:
                        stars_detected = True
                        x0 = max_x
                        y0 = max_y
                        u0 = u0
                        rotation_detected = True

                else:
                    stars_detected = False

                alignment_log(stars_detected, x0, y0, u0)
            # look for large field shift

            # look for large field rotation
            if not stars_detected and not skip_frame:

                alignment_log('Test large rotation')

                ustep = np.arcsin(float(star_std) / comparisons[int(len(comparisons) / 2)][0])

                angles = np.array([np.pi, 0])
                for ff in range(1, int(np.pi / ustep) + 1):
                    angles = np.append(angles, np.pi - ff * ustep)
                    angles = np.append(angles, np.pi + ff * ustep)
                    angles = np.append(angles, 0 - ff * ustep)
                    angles = np.append(angles, 0 + ff * ustep)

                if stars:

                    if show_progress.exit:
                        return None

                    tests = []

                    for rotation in angles:

                        test = 0

                        for star in stars:

                            max_x = star[0]
                            max_y = star[1]
                            circle.set_center((max_x, max_y))
                            canvas.draw()
                            show_progress.update()
                            alignment_log(science_file)
                            alignment_log('LOW SNR Testing star at: ', max_x, max_y, ', with rotation: ', rotation)

                            test = 0

                            for comp in comparisons_snr:

                                check_x = int(max_x + comp[0] * np.cos(rotation + comp[1]))
                                check_y = int(max_y + comp[0] * np.sin(rotation + comp[1]))
                                if 0 < check_x < x_length and 0 < check_y < y_length:
                                    check_sum = np.sum(fits[1].data[check_y - star_std:check_y + star_std + 1,
                                                       check_x - star_std:check_x + star_std + 1])
                                    check_lim = (fits[1].header[mean_key] +
                                                 2 * fits[1].header[std_key]) * ((2 * star_std + 1) ** 2)
                                    if check_sum > check_lim:
                                        test += 1
                                    else:
                                        test -= 1
                                    alignment_log('Check ref. star at: ', check_x, check_y, ', Test: ', test)

                                if abs(test) > check_num_snr:
                                    break

                            tests.append([test, max_x, max_y, rotation])
                            alignment_log([test, max_x, max_y, rotation])

                            if test > check_num_snr:
                                break

                        if test > check_num_snr:
                                break

                    tests.sort()
                    test, max_x, max_y, rotation = tests[-1]
                    if test < check_num_snr:
                        stars_detected = False

                    else:
                        stars_detected = True
                        x0 = max_x
                        y0 = max_y
                        u0 = rotation
                        rotation_detected = True

                else:
                    stars_detected = False

                alignment_log(stars_detected, x0, y0, u0)
            # look for large field rotation

            skip_time += time.time() - delta_skip_time

            if show_progress.exit:
                return None

            if stars_detected:

                if rotation_detected:

                    test_u0 = []
                    test_cos = []
                    test_sin = []

                    for ii in comparisons[:int(check_num + 0.5)]:
                        star = plc.find_single_star(fits[1].data,
                                                    x0 + ii[0] * np.cos(u0 + ii[1]),
                                                    y0 + ii[0] * np.sin(u0 + ii[1]),
                                                    mean=fits[1].header[mean_key],
                                                    std=fits[1].header[std_key],
                                                    burn_limit=burn_limit, star_std=star_std)
                        if star:
                            diff = plc.cartesian_to_polar(star[0], star[1], x0, y0)[1] - ii[1]
                            if diff < 0:
                                diff += 2 * np.pi
                            test_u0.append(diff)
                            test_cos.append(np.cos(diff))
                            test_sin.append(np.sin(diff))

                    if len(test_u0) > 0:
                        test_u0 = np.mean(test_u0)
                        test_cos = np.mean(test_cos)
                        test_sin = np.mean(test_sin)

                        u0 = np.arccos(test_cos)
                        if test_sin < 0:
                            u0 = np.pi + (np.pi - u0)

                fits[1].header.set(align_x0_key, x0)
                fits[1].header.set(align_y0_key, y0)
                fits[1].header.set(align_u0_key, u0)

                circle = mpatches.Circle((x0, y0), circles_diameter, ec='r', fill=False)
                ax.add_patch(circle)
                for ii in comparisons[:int(check_num + 0.5)]:
                    circle = mpatches.Circle((x0 + ii[0] * np.cos(u0 + ii[1]), y0 + ii[0] * np.sin(u0 + ii[1])),
                                             circles_diameter, ec='w', fill=False)
                    ax.add_patch(circle)

                canvas.draw()
                show_progress.update()

            else:

                fits[1].header.set(align_x0_key, False)
                fits[1].header.set(align_y0_key, False)
                fits[1].header.set(align_u0_key, False)

            fits.flush()
            fits.close()

            # counter
            new_percent = round(100 * (counter + 1) / float(len(science)), 1)
            if new_percent != percent:
                lt1 = time.time()
                rm_time = (100 - new_percent) * (lt1 - lt0 - skip_time) / new_percent
                hours = rm_time / 3600.0
                minutes = (hours - int(hours)) * 60
                seconds = (minutes - int(minutes)) * 60

                progress_bar_2['value'] = new_percent
                percent_label_2.configure(text='{0} % ({1}h {2}m {3}s left)'.format(new_percent, int(hours),
                                                                                    int(minutes), int(seconds)))

                percent = new_percent

                show_progress.update()

            if show_progress.exit:
                return None

    def check_visibility():

        star_std = log.read_local_log('alignment', 'star_std')
        all_stars_dict = plc.open_dict('all_stars.pickle')
        stars = np.array(all_stars_dict['all_stars'])

        polar_coords = []
        fits = pf.open(science[0])
        for star in all_stars_dict['all_stars']:
            polar_coords.append(plc.cartesian_to_polar(star[0], star[1], fits[1].header[align_x0_key],
                                                       fits[1].header[align_y0_key]))
        fits.close()

        in_fov = np.ones(len(polar_coords))

        percent = 0
        lt0 = time.time()
        for counter, science_file in enumerate(science):

            if show_progress.exit:
                return None

            in_fov_single = []
            fits = pf.open(science_file)
            if fits[1].header[align_x0_key]:
                ref_x_position = fits[1].header[align_x0_key]
                ref_y_position = fits[1].header[align_y0_key]
                ref_u_position = fits[1].header[align_u0_key]

                for star in polar_coords:

                    cartesian_x = ref_x_position + star[0] * np.cos(ref_u_position + star[1])
                    cartesian_y = ref_y_position + star[0] * np.sin(ref_u_position + star[1])

                    if (cartesian_x > 3 * star_std and cartesian_x < len(fits[1].data[0]) - 3 * star_std and cartesian_y > 3 * star_std
                            and cartesian_y < len(fits[1].data) - 3 * star_std):

                        in_fov_single.append(1)
                    else:
                        in_fov_single.append(0)

                in_fov *= in_fov_single

            fits.close()

            # counter
            new_percent = round(100 * (counter + 1) / float(len(science)), 1)
            if new_percent != percent:
                lt1 = time.time()
                rm_time = (100 - new_percent) * (lt1 - lt0) / new_percent
                hours = rm_time / 3600.0
                minutes = (hours - int(hours)) * 60
                seconds = (minutes - int(minutes)) * 60

                progress_bar_3['value'] = new_percent
                percent_label_3.configure(text='{0} % ({1}h {2}m {3}s left)'.format(new_percent, int(hours),
                                                                                    int(minutes), int(seconds)))

                percent = new_percent

                show_progress.update()

        all_stars_dict['in_fov'] = np.array(in_fov)

        visible_fov_x_min = np.min(stars[np.where(in_fov), 0]) - 3 * star_std
        visible_fov_x_max = np.max(stars[np.where(in_fov), 0]) + 3 * star_std
        visible_fov_y_min = np.min(stars[np.where(in_fov), 1]) - 3 * star_std
        visible_fov_y_max = np.max(stars[np.where(in_fov), 1]) + 3 * star_std

        log.write_local_log('alignment', float(visible_fov_x_min), 'min_x')
        log.write_local_log('alignment', float(visible_fov_y_min), 'min_y')
        log.write_local_log('alignment', float(visible_fov_x_max), 'max_x')
        log.write_local_log('alignment', float(visible_fov_y_max), 'max_y')

        plc.save_dict(all_stars_dict, 'all_stars.pickle')

    def run():

        if not show_progress.exit:
            find_all_stars()

        if log.read_local_log('pipeline', 'alignment_complete'):

            if not show_progress.exit:
                check_visibility()

        else:
            if not show_progress.exit:
                align()

            if not show_progress.exit:
                check_visibility()

        if not show_progress.exit:
            log.write_local_log('pipeline', True, 'alignment_complete')

        show_progress.close()

    # progress window

    show_progress = ProgressWindow('HOPS - Alignment', 0, 0, 5)

    f = Figure()
    f.patch.set_facecolor('white')
    ax = f.add_subplot(111)
    f.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    ax.axis('off')
    canvas = show_progress.FigureCanvasTkAgg(f)
    canvas.get_tk_widget().pack()

    fits = pf.open(science[0], memmap=False)
    ax.imshow(fits[1].data, origin='lower', cmap=cm.Greys_r,
              vmin=fits[1].header[mean_key] + frame_low_std * fits[1].header[std_key],
              vmax=fits[1].header[mean_key] + frame_upper_std * fits[1].header[std_key])
    fits.close()

    frame1 = show_progress.Frame()
    frame1.pack()

    label_1 = Label(frame1, text="Locating all stars in the FOV")
    progress_bar = ttk.Progressbar(frame1, orient=HORIZONTAL,
                                   length=300, maximum=100, mode='determinate', value=0)
    percent_label = Label(frame1, text='0.0 %')

    if log.read_local_log('pipeline', 'alignment_complete'):
        label_2 = Label(frame1, text='Running Alignment: --- Skipping ---')
    else:
        label_2 = Label(frame1, text='Running Alignment: ')
    progress_bar_2 = ttk.Progressbar(frame1, orient=HORIZONTAL,
                                     length=300, maximum=100, mode='determinate', value=0)
    percent_label_2 = Label(frame1, text='0.0 %')

    label_3 = Label(frame1, text=' ')

    label_4 = Label(frame1, text='Aligning all stars in the FOV')
    progress_bar_3 = ttk.Progressbar(frame1, orient=HORIZONTAL,
                                     length=300, maximum=100, mode='determinate', value=0)
    percent_label_3 = Label(frame1, text='0.0 %')

    setup_window(frame1, [
        [[label_1, 0, 2]],
        [[progress_bar, 0, 1, 1, (20, 0)], [percent_label, 1]],
        [[label_2, 0, 2]],
        [[label_3, 0, 2]],
        [[progress_bar_2, 0, 1, 1, (20, 0)], [percent_label_2, 1]],
        [[label_4, 0, 2]],
        [[progress_bar_3, 0, 1, 1, (20, 0)], [percent_label_3, 1]],
        []
    ], main_font='Courier')

    canvas.draw()
    show_progress.after(200, run)
    show_progress.loop()
