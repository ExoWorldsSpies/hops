
import sys

from tkinter import Tk, TclError
from tkinter import Label, Button, Entry, Checkbutton, Scrollbar, Listbox, PhotoImage, Radiobutton, Scale, Frame
from tkinter import StringVar, BooleanVar, DoubleVar, IntVar
from tkinter import DISABLED, NORMAL, END, RIGHT, LEFT, BOTH, Y, HORIZONTAL


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

from hops.hops_tools.windows import *
from hops.hops_tools.logs import log

import glob

def find_fits_files(fits_file):

    fits_list = glob.glob('*{0}*.f*t*'.format(fits_file)) + glob.glob('*{0}*.F*T*'.format(fits_file))
    fits_list = list(np.unique(fits_list))
    fits_list.sort()
    return fits_list


def test_fits_keyword(fits_file, keyword):

    if len(fits_file) == 0:
        return [False, 'No keyword found']

    else:
        try:
            fits_file = find_fits_files(fits_file)[0]

            fits = pf.open(fits_file)

            try:
                fits = [fits['SCI']]
            except KeyError:
                sci_id = 0
                for sci_id in range(len(fits)):
                    try:
                        if (fits[sci_id].data).all():
                            break
                    except:
                        pass
                fits = [fits[sci_id]]

            if fits[0].header[str(keyword)]:
                return [True, 'Keyword found', fits[0].header[str(keyword)]]

            else:
                return [False, 'No keyword found']

        except (KeyError, IndexError):
            return [False, 'No keyword found']


filter_map = {'Clear': 'V', 'Luminance': 'V',
              'U': 'U', 'B': 'B', 'V': 'V', 'R': 'R', 'I': 'I', 'H': 'H', 'J': 'J', 'K': 'K',
              'u': 'u', 'b': 'b', 'v': 'v', 'y': 'y',
              'u\'': 'u,', 'g\'': 'g,', 'r\'': 'r,', 'i\'': 'i,', 'z\'': 'z,',
              'Astrodon ExoPlanet-BB': 'R',
              'UV': 'U', 'Rc': 'R', 'Ic': 'I', 'Re': 'R', 'Ie': 'I', 'Y': 'y,', 'r': 'r,', 'z': 'z,', 'i': 'i,',
              }


def photometry():

    print('Photometry...')

    # get variables

    reduction_directory = log.read_local_log('pipeline', 'reduction_directory')
    light_curve_aperture_file = log.read_local_log('pipeline', 'light_curve_aperture_file')
    photometry_directory_base = log.read_local_log('pipeline', 'photometry_directory')
    photometry_file = log.read_local_log('pipeline', 'photometry_file')
    light_curve_gauss_file = log.read_local_log('pipeline', 'light_curve_gauss_file')
    results_figure = log.read_local_log('pipeline', 'results_figure')
    fov_figure = log.read_local_log('pipeline', 'fov_figure')
    mean_key = log.read_local_log('pipeline_keywords', 'mean_key')
    std_key = log.read_local_log('pipeline_keywords', 'std_key')
    frame_low_std = log.read_local_log('windows', 'frame_low_std')
    frame_upper_std = log.read_local_log('windows', 'frame_upper_std')
    align_x0_key = log.read_local_log('pipeline_keywords', 'align_x0_key')
    align_y0_key = log.read_local_log('pipeline_keywords', 'align_y0_key')
    align_u0_key = log.read_local_log('pipeline_keywords', 'align_u0_key')
    exposure_time_key = log.read_local_log('pipeline_keywords', 'exposure_time_key')
    observation_date_key = log.read_local_log('pipeline_keywords', 'observation_date_key')
    observation_time_key = log.read_local_log('pipeline_keywords', 'observation_time_key')
    mid_exposure = log.read_local_log('photometry', 'mid_exposure')
    star_std = log.read_local_log('alignment', 'star_std')
    star_psf = log.read_local_log('alignment', 'star_psf')
    sky_inner_aperture = log.read_local_log('photometry', 'sky_inner_aperture')
    sky_outer_aperture = log.read_local_log('photometry', 'sky_outer_aperture')
    max_comparisons = log.read_local_log('photometry', 'max_comparisons')
    bin_fits = int(log.read_local_log('reduction', 'bin_fits'))
    burn_limit = int(log.read_local_log('alignment', 'burn_limit')) * bin_fits * bin_fits
    targets_r_position = [log.read_local_log('photometry', 'target_r_position')]
    targets_u_position = [log.read_local_log('photometry', 'target_u_position')]
    targets_aperture = [log.read_local_log('photometry', 'target_aperture')]
    for comparison in range(max_comparisons):
        targets_r_position.append(log.read_local_log('photometry', 'comparison_{0}_r_position'.format(comparison + 1)))
        targets_u_position.append(log.read_local_log('photometry', 'comparison_{0}_u_position'.format(comparison + 1)))
        targets_aperture.append(log.read_local_log('photometry', 'comparison_{0}_aperture'.format(comparison + 1)))

    if star_psf == 0:
        star_psf = star_std

    science = find_fits_files(os.path.join(reduction_directory, '*'))

    def measure():

        gauss_targets_files = []
        gauss_targets_jd = []
        gauss_targets_x_position = []
        gauss_targets_y_position = []
        gauss_targets_x_std = []
        gauss_targets_y_std = []
        gauss_targets_flux = []
        gauss_targets_flux_error = []
        gauss_targets_sky = []
        gauss_targets_sky_error = []
        apperture_targets_files = []
        apperture_targets_jd = []
        apperture_targets_x_position = []
        apperture_targets_y_position = []
        apperture_targets_flux = []
        apperture_targets_flux_error = []
        apperture_targets_sky = []
        apperture_targets_sky_error = []

        # TODO exclude live points

        # for each science_file
        percent = 0
        lt0 = time.time()

        for counter, science_file in enumerate(science):

            fits = pf.open(science_file)

            label_1.configure(text='Running Photometry: {0}'.format(science_file.split(os.sep)[-1]))
            label_1.update()

            if fits[1].header[align_x0_key]:

                # calculate heliocentric julian date

                if observation_date_key == observation_time_key:
                    observation_time = ' '.join(fits[1].header[observation_date_key].split('T'))
                else:
                    observation_time = ' '.join([fits[1].header[observation_date_key].split('T')[0],
                                                 fits[1].header[observation_time_key]])

                exp_time = fits[1].header[exposure_time_key]

                julian_date = plc.UTC(observation_time).jd - 0.5 * float(mid_exposure) * exp_time / 60.0 / 60.0 / 24.0

                ref_x_position = fits[1].header[align_x0_key]
                ref_y_position = fits[1].header[align_y0_key]
                ref_u_position = fits[1].header[align_u0_key]

                gauss_targets_files_test = [science_file]
                gauss_targets_jd_test = [julian_date]
                gauss_targets_x_position_test = []
                gauss_targets_y_position_test = []
                gauss_targets_x_std_test = []
                gauss_targets_y_std_test = []
                gauss_targets_flux_test = []
                gauss_targets_flux_error_test = []
                gauss_targets_sky_test = []
                gauss_targets_sky_error_test = []

                apperture_targets_files_test = [science_file]
                apperture_targets_jd_test = [julian_date]
                apperture_targets_x_position_test = []
                apperture_targets_y_position_test = []
                apperture_targets_flux_test = []
                apperture_targets_flux_error_test = []
                apperture_targets_sky_test = []
                apperture_targets_sky_error_test = []

                skip_gauss = False
                skip_aperture = False

                plotx = ref_x_position + targets_r_position[0] * np.cos(ref_u_position + targets_u_position[0])
                ploty = ref_y_position + targets_r_position[0] * np.sin(ref_u_position + targets_u_position[0])

                for target in range(max_comparisons + 1):

                    if targets_aperture[target] > 0:

                        expected_x = (ref_x_position + targets_r_position[target] *
                         np.cos(ref_u_position + targets_u_position[target]))

                        expected_y = (ref_y_position + targets_r_position[target] *
                         np.sin(ref_u_position + targets_u_position[target]))

                        if (expected_x > 0 and expected_y > 0 and expected_x < len(fits[1].data[0]) and
                                expected_y < len(fits[1].data)):

                            star = plc.find_single_star(fits[1].data,
                                                        (ref_x_position + targets_r_position[target] *
                                                                     np.cos(ref_u_position + targets_u_position[target])),
                                                        (ref_y_position + targets_r_position[target] *
                                                                     np.sin(ref_u_position + targets_u_position[target])),
                                                        mean=fits[1].header[mean_key], std=fits[1].header[std_key],
                                                        burn_limit=burn_limit, star_std=star_std
                                                        )
                            # print(star)

                            if star:

                                x_mean, y_mean, norm, floor, x_sigma, y_sigma, centroid_x, centroid_y = star

                                if target:
                                    plotx = x_mean
                                    ploty = y_mean

                                gauss_targets_x_position_test.append(x_mean)
                                gauss_targets_y_position_test.append(y_mean)
                                gauss_targets_x_std_test.append(x_sigma)
                                gauss_targets_y_std_test.append(y_sigma)
                                gauss_targets_flux_test.append(2 * np.pi * norm * x_sigma * y_sigma)
                                gauss_targets_flux_error_test.append(
                                    np.sqrt(
                                        np.abs(2 * np.pi * norm * x_sigma * y_sigma) +
                                        2 * np.abs(9 * floor * x_sigma * y_sigma)
                                    )
                                )
                                gauss_targets_sky_test.append(floor)
                                gauss_targets_sky_error_test.append(np.sqrt(floor))

                                try:

                                    sky_area_1 = int(round(sky_inner_aperture * targets_aperture[target]))
                                    sky_area_2 = int(round(sky_outer_aperture * targets_aperture[target]))
                                    sky_area_2 = max(sky_area_2, sky_area_1 + 3)

                                    sky_area = fits[1].data[int(y_mean) - sky_area_2:int(y_mean) + sky_area_2 + 1,
                                                            int(x_mean) - sky_area_2:int(x_mean) + sky_area_2 + 1]

                                    sky_area = np.ones_like(sky_area) * sky_area

                                    sky_center = int(len(sky_area) / 2)

                                    sky_area[sky_center - sky_area_1:sky_center + sky_area_1 + 1,
                                             sky_center - sky_area_1:sky_center + sky_area_1 + 1] = -100000

                                    sky_area = sky_area[np.where((sky_area!=-100000) * (sky_area < fits[1].header[mean_key] + 3 * fits[1].header[
                                                                 std_key]))]

                                    sky_total = np.sum(sky_area)
                                    sky = np.sum(sky_total) / sky_area.size
                                    sky_err = np.sqrt(np.abs(sky_total)) / sky_area.size

                                    flux_area = fits[1].data[int(y_mean - targets_aperture[target]) - 2:
                                                             int(y_mean + targets_aperture[target]) + 3,
                                                             int(x_mean - targets_aperture[target]) - 2:
                                                             int(x_mean + targets_aperture[target]) + 3]

                                    flux_area_x, flux_area_y = np.meshgrid(
                                        np.arange(max(0, int(x_mean - targets_aperture[target]) - 2),
                                                  min(len(fits[1].data[0]), int(x_mean + targets_aperture[target]) + 3), 1) + 0.5,
                                        np.arange(max(0, int(y_mean - targets_aperture[target]) - 2),
                                                  min(len(fits[1].data), int(y_mean + targets_aperture[target]) + 3), 1) + 0.5)

                                    flux_pixels = np.concatenate(np.swapaxes([flux_area, flux_area_x, flux_area_y], 0, 2))

                                    flux = 0
                                    flux_err = 0
                                    for pixel in flux_pixels:
                                        overlap = plc.pixel_to_aperture_overlap(pixel[1], pixel[2], x_mean, y_mean,
                                                                                targets_aperture[target])
                                        flux += (pixel[0] - sky) * overlap
                                        flux_err = np.sqrt(flux_err ** 2 + (np.abs(pixel[0]) + sky_err**2) * (overlap**2))

                                    apperture_targets_x_position_test.append(x_mean)
                                    apperture_targets_y_position_test.append(y_mean)
                                    apperture_targets_flux_test.append(flux)
                                    apperture_targets_flux_error_test.append(flux_err)
                                    apperture_targets_sky_test.append(sky)
                                    apperture_targets_sky_error_test.append(sky_err)

                                except:
                                    skip_aperture = True

                            else:
                                skip_gauss = True
                                skip_aperture = True

                        else:
                            skip_gauss = True
                            skip_aperture = True

                if not skip_gauss:
                    gauss_targets_files += gauss_targets_files_test
                    gauss_targets_jd += gauss_targets_jd_test
                    gauss_targets_x_position += gauss_targets_x_position_test
                    gauss_targets_y_position += gauss_targets_y_position_test
                    gauss_targets_x_std += gauss_targets_x_std_test
                    gauss_targets_y_std += gauss_targets_y_std_test
                    gauss_targets_flux += gauss_targets_flux_test
                    gauss_targets_flux_error += gauss_targets_flux_error_test
                    gauss_targets_sky += gauss_targets_sky_test
                    gauss_targets_sky_error += gauss_targets_sky_error_test
                else:
                    print('Skipping Gauss for:', science_file)

                if not skip_aperture:
                    apperture_targets_files += apperture_targets_files_test
                    apperture_targets_jd += apperture_targets_jd_test
                    apperture_targets_x_position += apperture_targets_x_position_test
                    apperture_targets_y_position += apperture_targets_y_position_test
                    apperture_targets_flux += apperture_targets_flux_test
                    apperture_targets_flux_error += apperture_targets_flux_error_test
                    apperture_targets_sky += apperture_targets_sky_test
                    apperture_targets_sky_error += apperture_targets_sky_error_test
                else:
                    print('Skipping aperture for:', science_file)

            # counter

                ax.cla()
                ax.imshow(fits[1].data[int(ploty - 3 * targets_aperture[0]): int(ploty + 3 * targets_aperture[0]),
                                       int(plotx - 3 * targets_aperture[0]): int(plotx + 3 * targets_aperture[0])],
                          origin='lower',
                          extent=(int(plotx - 3 * targets_aperture[0]), int(plotx + 3 * targets_aperture[0]),
                                  int(ploty - 3 * targets_aperture[0]), int(ploty + 3 * targets_aperture[0])),
                          cmap=cm.Greys_r,
                          vmin=fits[1].header[mean_key] + frame_low_std * fits[1].header[std_key],
                          vmax=fits[1].header[mean_key] + frame_upper_std * fits[1].header[std_key])

                ax.set_xlim(plotx - 2 * targets_aperture[0], plotx + 2 * targets_aperture[0])
                ax.set_ylim(ploty - 2 * targets_aperture[0], ploty + 2 * targets_aperture[0])
                ax.axis('off')

                if not skip_gauss:
                    circle = mpatches.Circle((plotx, ploty), targets_aperture[0], ec='r', fill=False)
                    ax.add_patch(circle)

                canvas.draw()

            new_percent = round(100 * (counter + 1) / float(len(science)), 1)
            if new_percent != percent:
                lt1 = time.time()
                rm_time = (100 - new_percent) * (lt1 - lt0) / new_percent
                hours = rm_time / 3600.0
                minutes = (hours - int(hours)) * 60
                seconds = (minutes - int(minutes)) * 60

                progress_bar_1['value'] = new_percent
                percent_label_1.configure(text='{0} % ({1}h {2}m {3}s left)'.format(new_percent, int(hours),
                                                                                    int(minutes), int(seconds)))

                percent = new_percent

            if show_progress.exit:
                break

            if counter + 1 == len(science):
                log.write_local_log('pipeline', True, 'photometry_complete')

            show_progress.update()

        return [
            gauss_targets_files,
            gauss_targets_jd,
            gauss_targets_x_position,
            gauss_targets_y_position,
            gauss_targets_x_std,
            gauss_targets_y_std,
            gauss_targets_flux,
            gauss_targets_flux_error,
            gauss_targets_sky,
            gauss_targets_sky_error,
            apperture_targets_files,
            apperture_targets_jd,
            apperture_targets_x_position,
            apperture_targets_y_position,
            apperture_targets_flux,
            apperture_targets_flux_error,
            apperture_targets_sky,
            apperture_targets_sky_error,
        ]

    def plot_lcs(measurements):

        (gauss_targets_files, gauss_targets_jd, gauss_targets_x_position, gauss_targets_y_position, gauss_targets_x_std,
         gauss_targets_y_std, gauss_targets_flux, gauss_targets_flux_error, gauss_targets_sky, gauss_targets_sky_error,
         apperture_targets_files, apperture_targets_jd, apperture_targets_x_position, apperture_targets_y_position,
         apperture_targets_flux, apperture_targets_flux_error, apperture_targets_sky, apperture_targets_sky_error) = measurements

        # save results, create photometry directory and move results there
        comparisons_number = len(gauss_targets_x_position) // len(gauss_targets_files) - 1

        measurements_number = len(gauss_targets_files)
        targets_number = len(gauss_targets_x_position) // len(gauss_targets_files)

        gauss_targets_jd = np.array(gauss_targets_jd)
        targets_jd = np.array(gauss_targets_jd)

        targets_x_position = np.swapaxes(np.reshape(gauss_targets_x_position,
                                                    (measurements_number, targets_number)), 0, 1)
        targets_y_position = np.swapaxes(np.reshape(gauss_targets_y_position,
                                                    (measurements_number, targets_number)), 0, 1)
        targets_x_std = np.swapaxes(np.reshape(gauss_targets_x_std, (measurements_number, targets_number)), 0, 1)
        targets_y_std = np.swapaxes(np.reshape(gauss_targets_y_std, (measurements_number, targets_number)), 0, 1)
        targets_flux = np.swapaxes(np.reshape(gauss_targets_flux, (measurements_number, targets_number)), 0, 1)
        targets_gauss_flux = np.swapaxes(np.reshape(gauss_targets_flux, (measurements_number, targets_number)), 0, 1)
        targets_flux_error = np.swapaxes(np.reshape(gauss_targets_flux_error, (measurements_number, targets_number)), 0, 1)
        targets_sky = np.swapaxes(np.reshape(gauss_targets_sky, (measurements_number, targets_number)), 0, 1)
        targets_sky_error = np.swapaxes(np.reshape(gauss_targets_sky_error, (measurements_number, targets_number)), 0, 1)

        targets_results = [targets_jd] + (list(targets_x_position) + list(targets_y_position) +
                                          list(targets_x_std) + list(targets_y_std) + list(targets_flux) +
                                          list(targets_flux_error) + list(targets_sky) + list(targets_sky_error))

        np.savetxt(photometry_file.replace('.txt', '_g.txt'), np.swapaxes(targets_results, 0, 1))

        np.savetxt(light_curve_gauss_file,
                   np.swapaxes([targets_jd,
                                targets_flux[0] / np.sum(targets_flux[1:], 0),
                                np.sqrt(
                                    (targets_flux_error[0] / targets_flux[0]) ** 2 +
                                    (np.sqrt(np.sum(targets_flux_error[1:]**2, 0)) / np.sum(targets_flux[1:], 0)) ** 2
                                ) * np.abs(targets_flux[0] / np.sum(targets_flux[1:], 0))
                                ], 0, 1))

        gflux = targets_flux[0] / np.sum(targets_flux[1:], 0)
        diff = np.abs(gflux[1:] - gflux[:-1])
        gauss_scatter = np.std(diff)

        measurements_number = len(apperture_targets_files)
        targets_number = len(apperture_targets_x_position) // len(apperture_targets_files)

        targets_jd = np.array(apperture_targets_jd)

        targets_x_position = np.swapaxes(np.reshape(apperture_targets_x_position,
                                                    (measurements_number, targets_number)), 0, 1)
        targets_y_position = np.swapaxes(np.reshape(apperture_targets_y_position,
                                                    (measurements_number, targets_number)), 0, 1)
        targets_flux = np.swapaxes(np.reshape(apperture_targets_flux, (measurements_number, targets_number)), 0, 1)
        targets_flux_error = np.swapaxes(np.reshape(apperture_targets_flux_error, (measurements_number, targets_number)), 0, 1)
        targets_sky = np.swapaxes(np.reshape(apperture_targets_sky, (measurements_number, targets_number)), 0, 1)
        targets_sky_error = np.swapaxes(np.reshape(apperture_targets_sky_error, (measurements_number, targets_number)), 0, 1)

        targets_results = [targets_jd] + (list(targets_x_position) + list(targets_y_position) + list(targets_flux) +
                                          list(targets_flux_error) + list(targets_sky) + list(targets_sky_error))

        np.savetxt(photometry_file.replace('.txt', '_a.txt'), np.swapaxes(targets_results, 0, 1))

        np.savetxt(light_curve_aperture_file,
                   np.swapaxes([targets_jd,
                                targets_flux[0] / np.sum(targets_flux[1:], 0),
                                np.sqrt(
                                    (targets_flux_error[0] / targets_flux[0]) ** 2 +
                                    (np.sqrt(np.sum(targets_flux_error[1:] ** 2, 0)) / np.sum(targets_flux[1:], 0)) ** 2
                                ) * np.abs(targets_flux[0] / np.sum(targets_flux[1:], 0))
                                ], 0, 1))

        aflux = targets_flux[0] / np.sum(targets_flux[1:], 0)
        diff = np.abs(aflux[1:] - aflux[:-1])
        apperture_scatter = np.std(diff)

        photometry_directory = photometry_directory_base

        if not os.path.isdir(photometry_directory):
            os.mkdir(photometry_directory)
        else:
            fi = 2
            while os.path.isdir('{0}_{1}'.format(photometry_directory, str(fi))):
                fi += 1
            photometry_directory = '{0}_{1}'.format(photometry_directory, str(fi))
            os.mkdir(photometry_directory)

        # plot

        root = ProgressWindow('HOPS - Photometry')
        f = Figure()
        f.set_figwidth(9)
        f.set_figheight(0.8 * root.root.winfo_screenheight() / f.get_dpi())
        if comparisons_number > 1:
            ax = f.add_subplot(comparisons_number + 1, 1, 1)
        else:
            ax = f.add_subplot(1, 1, 1)

        f.patch.set_facecolor('white')
        canvas = root.FigureCanvasTkAgg(f)
        canvas.get_tk_widget().pack()
        root.NavigationToolbar2Tk(canvas)

        ax.plot(targets_jd - targets_jd[0], targets_flux[0] / np.sum(targets_flux[1:], 0)
                / np.median(targets_flux[0] / np.sum(targets_flux[1:], 0)), 'ko', ms=3, label='Aperture')
        ax.plot(gauss_targets_jd - gauss_targets_jd[0], targets_gauss_flux[0] / np.sum(targets_gauss_flux[1:], 0)
                / np.median(targets_gauss_flux[0] / np.sum(targets_gauss_flux[1:], 0)), 'ro', ms=3, mec='r', label='PSF')
        ax.tick_params(labelbottom=False)
        ax.set_title(r'$\mathrm{Target}$')
        ax.legend(loc=(0, 1))

        if comparisons_number > 1:

            all_relative_gauss = []
            all_relative_aperture = []

            for comp in range(comparisons_number):
                test_aperture_flux = list(targets_flux[1:])
                test_gauss_flux = list(targets_gauss_flux[1:])
                del test_aperture_flux[comp]
                del test_gauss_flux[comp]
                all_relative_aperture.append(targets_flux[1:][comp] / np.sum(test_aperture_flux, 0))
                all_relative_gauss.append(targets_gauss_flux[1:][comp] / np.sum(test_gauss_flux, 0))

            for comp in range(comparisons_number):

                ax = f.add_subplot(comparisons_number + 1, 1, comp + 2)
                ax.plot(targets_jd - targets_jd[0],
                        all_relative_aperture[comp] / np.median(all_relative_aperture[comp]), 'ko', ms=3)
                ax.plot(gauss_targets_jd - gauss_targets_jd[0],
                        all_relative_gauss[comp] / np.median(all_relative_gauss[comp]), 'ro', ms=3, mec='r')
                ax.tick_params(labelbottom=False)
                f.text(0.881, 0.07 + (comparisons_number - comp - 0.5) * (1 - 0.07 - 0.12) / (comparisons_number + 1),
                       r'${0}{1}{2}$'.format('\mathrm{', "Comparison \, {0}".format(comp + 1), '}'),
                       ha='left', va='center')

        ax.tick_params(labelbottom=True)
        ax.set_xlabel(r'$\Delta t \ \mathrm{[days]}$', fontsize=12)
        f.text(0.03, 0.5, r'$\mathrm{relative} \ \mathrm{flux}$', fontsize=12,
               ha='center', va='center', rotation='vertical')
        f.subplots_adjust(0.1, 0.07, 1 - 0.12, 1 - 0.1, 0, 0)
        f.set_figheight(2 * (comparisons_number + 1))
        f.savefig(results_figure, dpi=200)

        shutil.move(fov_figure, '{0}/{1}'.format(photometry_directory, fov_figure))
        shutil.move(results_figure, '{0}/{1}'.format(photometry_directory, results_figure))
        photometry_file_g = photometry_file.replace('.txt', '_g.txt')
        photometry_file_a = photometry_file.replace('.txt', '_a.txt')
        shutil.move(photometry_file_g, '{0}/{1}'.format(photometry_directory, photometry_file_g))
        shutil.move(photometry_file_a, '{0}/{1}'.format(photometry_directory, photometry_file_a))
        shutil.move(light_curve_gauss_file, '{0}/{1}'.format(photometry_directory, light_curve_gauss_file))
        shutil.move(light_curve_aperture_file, '{0}/{1}'.format(photometry_directory, light_curve_aperture_file))
        shutil.copy('log.yaml', '{0}/log.yaml'.format(photometry_directory))

        fits = pf.open(science[np.random.randint(len(science))])
        exp_time = round(fits[1].header[exposure_time_key], 1)

        ra_dec_string = log.read_local_log('photometry', 'target_ra_dec')
        ra_dec_string = ra_dec_string.replace(':', ' ').split(' ')
        target = plc.Target(plc.Hours(*ra_dec_string[:3]), plc.Degrees(*ra_dec_string[3:]))

        ecc_planet = plc.find_nearest(target)

        phot_filter = 'None detected'
        for key in log.read_local_log_profile('filter_key').split(','):
            check_filter = test_fits_keyword(fits, key)
            if check_filter[0]:
                if check_filter[2] in filter_map:
                    phot_filter = check_filter[2]
                    break
            if log.read_local_log_profile('filter') != '':
                phot_filter = log.read_local_log_profile('filter')

        if gauss_scatter < apperture_scatter:
            files_to_upload = ['PHOTOMETRY_GAUSS.txt', 'PHOTOMETRY_APERTURE.txt']
        else:
            files_to_upload = ['PHOTOMETRY_APERTURE.txt', 'PHOTOMETRY_GAUSS.txt']

        w = open('{0}/ExoClock_info.txt'.format(photometry_directory), 'w')
        w.write('\n'.join([
            'The ExoClock Project is an effort to keep the ephemerides of exoplanets as precise as \n'
            'possible for planning future observations. If you have observed an exoplanet you can\n'
            'contribute your observation at: \n\nhttps://www.exoclock.space\n\n'
            'File to upload: {0} \n(this is a suggestion based on the scatter \nof your light curves, '
            'you can also try \nuploading {1})'.format(*files_to_upload),
            '',
            'Planet: {0} \n(this is the closest known exoplanet found \nin the catalogue, if this is not the '
            'target \nyou observed, please ignore)'.format(ecc_planet.planet.name),
            '',
            'Time format: JD_UTC \n(UTC-based Julian date)',
            '',
            'Time stamp: Exposure start \n(the time produced refers to the beginning of each exposure)',
            '',
            'Flux format: Flux \n(flux of target over summed flux of comparisons)',
            '',
            'Filter: {0}'.format(phot_filter),
            '',
            'Exposure time in seconds: {0}'.format(exp_time),
        ]))
        w.close()

        shutil.copy(log.photometry_output_description, photometry_directory)

        root.show()

        while not root.exit:
            root.update()

        root.close()

    def run():

        measurements = measure()
        if not show_progress.exit:
            plot_lcs(measurements)
            if not show_progress.exit:
                log.write_local_log('pipeline', True, 'photometry_complete')

        show_progress.close()

    # progress window

    show_progress = ProgressWindow('HOPS - Photometry', 0, 0, 5)

    f = Figure(figsize=(3, 3))
    f.patch.set_facecolor('white')
    ax = f.add_subplot(111)
    f.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    ax.axis('off')
    canvas = show_progress.FigureCanvasTkAgg(f)
    canvas.get_tk_widget().pack()

    fits_file = pf.open(science[0], memmap=False)
    try:
        fits = [fits_file['SCI']]
    except KeyError:
        sci_id = 0
        for sci_id in range(len(fits_file)):
            try:
                if (fits_file[sci_id].data).all():
                    break
            except:
                pass
        fits = [0, fits_file[sci_id]]

    ref_x_position = fits[1].header[align_x0_key]
    ref_y_position = fits[1].header[align_y0_key]
    ref_u_position = fits[1].header[align_u0_key]

    plotx = ref_x_position + targets_r_position[0] * np.cos(ref_u_position + targets_u_position[0])
    ploty = ref_y_position + targets_r_position[0] * np.sin(ref_u_position + targets_u_position[0])

    ax.cla()
    ax.imshow(fits[1].data[int(ploty - 3 * targets_aperture[0]): int(ploty + 3 * targets_aperture[0]),
              int(plotx - 3 * targets_aperture[0]): int(plotx + 3 * targets_aperture[0])],
              origin='lower',
              extent=(int(plotx - 3 * targets_aperture[0]), int(plotx + 3 * targets_aperture[0]),
                      int(ploty - 3 * targets_aperture[0]), int(ploty + 3 * targets_aperture[0])),
              cmap=cm.Greys_r,
              vmin=fits[1].header[mean_key] + frame_low_std * fits[1].header[std_key],
              vmax=fits[1].header[mean_key] + frame_upper_std * fits[1].header[std_key])

    ax.set_xlim(plotx - 2 * targets_aperture[0], plotx + 2 * targets_aperture[0])
    ax.set_ylim(ploty - 2 * targets_aperture[0], ploty + 2 * targets_aperture[0])
    ax.axis('off')

    circle = mpatches.Circle((plotx, ploty), targets_aperture[0], ec='r', fill=False)
    ax.add_patch(circle)

    fits_file.close()

    frame1 = show_progress.Frame()
    frame1.pack()

    label_1 = Label(frame1, text="Running Photometry: {0}".format(os.path.split(science[0])[1]))
    progress_bar_1 = ttk.Progressbar(frame1, orient=HORIZONTAL,
                                     length=200, maximum=100, mode='determinate', value=0)
    percent_label_1 = Label(frame1, text='0.0 % (0h 0m 0s left)')

    setup_window(frame1, [
        [],
        [[label_1, 0]],
        [[progress_bar_1, 0]],
        [[percent_label_1, 0]],
        []
    ], main_font='Courier')

    canvas.draw()
    show_progress.after(200, run)
    show_progress.loop()

