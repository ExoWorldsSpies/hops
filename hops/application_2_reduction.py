
import os
import time
import numpy as np
import shutil
import hops.pylightcurve41 as plc
import sys

from astropy.io import fits as pf

from hops.hops_tools.fits import *
from hops.hops_tools.image_analysis import image_mean_std, image_burn_limit, image_psf, bin_frame
from hops.application_windows import MainWindow


class ReductionWindow(MainWindow):

    def __init__(self, log):

        MainWindow.__init__(self, log, name='HOPS - Reduction', position=2)

        self.file_info = []

        self.bias_files = find_fits_files(self.log.get_param('bias_files'))
        self.bias_frames = []
        self.bias_frames_exp = []
        self.bias_counter = 0

        self.dark_files = find_fits_files(self.log.get_param('dark_files'))
        self.dark_frames = []
        self.dark_counter = 0

        self.darkf_files = find_fits_files(self.log.get_param('darkf_files'))
        self.darkf_frames = []
        self.darkf_counter = 0

        self.flat_files = find_fits_files(self.log.get_param('flat_files'))
        self.flat_frames = []
        self.flat_counter = 0

        self.science_files = find_fits_files(self.log.get_param('observation_files'))
        self.all_frames = {}
        self.science_counter = 0
        self.ref_stars = []
        self.psf = 10

        location_string = self.log.get_param('location').split(' ')
        self.observatory = plc.Observatory(plc.Degrees(location_string[0]), plc.Degrees(location_string[1]))
        ra_dec_string = self.log.get_param('target_ra_dec').split(' ')
        self.target = plc.FixedTarget(plc.Hours(ra_dec_string[0]), plc.Degrees(ra_dec_string[1]))

        self.filter = self.log.get_param('filter')

        test_fits_name = self.science_files[0]
        test_fits_data, test_fits_header = self.get_fits_data_and_header(test_fits_name)

        t0 = time.time()
        _ = plc.mean_std_from_median_mad(test_fits_data)
        self.fr_time = int(2000 * (time.time()-t0))

        self.progress_figure = self.FitsWindow(fits_data=test_fits_data, fits_header=test_fits_header,
                                               input_name=self.science_files[0])
        self.progress_bias = self.Progressbar(task="Loading bias frames")
        self.progress_dark = self.Progressbar(task="Loading dark frames")
        self.progress_darkf = self.Progressbar(task="Loading dark-flat frames")
        self.progress_flat = self.Progressbar(task="Loading flat frames")
        self.progress_science = self.Progressbar(task="Reducing science frames and calculating statistics")
        self.progress_science_loop = self.CheckButton(text='Show all frames', initial=0)

        setup_window = [
            [[self.progress_figure, 0, 2]]
        ]

        if len(self.bias_files) > 0:
            setup_window.append([[self.progress_bias, 0, 2]])
        if len(self.dark_files) > 0:
            setup_window.append([[self.progress_dark, 0, 2]])
        if len(self.darkf_files) > 0:
            setup_window.append([[self.progress_darkf, 0, 2]])
        if len(self.flat_files) > 0:
            setup_window.append([[self.progress_flat, 0, 2]])

        setup_window += [
            [[self.progress_science, 0], [self.progress_science_loop, 1]],
            [[self.Button(text='STOP REDUCTION & RETURN TO MAIN MENU', command=self.trigger_exit), 0, 2]],
            []
        ]

        self.setup_window(setup_window)

        self.set_close_button_function(self.trigger_exit)

    def get_fits_data_and_header(self, fits_name):

        fits_data, fits_header = get_fits_data_and_header(fits_name)

        if self.file_info == []:
            self.file_info = [len(fits_data[0]), len(fits_data), fits_header['BITPIX']]
        else:
            if len(fits_data[0]) != self.file_info[0] or self.file_info[1] != len(fits_data) or self.file_info[2] != fits_header['BITPIX']:
                self.showinfo(
                    'Inconsistent image size',
                    'File {0} has different size from the other files in the dataset:\n'
                    'x-size: {1} pixels\ny-size: {2} pixels\n{3}-bit\nReduction will terminate, please check your files!'.format(fits_name, len(fits_data[0]), len(fits_data), fits_header['BITPIX']))
                self.show()
                self.exit = True
                self.after(self.reduce_science)

        return fits_data, fits_header

    # define functions

    def run_reduction(self):

        if self.log.get_param('reduction_complete'):
            if self.askyesno('Overwrite files', 'Reduction has been completed, do you want to run again?'):
                self.log.set_param('reduction_complete', False)
                self.log.save_local_log()
            else:
                self.log.set_param('proceed', True)

        if not self.log.get_param('reduction_complete'):

            if not os.path.isdir(self.log.reduction_directory):
                os.mkdir(self.log.reduction_directory)
            else:
                shutil.rmtree(self.log.reduction_directory)
                os.mkdir(self.log.reduction_directory)

            self.after(self.get_bias)

        else:
            self.close()

    # reduction routines

    def get_bias(self):

        if self.exit or len(self.bias_files) == 0:
            self.after(self.get_master_bias)

        else:

            if self.bias_counter == 0:
                self.progress_bias.initiate(len(self.bias_files))

            fits_data, fits_header = self.get_fits_data_and_header(self.bias_files[self.bias_counter])
            self.bias_frames.append(np.ones_like(fits_data) * fits_data)

            try:
                bias_exptime = fits_header[self.log.get_param('exposure_time_key')]
            except:
                bias_exptime = 0
            self.bias_frames_exp.append(bias_exptime)

            print('{0}: median = {1}, exp.time = {2}'.format(self.bias_files[self.bias_counter], np.nanmedian(self.bias_frames[-1]), bias_exptime))

            self.progress_bias.update()
            self.bias_counter += 1

            if self.bias_counter >= len(self.bias_files):
                self.progress_bias.show_message('Calculating master bias...')
                self.after(self.get_master_bias)
            else:
                self.after(self.get_bias)

    def get_master_bias(self):

        if self.exit:
            self.after(self.get_dark)

        else:

            if len(self.bias_frames) > 0:

                consistent_exp_time = self.bias_frames_exp == np.median(self.bias_frames_exp)

                self.bias_frames = [self.bias_frames[ff] for ff in range(len(self.bias_frames))  if consistent_exp_time[ff]]

                if self.log.get_param('master_bias_method') == 'median':
                    self.master_bias = np.array([np.nanmedian([xx[ff] for xx in self.bias_frames], 0) for ff in range(len(self.bias_frames[0]))])
                elif self.log.get_param('master_bias_method') == 'mean':
                    self.master_bias = np.array([np.nanmean([xx[ff] for xx in self.bias_frames], 0) for ff in range(len(self.bias_frames[0]))])
                else:
                    self.master_bias = np.array([np.nanmedian([xx[ff] for xx in self.bias_frames], 0) for ff in range(len(self.bias_frames[0]))])

                self.bias_frames_exp = np.median(self.bias_frames_exp)

            else:
                self.master_bias = 0.0
                self.bias_frames_exp = 0.0

            print('Median Bias: ', round(np.nanmedian(self.master_bias), 3))
            print('Bias exp. time: ', self.bias_frames_exp)
            self.progress_bias.show_message('Calculating master bias... Completed!')

            self.after(self.get_dark)

    def get_dark(self):

        if self.exit or len(self.dark_files) == 0:
            self.after(self.get_master_dark)

        else:

            if self.dark_counter == 0:
                self.progress_dark.initiate(len(self.dark_files))

            fits_data, fits_header = self.get_fits_data_and_header(self.dark_files[self.dark_counter])
            dark_frame = np.ones_like(fits_data) * fits_data
            self.dark_frames.append((dark_frame - self.master_bias) / (fits_header[self.log.get_param('exposure_time_key')] - self.bias_frames_exp))

            print('{0}: median = {1}'.format(self.dark_files[self.dark_counter], np.nanmedian(self.dark_frames[-1])))

            self.progress_dark.update()
            self.dark_counter += 1

            if self.dark_counter >= len(self.dark_files):
                self.progress_dark.show_message('Calculating master dark...')
                self.after(self.get_master_dark)
            else:
                self.after(self.get_dark)

    def get_master_dark(self):

        if self.exit:
            self.after(self.get_flat)
        else:

            if len(self.dark_frames) > 0:
                if self.log.get_param('master_dark_method') == 'median':
                    self.master_dark = np.array([np.nanmedian([xx[ff] for xx in self.dark_frames], 0) for ff in range(len(self.dark_frames[0]))])
                elif self.log.get_param('master_dark_method') == 'mean':
                    self.master_dark = np.array([np.nanmean([xx[ff] for xx in self.dark_frames], 0) for ff in range(len(self.dark_frames[0]))])
                else:
                    self.master_dark = np.array([np.nanmedian([xx[ff] for xx in self.dark_frames], 0) for ff in range(len(self.dark_frames[0]))])
            else:
                self.master_dark = 0.0

            print('Median Dark: ', round(np.nanmedian(self.master_dark), 3))
            self.progress_dark.show_message('Calculating master dark... Completed!')

            self.after(self.get_darkf)

    def get_darkf(self):

        if self.exit or len(self.darkf_files) == 0:
            self.after(self.get_master_darkf)

        else:

            if self.darkf_counter == 0:
                self.progress_darkf.initiate(len(self.darkf_files))

            fits_data, fits_header = self.get_fits_data_and_header(self.darkf_files[self.darkf_counter])
            darkf_frame = np.ones_like(fits_data) * fits_data
            self.darkf_frames.append((darkf_frame - self.master_bias) / (fits_header[self.log.get_param('exposure_time_key')] - self.bias_frames_exp))

            print('{0}: median = {1}'.format(self.darkf_files[self.darkf_counter], np.nanmedian(self.darkf_frames[-1])))

            self.progress_darkf.update()
            self.darkf_counter += 1

            if self.darkf_counter >= len(self.darkf_files):
                self.progress_darkf.show_message('Calculating master dark-flat...')
                self.after(self.get_master_darkf)
            else:
                self.after(self.get_darkf)

    def get_master_darkf(self):

        if self.exit:
            self.after(self.get_flat)
        else:

            if len(self.darkf_frames) > 0:
                if self.log.get_param('master_darkf_method') == 'median':
                    self.master_darkf = np.array([np.nanmedian([xx[ff] for xx in self.darkf_frames], 0) for ff in range(len(self.darkf_frames[0]))])
                elif self.log.get_param('master_darkf_method') == 'mean':
                    self.master_darkf = np.array([np.nanmean([xx[ff] for xx in self.darkf_frames], 0) for ff in range(len(self.darkf_frames[0]))])
                else:
                    self.master_darkf = np.array([np.nanmedian([xx[ff] for xx in self.darkf_frames], 0) for ff in range(len(self.darkf_frames[0]))])
            else:
                self.master_darkf = self.master_dark

            print('Median Dark-Flat: ', round(np.nanmedian(self.master_darkf), 3))
            self.progress_darkf.show_message('Calculating master dark-flat... Completed!')

            self.after(self.get_flat)

    def get_flat(self):

        if self.exit or len(self.flat_files) == 0:
            self.after(self.get_master_flat)

        else:

            if self.flat_counter == 0:
                self.progress_flat.initiate(len(self.flat_files))

            fits_data, fits_header = self.get_fits_data_and_header(self.flat_files[self.flat_counter])
            flat_frame = np.ones_like(fits_data) * fits_data

            self.flat_frames.append(
                flat_frame - self.master_bias - (fits_header[self.log.get_param('exposure_time_key')] - self.bias_frames_exp) * self.master_darkf)

            print('{0}: median = {1}'.format(self.flat_files[self.flat_counter], np.nanmedian(self.flat_frames[-1])))

            self.progress_flat.update()
            self.flat_counter += 1

            if self.flat_counter >= len(self.flat_files):
                self.progress_flat.show_message('Calculating master flat...')
                self.after(self.get_master_flat)
            else:
                self.after(self.get_flat)

    def get_master_flat(self):

        if self.exit:
            self.after(self.reduce_science)

        else:

            if len(self.flat_frames) > 0:
                if self.log.get_param('master_flat_method') == 'mean':
                    flat_frames = [ff / np.nanmean(ff) for ff in self.flat_frames]
                    self.master_flat = np.array([np.nanmean([xx[ff] for xx in flat_frames], 0) for ff in range(len(flat_frames[0]))])
                else:
                    flat_frames = [ff / np.nanmedian(ff) for ff in self.flat_frames]
                    self.master_flat = np.array([np.nanmedian([xx[ff] for xx in flat_frames], 0) for ff in range(len(flat_frames[0]))])
                print('Median Flat: ', round(np.nanmedian(self.master_flat), 3))
                self.master_flat = self.master_flat / np.nanmedian(self.master_flat)
                self.master_flat = np.where(self.master_flat == 0, 1, self.master_flat)
            else:
                self.master_flat = 1.0

            self.progress_flat.show_message('Calculating master flat... Completed!')
            sys.setrecursionlimit(100 * len(self.science_files))
            self.after(self.reduce_science)

    def reduce_science(self):

        timing = False
        # timing = True

        # correct each observation_files file

        if self.exit:
            self.after(self.save)

        else:

            if self.science_counter == 0:
                self.progress_science.initiate(len(self.science_files))

            science_file = self.science_files[self.science_counter]

            # correct it with master bias_files, master dark_files and master flat_files
            t00 = time.time()
            t0 = time.time()
            fits_data, fits_header = self.get_fits_data_and_header(science_file)

            if timing:
                print('Loading: ', time.time()-t0)

            t0 = time.time()

            saturation = image_burn_limit(fits_header, key=self.log.hops_saturation_key)
            exp_time = float(fits_header[self.log.get_param('exposure_time_key')])
            data_frame = np.ones_like(fits_data) * fits_data
            dq_frame = np.where(data_frame == saturation, 1, 0)
            data_frame = (data_frame - self.master_bias - (exp_time - self.bias_frames_exp) * self.master_dark) / self.master_flat
            data_frame[np.where(np.isnan(data_frame))] = 0

            if timing:
                print('Reduction: ', time.time()-t0)

            t0 = time.time()

            crop_x1 = int(max(0, self.log.get_param('crop_x1')))
            crop_x2 = int(min(self.log.get_param('crop_x2'), len(data_frame[0])))
            crop_y1 = int(max(0, self.log.get_param('crop_y1')))
            crop_y2 = int(min(self.log.get_param('crop_y2'), len(data_frame)))

            if crop_x2 == 0:
                crop_x2 = len(data_frame[0])
            if crop_y2 == 0:
                crop_y2 = len(data_frame)

            if not (np.array([crop_x1, crop_x2, crop_y1, crop_y2]) == np.array([0, len(data_frame[0]), 0, len(data_frame)])).all():
                data_frame = data_frame[crop_y1: crop_y2]
                data_frame = data_frame[:, crop_x1: crop_x2]
                dq_frame = dq_frame[crop_y1: crop_y2]
                dq_frame = dq_frame[:, crop_x1: crop_x2]

            crop_edge_pixels = int(self.log.get_param('crop_edge_pixels'))
            if crop_edge_pixels > 0:
                data_frame = data_frame[crop_edge_pixels: -crop_edge_pixels, crop_edge_pixels: -crop_edge_pixels]
                dq_frame = dq_frame[crop_edge_pixels: -crop_edge_pixels, crop_edge_pixels: -crop_edge_pixels]

            bin_fits = self.log.get_param('bin_fits')
            if bin_fits > 1:
                data_frame = bin_frame(data_frame, bin_fits)
                saturation = saturation * bin_fits * bin_fits
                dq_frame = bin_frame(dq_frame, bin_fits)
                dq_frame = np.where(dq_frame > 0, 1, 0)

            data_frame[np.where(dq_frame > 0)] = saturation

            if timing:
                print('Binning and cropping: ', time.time()-t0)

            if self.science_counter == 0:
                t0 = time.time()
                _ = plc.mean_std_from_median_mad(data_frame)
                self.fr_time = int(1000 * (time.time()-t0))

            t0 = time.time()
            mean, std = image_mean_std(data_frame, samples=10000, mad_filter=5.0)
            if timing:
                print('SKY: ', time.time()-t0)

            t0 = time.time()
            psf = image_psf(data_frame, fits_header, mean, std, saturation)
            if np.isnan(psf):
                psf = 10
                skip = True
            else:
                skip = False
            if timing:
                print('PSF: ', time.time()-t0)

            t0 = time.time()
            if self.log.get_param('observation_date_key') == self.log.get_param('observation_time_key'):
                observation_time = ' '.join(fits_header[self.log.get_param('observation_date_key')].split('T'))
            else:
                observation_time = ' '.join([fits_header[self.log.get_param('observation_date_key')].split('T')[0],
                                             fits_header[self.log.get_param('observation_time_key')]])

            observation_time = plc.UTC(observation_time)
            if self.log.get_param('time_stamp') == 'exposure start':
                pass
            elif self.log.get_param('time_stamp') == 'mid-exposure':
                observation_time = observation_time - plc.DTime(seconds=0.5 * exp_time)
            elif self.log.get_param('time_stamp') == 'exposure end':
                observation_time = observation_time - plc.DTime(seconds=exp_time)
            else:
                raise RuntimeError('Not acceptable time stamp.')

            julian_date = observation_time.jd()
            airmass = self.observatory.airmass(self.target, observation_time)

            # write the new fits file
            # important to keep it like this for windows!

            time_in_file = observation_time.utc.isoformat()
            time_in_file = time_in_file.split('.')[0]
            time_in_file = time_in_file.replace('-', '_').replace(':', '_').replace('T', '_')

            new_name = '{0}{1}_{2}'.format(self.log.reduction_prefix, time_in_file, science_file.split(os.sep)[-1])
            primary = pf.PrimaryHDU()
            image = pf.CompImageHDU()
            image.data = np.array(data_frame, dtype=np.int32)
            image.header.set('BITPIX', fits_header['BITPIX'])
            image.header.set('NAXIS1', len(fits_data[0]))
            image.header.set('NAXIS2', len(fits_data))
            image.header.set('XBINNING', bin_fits)
            image.header.set('YBINNING', bin_fits)
            image.header.set('BZERO',  0)
            image.header.set('BSCALE', 1)
            image.header.set(self.log.hops_observatory_latitude_key, self.observatory.latitude.deg_coord())
            image.header.set(self.log.hops_observatory_longitude_key, self.observatory.longitude.deg())
            image.header.set(self.log.hops_target_ra_key, self.target.ra.deg())
            image.header.set(self.log.hops_target_dec_key, self.target.dec.deg_coord())
            image.header.set(self.log.hops_datetime_key, observation_time.iso())
            image.header.set(self.log.hops_exposure_key, exp_time)
            image.header.set(self.log.hops_filter_key, self.filter)
            image.header.set(self.log.time_key, julian_date)
            image.header.set(self.log.airmass_key, airmass)
            image.header.set(self.log.mean_key, mean)
            image.header.set(self.log.std_key, std)
            image.header.set(self.log.hops_saturation_key, saturation)
            image.header.set(self.log.psf_key, psf)
            image.header.set(self.log.skip_key, skip)
            image.header.set(self.log.align_x0_key, False)
            image.header.set(self.log.align_y0_key, False)
            image.header.set(self.log.align_u0_key, False)
            image.header.set(self.log.align_u0_key, False)
            fits_header[self.log.mean_key] = mean
            fits_header[self.log.std_key] = std

            plc.save_fits(pf.HDUList([primary, image]), os.path.join(self.log.reduction_directory, new_name))

            self.all_frames[new_name] = {
                self.log.mean_key: mean,
                self.log.std_key: std,
                self.log.psf_key: psf,
                self.log.time_key: julian_date,
                self.log.airmass_key: airmass,
                self.log.get_param('exposure_time_key'): exp_time,
                self.log.skip_key: skip,
                self.log.align_x0_key: False,
                self.log.align_y0_key: False,
                self.log.align_u0_key: False,
            }

            if timing:
                print('Saving: ', time.time()-t0)
                print('Total: ', time.time()-t00)

            self.progress_science.update()
            self.science_counter += 1

            if self.science_counter >= len(self.science_files):
                self.after(self.save)
            else:
                if self.progress_science_loop.get() or self.science_counter == 1:
                    self.progress_figure.load_fits(data_frame, fits_header, new_name)
                    self.progress_figure.draw()

                if len(self.jobs) > self.jobs_completed + 1:
                    self.fr_time += 10

                self.after(self.reduce_science, time=self.fr_time)

    def save(self):

        if self.exit:
            self.close()
        else:
            plc.save_dict(self.all_frames, self.log.all_frames)
            self.log.set_param('reduction_complete', True)
            self.log.set_param('reduction_version', self.log.version)
            self.log.save_local_log()
            self.log.set_param('proceed', True)
            self.close()
