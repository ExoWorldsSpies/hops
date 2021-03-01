
import os
import time
import numpy as np
import shutil
import warnings

from astropy.io import fits as pf

import hops.pylightcurve3 as plc

from hops.hops_tools.fits import find_fits_files, get_fits_data
from hops.application_windows import MainWindow


class ReductionWindow(MainWindow):

    def __init__(self, log):

        MainWindow.__init__(self, log, name='HOPS - Reduction', position=2)

        fits_name = find_fits_files(self.log.get_param('observation_files'))[0]
        fits = get_fits_data(fits_name)

        y_scale = (self.root.winfo_screenheight() - 500) / self.root.winfo_screenheight()

        self.progress_figure = self.FitsWindow(figsize=(0.5, y_scale, 10, 10, len(fits[0].data[0])/len(fits[0].data)))
        self.progress_figure.load_fits(fits[0], input_name=os.path.split(fits_name)[1])
        self.progress_bias = self.Progressbar(task="Creating master bias")
        self.progress_dark = self.Progressbar(task="Creating master dark")
        self.progress_flat = self.Progressbar(task="Creating master flat")
        self.progress_science = self.Progressbar(task="Reducing data and calculating statistics")
        self.progress_science_loop = self.CheckButton(text='Show all frames', initial=0)

        self.setup_window([
            [[self.progress_figure, 0, 2]],
            [[self.progress_bias, 0, 2]],
            [[self.progress_dark, 0, 2]],
            [[self.progress_flat, 0, 2]],
            [[self.progress_science, 0], [self.progress_science_loop, 1]],
            [[self.Button(text='STOP REDUCTION & RETURN TO MAIN MENU', command=self.trigger_exit), 0, 2]],
            []
        ])

        self.set_close_button_function(self.trigger_exit)

    # define functions

    def run_reduction(self):

        if self.log.get_param('reduction_complete'):
            if self.askyesno('Overwrite files', 'Reduction has been completed, do you want to run again?'):
                self.log.set_param('reduction_complete', False)
                self.log.save_local_log()
            else:
                self.log.set_param('proceed', True)
            self.show()

        if not self.log.get_param('reduction_complete'):

            if not os.path.isdir(self.log.reduction_directory):
                os.mkdir(self.log.reduction_directory)
            else:
                shutil.rmtree(self.log.reduction_directory)
                os.mkdir(self.log.reduction_directory)

            self.bias_files = find_fits_files(self.log.get_param('bias_files'))
            self.bias_frames = []
            self.bias_counter = 0

            self.dark_files = find_fits_files(self.log.get_param('dark_files'))
            self.dark_frames = []
            self.dark_counter = 0

            self.flat_files = find_fits_files(self.log.get_param('flat_files'))
            self.flat_frames = []
            self.flat_counter = 0

            self.science_files = find_fits_files(self.log.get_param('observation_files'))
            self.all_frames = {}
            self.science_counter = 0
            self.psf = 10

            self.after(self.get_bias)

        else:
            self.close()

    # reduction routines

    def get_bias(self):

        if self.exit or len(self.bias_files) == 0:
            self.after(self.get_master_bias)

        else:

            if self.bias_counter == 0:
                self.progress_bias.initiate(10 * len(self.bias_files))

            fits = get_fits_data(self.bias_files[self.bias_counter])

            self.bias_frames.append(np.ones_like(fits[0].data) * fits[0].data)

            self.progress_bias.update()
            self.bias_counter += 1

            if self.bias_counter >= len(self.bias_files):
                self.after(self.get_master_bias)
            else:
                self.after(self.get_bias)

    def get_master_bias(self):

        if self.exit:
            self.after(self.get_dark)

        else:

            if len(self.bias_frames) > 0:
                if self.log.get_param('master_bias_method') == 'median':
                    self.master_bias = np.nanmedian(self.bias_frames, 0)
                elif self.log.get_param('master_bias_method') == 'mean':
                    self.master_bias = np.nanmean(self.bias_frames, 0)
                else:
                    self.master_bias = np.nanmedian(self.bias_frames, 0)

                self.progress_bias.update(step=9 * len(self.bias_files))
            else:
                self.master_bias = 0.0

            print('Median Bias: ', round(np.nanmedian(self.master_bias), 3))

            self.after(self.get_dark)

    def get_dark(self):

        if self.exit or len(self.dark_files) == 0:
            self.after(self.get_master_dark)

        else:

            if self.dark_counter == 0:
                self.progress_dark.initiate(10 * len(self.dark_files))

            fits = get_fits_data(self.dark_files[self.dark_counter])

            dark_frame = np.ones_like(fits[0].data) * fits[0].data
            self.dark_frames.append((dark_frame - self.master_bias) / fits[0].header[self.log.get_param('exposure_time_key')])

            self.progress_dark.update()
            self.dark_counter += 1

            if self.dark_counter >= len(self.dark_files):
                self.after(self.get_master_dark)
            else:
                self.after(self.get_dark)

    def get_master_dark(self):

        if self.exit:
            self.after(self.get_flat)
        else:

            if len(self.dark_frames) > 0:
                if self.log.get_param('master_dark_method') == 'median':
                    self.master_dark = np.nanmedian(self.dark_frames, 0)
                elif self.log.get_param('master_dark_method') == 'mean':
                    self.master_dark = np.nanmean(self.dark_frames, 0)
                else:
                    self.master_dark = np.nanmedian(self.dark_frames, 0)

                self.progress_dark.update(step=9 * len(self.dark_files))
            else:
                self.master_dark = 0.0

            print('Median Dark: ', round(np.nanmedian(self.master_dark), 3))

            self.after(self.get_flat)

    def get_flat(self):

        if self.exit or len(self.flat_files) == 0:
            self.after(self.get_master_flat)

        else:

            if self.flat_counter == 0:
                self.progress_flat.initiate(10 * len(self.flat_files))

            fits = get_fits_data(self.flat_files[self.flat_counter])

            flat_frame = np.ones_like(fits[0].data) * fits[0].data
            self.flat_frames.append(
                flat_frame - self.master_bias - fits[0].header[self.log.get_param('exposure_time_key')] * self.master_dark)

            self.progress_flat.update()
            self.flat_counter += 1

            if self.flat_counter >= len(self.flat_files):
                self.after(self.get_master_flat)
            else:
                self.after(self.get_flat)

    def get_master_flat(self):

        if self.exit:
            self.after(self.reduce_science)

        else:

            if len(self.flat_frames) > 0:
                print('Median of each Flat: ', ' '.join([str(round(np.nanmedian(ff))) for ff in self.flat_frames]))
                if self.log.get_param('master_flat_method') == 'mean':
                    flat_frames = [ff / np.nanmean(ff) for ff in self.flat_frames]
                    self.master_flat = np.nanmean(flat_frames, 0)
                else:
                    flat_frames = [ff / np.nanmedian(ff) for ff in self.flat_frames]
                    self.master_flat = np.nanmedian(flat_frames, 0)
                print('Median Flat: ', round(np.nanmedian(self.master_flat), 3))
                self.master_flat = self.master_flat / np.nanmedian(self.master_flat)
                self.master_flat = np.where(self.master_flat == 0, 1, self.master_flat)

                self.progress_flat.update(step=9 * len(self.flat_files))
            else:
                self.master_flat = 1.0

            self.after(self.reduce_science)

    def reduce_science(self):

        # correct each observation_files file

        if self.exit:
            self.after(self.save)

        else:
            xx = time.time()

            if self.science_counter == 0:
                self.progress_science.initiate(len(self.science_files))

            science_file = self.science_files[self.science_counter]

            # correct it with master bias_files, master dark_files and master flat_files

            fits = get_fits_data(science_file)

            exp_time = fits[0].header[self.log.get_param('exposure_time_key')]
            data_frame = np.ones_like(fits[0].data) * fits[0].data
            data_frame = (data_frame - self.master_bias - exp_time * self.master_dark) / self.master_flat
            data_frame[np.where(np.isnan(data_frame))] = 0

            if self.log.get_param('bin_fits') > 1:
                data_frame = plc.bin_frame(data_frame, self.log.get_param('bin_fits'))

            try:
                distribution = plc.one_d_distribution(data_frame.flatten()[::int(200000.0/self.log.bin_to)],
                                                      gaussian_fit=True, mad_filter=5.0)
                mean = distribution[2]
                std = distribution[3]
            except:
                mean = np.median(data_frame)
                std = plc.mad(data_frame) * 1.5

            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                psf = plc.fast_psf_find(data_frame, mean, std, 0.95 * self.log.get_param('burn_limit'))
            if np.isnan(psf):
                psf = self.psf + 10
                skip = True
            else:
                self.psf = psf
                skip = False

            if self.log.get_param('observation_date_key') == self.log.get_param('observation_time_key'):
                observation_time = ' '.join(fits[0].header[self.log.get_param('observation_date_key')].split('T'))
            else:
                observation_time = ' '.join([fits[0].header[self.log.get_param('observation_date_key')].split('T')[0],
                                             fits[0].header[self.log.get_param('observation_time_key')]])

            observation_time = plc.UTC(observation_time)
            if self.log.get_param('time_stamp') == 'exposure start':
                julian_date = observation_time.jd
            elif self.log.get_param('time_stamp') == 'mid-exposure':
                julian_date = observation_time.jd - 0.5 * exp_time / 60.0 / 60.0 / 24.0
            elif self.log.get_param('time_stamp') == 'exposure end':
                julian_date = observation_time.jd - exp_time / 60.0 / 60.0 / 24.0
            else:
                raise RuntimeError('Not acceptable time stamp.')

            fits[0].header.set(self.log.mean_key, mean)
            fits[0].header.set(self.log.std_key, std)
            fits[0].header.set(self.log.psf_key, psf)
            fits[0].header.set(self.log.time_key, julian_date)

            # write the new fits file
            # important to keep it like this for windows!

            time_in_file = observation_time.utc.isoformat()
            time_in_file = time_in_file.split('.')[0]
            time_in_file = time_in_file.replace('-', '_').replace(':', '_').replace('T', '_')

            new_name = '{0}{1}_{2}'.format(self.log.reduction_prefix, time_in_file, science_file.split(os.sep)[-1])

            hdu = pf.ImageHDU(header=fits[0].header, data=np.array(data_frame, dtype=np.float32))

            plc.save_fits(pf.HDUList([pf.PrimaryHDU(), hdu]), os.path.join(self.log.reduction_directory, new_name))

            self.all_frames[new_name] = {self.log.mean_key: mean, self.log.std_key: std, self.log.psf_key: psf,
                                    self.log.time_key: julian_date,
                                    self.log.get_param('exposure_time_key'): fits[0].header[self.log.get_param('exposure_time_key')],
                                    self.log.skip_key: skip,
                                    self.log.align_x0_key: False,
                                    self.log.align_y0_key: False,
                                    self.log.align_u0_key: False}

            if self.progress_science_loop.get() or self.science_counter == 0:
                self.progress_figure.load_fits(hdu, new_name)

            # counter

            self.progress_science.update()
            self.science_counter += 1

            if self.science_counter >= len(self.science_files):
                self.after(self.save)
            else:
                self.after(self.reduce_science)

    def save(self):

        if self.exit:
            self.def_close()
        else:
            plc.save_dict(self.all_frames, self.log.all_frames)
            self.log.set_param('reduction_complete', True)
            self.log.set_param('reduction_version', self.log.version)
            self.log.save_local_log()
            self.log.set_param('proceed', True)
            self.def_close()
