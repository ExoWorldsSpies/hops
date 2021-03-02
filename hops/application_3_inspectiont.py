
import os
import numpy as np
import matplotlib.patches as mpatches

from matplotlib.backend_bases import MouseEvent as mpl_MouseEvent

import hops.pylightcurve3 as plc

from hops.hops_tools.fits import get_fits_data
from hops.application_windows import MainWindow


class InspectiontWindow(MainWindow):

    def __init__(self, log):

        MainWindow.__init__(self, log, name='HOPS - Inspection', position=2)

        # set variables, create and place widgets

        self.all_frames = plc.open_dict(self.log.all_frames)

        self.science_files = []
        for science_file in self.all_frames:
            self.science_files.append([self.all_frames[science_file][self.log.time_key], science_file])

        self.science_files.sort()
        self.science_files = [ff[1] for ff in self.science_files]

        self.science_keys = []
        self.science = []
        self.time_array = []
        self.sky_mean_array = []
        self.sky_std_array = []
        self.psf_array = []
        self.skip_array = []

        for science_file in self.science_files:
            self.science_keys.append(science_file)
            self.science.append(os.path.join(self.log.reduction_directory, science_file))
            self.time_array.append(self.all_frames[science_file][self.log.time_key])
            self.sky_mean_array.append(self.all_frames[science_file][self.log.mean_key] / self.all_frames[science_file][self.log.get_param('exposure_time_key')])
            self.sky_std_array.append(self.all_frames[science_file][self.log.std_key] / self.all_frames[science_file][self.log.get_param('exposure_time_key')])
            self.psf_array.append(self.all_frames[science_file][self.log.psf_key] * 2.355 / 2)
            self.skip_array.append(self.all_frames[science_file][self.log.skip_key])

        self.time_array = (np.array(self.time_array) - np.min(self.time_array)) * 24
        self.sky_mean_array = np.array(self.sky_mean_array)
        self.psf_array = np.array(self.psf_array)
        self.skip_array = np.array(self.skip_array)

        # main window

        y_scale = (self.root.winfo_screenheight() - 375) / self.root.winfo_screenheight()

        data = get_fits_data(self.science[0])[0].data

        self.fits_figure = self.FitsWindow(figsize=(0.5, y_scale, 20, 20, len(data[0])/len(data)), show_controls=True)
        self.fits_to_plot = np.argmin(self.time_array)
        self.fits_figure.load_fits(self.science[self.fits_to_plot])

        self.sky_threshold = self.Entry(value=self.log.get_param('sky_threshold'), instance=float, command=self.apply_thresholds)
        self.psf_threshold = self.Entry(value=self.log.get_param('psf_threshold'), instance=float, command=self.apply_thresholds)

        self.inspection_figure = self.FigureWindow(figsize=(0.5, 0.5, 10, 6, 1.1), show_nav=True)
        self.inspection_figure_ax1 = self.inspection_figure.figure.add_subplot(211)
        self.inspection_figure_ax2 = self.inspection_figure.figure.add_subplot(212)
        self.inspection_figure.figure.subplots_adjust(left=0.15, right=0.99, bottom=0.1, top=0.99)
        self.inspection_figure.figure.canvas.callbacks.connect('button_press_event', self.update_window)
        self.arrow1 = 0
        self.arrow2 = 0

        self.inspection_figure_ax1.plot(self.time_array, self.sky_mean_array, 'ko', ms=3)
        self.inspection_figure_ax2.plot(self.time_array, self.psf_array, 'ko', ms=3)

        self.inspection_figure_ax1.set_xlim(np.min(self.time_array) - 0.1, np.max(self.time_array) + 0.1)
        self.inspection_figure_ax2.set_xlim(np.min(self.time_array) - 0.1, np.max(self.time_array) + 0.1)

        self.yspil = self.inspection_figure.figure.get_size_inches()[1] * self.inspection_figure.figure.dpi * 0.55
        box1 = self.inspection_figure_ax1.get_window_extent()
        self.ax1_width = box1.width
        self.ax1_height = box1.height
        box2 = self.inspection_figure_ax2.get_window_extent()
        self.ax2_width = box2.width
        self.ax2_height = box2.height

        self.dbclick = False
        self.dbclick_xy = (0, 0)

        self.replot()

        self.setup_window([
            [[self.fits_figure, 0, 1, 10], [self.inspection_figure, 1, 4]],
            [],
            [[self.Label(text='On the time-sky or PSF-sky graph above'
                              '\n'
                              'double-click on a point to see the frame on the left panel.'
                              '\n'
                              'To mark this point as faulty, use the right double-click.'
                              '\n'
                              'To undo, use the right double-click again.'), 1, 4]],
            [],
            [[self.Label(text='Sky Threshold'), 1], [self.sky_threshold, 2],
             [self.Label(text='PSF Threshold'), 3], [self.psf_threshold, 4]],
            [],
            [[self.Button(text='RETURN TO MAIN MENU', command=self.close), 1, 4]],
            [[self.Button(text='SAVE OPTIONS & RETURN TO MAIN MENU', command=self.save_and_return), 1, 4]],
            [[self.Button(text='SAVE OPTIONS & PROCEED', command=self.save_and_proceed,
                          bg='green', highlightbackground='green'), 1, 4]],
            []
        ])

        self.replot()

    def update_window(self, event=None):

        if isinstance(event, mpl_MouseEvent):

            if event.inaxes is None:
                return None

        if (event.x, event.y) != self.dbclick_xy:
            self.dbclick_xy = (event.x, event.y)
            return None
        else:
            self.dbclick_xy = (0, 0)

            if event.y < self.yspil:

                time_norm = self.ax2_width / (self.inspection_figure_ax2.get_xlim()[1] - self.inspection_figure_ax2.get_xlim()[0])
                psf_norm = self.ax2_height / (self.inspection_figure_ax2.get_ylim()[1] - self.inspection_figure_ax2.get_ylim()[0])

                self.fits_to_plot = np.argmin(np.sqrt(
                    ((event.xdata - np.array(self.time_array)) * time_norm) ** 2 + ((event.ydata - np.array(self.psf_array)) * psf_norm)**2))

            else:

                time_norm = self.ax1_width / (self.inspection_figure_ax1.get_xlim()[1] - self.inspection_figure_ax1.get_xlim()[0])
                sky_norm = self.ax1_height / (self.inspection_figure_ax1.get_ylim()[1] - self.inspection_figure_ax1.get_ylim()[0])

                self.fits_to_plot = np.argmin(np.sqrt(
                    ((event.xdata - np.array(self.time_array)) * time_norm) ** 2 + ((event.ydata - np.array(self.sky_mean_array)) * sky_norm) ** 2))

            self.fits_figure.load_fits(self.science[self.fits_to_plot])

            if event.button == 3:

                if self.skip_array[self.fits_to_plot]:
                    self.skip_array[self.fits_to_plot] = False
                else:
                    self.skip_array[self.fits_to_plot] = True

            self.replot()

    def replot(self):

        ax1_xlim = self.inspection_figure_ax1.get_xlim()
        ax1_ylim = self.inspection_figure_ax1.get_ylim()
        ax2_xlim = self.inspection_figure_ax2.get_xlim()
        ax2_ylim = self.inspection_figure_ax2.get_ylim()

        del self.arrow1
        del self.arrow2
        self.inspection_figure_ax1.cla()
        self.inspection_figure_ax2.cla()

        self.inspection_figure_ax1.set_xlim(ax1_xlim)
        self.inspection_figure_ax1.set_ylim(ax1_ylim)
        self.inspection_figure_ax2.set_xlim(ax2_xlim)
        self.inspection_figure_ax2.set_ylim(ax2_ylim)

        if self.sky_threshold.get():
            self.inspection_figure_ax1.axhline(self.sky_threshold.get(), ls='--', c='r', lw=0.75)

        if self.psf_threshold.get():
            self.inspection_figure_ax2.axhline(self.psf_threshold.get(), ls='--', c='r', lw=0.75)

        self.inspection_figure_ax2.set_xlabel('Time (hours in observation)')

        self.inspection_figure_ax1.set_ylabel('Sky (counts/pix/s)')
        self.inspection_figure_ax2.set_ylabel('PSF max. HWHM (pix)')

        self.inspection_figure_ax1.get_xlim()
        self.inspection_figure_ax2.get_xlim()

        self.inspection_figure_ax1.plot(self.time_array, self.sky_mean_array, 'ko', ms=3)
        self.inspection_figure_ax2.plot(self.time_array, self.psf_array, 'ko', ms=3)

        self.inspection_figure_ax1.plot(self.time_array[np.where(self.skip_array)], self.sky_mean_array[np.where(self.skip_array)], 'ro', ms=3)
        self.inspection_figure_ax2.plot(self.time_array[np.where(self.skip_array)], self.psf_array[np.where(self.skip_array)], 'ro', ms=3)

        ymin1 = self.inspection_figure_ax1.get_ylim()[0]
        ymin2 = self.inspection_figure_ax2.get_ylim()[0]

        self.arrow1 = mpatches.Arrow(self.time_array[self.fits_to_plot], ymin1, 0, self.sky_mean_array[self.fits_to_plot] - ymin1,
                                     width=self.time_array[1] - self.time_array[0], fc='r')

        self.arrow2 = mpatches.Arrow(self.time_array[self.fits_to_plot], ymin2, 0, self.psf_array[self.fits_to_plot] - ymin2,
                                     width=self.time_array[1] - self.time_array[0], fc='r')

        self.inspection_figure_ax1.add_patch(self.arrow1)
        self.inspection_figure_ax2.add_patch(self.arrow2)

        self.inspection_figure.draw()

    def apply_thresholds(self, *event):

        try:

            for index in range(len(self.sky_mean_array)):

                keep = True

                if self.sky_threshold.get():
                    if self.sky_mean_array[index] > self.sky_threshold.get():
                        keep = False

                if self.psf_threshold.get():
                    if self.psf_array[index] > self.psf_threshold.get():
                        keep = False

                if not keep:
                    if not self.skip_array[index]:
                        self.skip_array[index] = True
                else:
                    if self.skip_array[index]:
                        self.skip_array[index] = False

            self.replot()

        except:
            pass

    def save(self):

        self.log.set_param('sky_threshold', self.sky_threshold.get())
        self.log.set_param('psf_threshold', self.psf_threshold.get())

        for science_file in range(len(self.science_files)):
            self.all_frames[self.science_files[science_file]][self.log.skip_key] = self.skip_array[science_file]

        plc.save_dict(self.all_frames, self.log.all_frames)

        self.log.save_local_log()

    def save_and_return(self):

        self.save()
        self.log.set_param('proceed', False)
        self.close()

    def save_and_proceed(self):

        self.save()
        self.log.set_param('proceed', True)
        self.close()
