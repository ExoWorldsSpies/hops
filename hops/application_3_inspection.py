
import os
import numpy as np
import hops.pylightcurve41 as plc
import matplotlib.patches as mpatches
from astropy.coordinates import SkyCoord

from matplotlib.backend_bases import MouseEvent as mpl_MouseEvent

from hops.hops_tools.fits import get_fits_data_and_header
from hops.application_windows import MainWindow


class InspectiontWindow(MainWindow):

    def __init__(self, log):

        MainWindow.__init__(self, log, name='HOPS - Inspection', position=2)

        ##############################################################################
        # set up time series

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
        self.airmass_array = []
        self.psf_array = []
        self.skip_array = []

        for science_file in self.science_files:
            self.science_keys.append(science_file)
            self.science.append(os.path.join(self.log.reduction_directory, science_file))
            self.time_array.append(self.all_frames[science_file][self.log.time_key])
            self.sky_mean_array.append(self.all_frames[science_file][self.log.mean_key] / self.all_frames[science_file][self.log.get_param('exposure_time_key')])
            self.sky_std_array.append(self.all_frames[science_file][self.log.std_key] / self.all_frames[science_file][self.log.get_param('exposure_time_key')])
            self.psf_array.append(self.all_frames[science_file][self.log.psf_key] * 2.355 / 2)
            try:
                self.airmass_array.append(self.all_frames[science_file][self.log.airmass_key])
            except:
                self.airmass_array.append(1)
            self.skip_array.append(self.all_frames[science_file][self.log.skip_key])

        self.time_array = (np.array(self.time_array) - np.min(self.time_array)) * 24
        self.sky_mean_array = np.array(self.sky_mean_array)
        self.airmass_array = np.array(self.airmass_array)
        self.psf_array = np.array(self.psf_array)
        self.skip_array = np.array(self.skip_array)
        ##############################################################################

        ##############################################################################
        # set up plot of first frame

        test_fits_name = self.science[0]
        test_fits_data, test_fits_header = get_fits_data_and_header(test_fits_name)

        self.fits_figure = self.FitsWindow(fits_data=test_fits_data, fits_header=test_fits_header,
                                           input_name=self.science_files[0],
                                           show_controls=True, show_axes=True,
                                           subplots_adjust=(0.07, 0.99, 0.05, 0.99))
        self.fits_to_plot = np.argmin(self.time_array)
        self.replot_fits()
        ##############################################################################

        ##############################################################################
        # set up inspection plot

        # self.sky_threshold = self.Entry(value=self.log.get_param('sky_threshold'), instance=float,
        #                                 command=self.apply_thresholds)
        # self.psf_threshold = self.Entry(value=self.log.get_param('psf_threshold'), instance=float,
        #                                 command=self.apply_thresholds)

        self.inspection_figure = self.FigureWindow(figsize=(1, 1), show_nav=True)

        self.inspection_figure_ax1 = self.inspection_figure.figure.add_subplot(211)
        self.inspection_figure_ax1.plot(self.time_array, self.sky_mean_array, 'ko', ms=3)
        self.excluded1,  = self.inspection_figure_ax1.plot(self.time_array[np.where(self.skip_array)], self.sky_mean_array[np.where(self.skip_array)], 'ro', ms=3)
        self.inspection_figure_ax1.set_xlim(np.min(self.time_array) - 0.1, np.max(self.time_array) + 0.1)
        self.inspection_figure_ax1.set_ylabel('Sky (counts/pix/s)')
        self.inspection_figure_ax1_airmass = self.inspection_figure_ax1.twinx()
        self.inspection_figure_ax1_airmass.plot(self.time_array, self.airmass_array, 'k--', lw=1)
        self.inspection_figure_ax1_airmass.set_ylim(1, max(self.airmass_array))
        self.inspection_figure_ax1_airmass.set_ylabel('Airmass')
        self.inspection_figure_ax1_airmass.set_zorder(0)
        self.inspection_figure_ax1.set_zorder(1)
        self.inspection_figure_ax1.patch.set_alpha(0.0)
        if self.sky_mean_array[self.fits_to_plot] - min(self.sky_mean_array) > max(self.sky_mean_array) - self.sky_mean_array[self.fits_to_plot]:
            arrow1 = mpatches.Arrow(self.time_array[self.fits_to_plot],
                                    min(self.sky_mean_array),
                                    0, self.sky_mean_array[self.fits_to_plot] - min(self.sky_mean_array),
                                    width=0.02 * (max(self.inspection_figure_ax1.get_xlim()) - min(self.inspection_figure_ax1.get_xlim())), fc='r')
        else:
            arrow1 = mpatches.Arrow(self.time_array[self.fits_to_plot],
                                    max(self.sky_mean_array),
                                    0, self.sky_mean_array[self.fits_to_plot] - max(self.sky_mean_array),
                                    width=0.02 * (max(self.inspection_figure_ax1.get_xlim()) - min(self.inspection_figure_ax1.get_xlim())), fc='r')
        self.arrow1 = self.inspection_figure_ax1_airmass.add_patch(arrow1)

        self.inspection_figure_ax2 = self.inspection_figure.figure.add_subplot(212)
        self.inspection_figure_ax2.plot(self.time_array, self.psf_array, 'ko', ms=3)
        self.excluded2,  = self.inspection_figure_ax2.plot(self.time_array[np.where(self.skip_array)], self.psf_array[np.where(self.skip_array)], 'ro', ms=3)
        self.inspection_figure_ax2.set_xlabel('Time (hours in observation)')
        self.inspection_figure_ax2.set_ylabel('PSF max.\nHWHM (pix)')
        self.inspection_figure_ax2.set_xlim(np.min(self.time_array) - 0.1, np.max(self.time_array) + 0.1)
        self.inspection_figure_ax2.xaxis.tick_top()
        self.inspection_figure_ax2.xaxis.set_label_position('top')
        self.inspection_figure_ax2_airmass = self.inspection_figure_ax2.twinx()
        self.inspection_figure_ax2_airmass.plot(self.time_array, self.airmass_array, 'k--', lw=1)
        self.inspection_figure_ax2_airmass.set_ylim(1, max(self.airmass_array))
        self.inspection_figure_ax2_airmass.set_ylabel('Airmass')
        self.inspection_figure_ax2_airmass.set_zorder(0)
        self.inspection_figure_ax2.set_zorder(1)
        self.inspection_figure_ax2.patch.set_alpha(0.0)
        if self.psf_array[self.fits_to_plot] - min(self.psf_array) > max(self.psf_array) - self.psf_array[self.fits_to_plot]:
            arrow2 = mpatches.Arrow(self.time_array[self.fits_to_plot],
                                    min(self.psf_array),
                                    0, self.psf_array[self.fits_to_plot] - min(self.psf_array),
                                    width=0.02 * (max(self.time_array) - min(self.time_array)), fc='r')
        else:
            arrow2 = mpatches.Arrow(self.time_array[self.fits_to_plot],
                                    max(self.psf_array),
                                    0, self.psf_array[self.fits_to_plot] - max(self.psf_array),
                                    width=0.02 * (max(self.time_array) - min(self.time_array)), fc='r')
        self.arrow2 = self.inspection_figure_ax2.add_patch(arrow2)

        self.inspection_figure.figure.canvas.callbacks.connect('button_press_event', self.select_frame)
        self.inspection_figure.figure.subplots_adjust(left=0.2, right=0.85, bottom=0.03, top=0.97, hspace=0.5, wspace=0.5)

        box1 = self.inspection_figure_ax1.get_window_extent()
        self.ax1_width = box1.width
        self.ax1_height = box1.height
        box2 = self.inspection_figure_ax2.get_window_extent()
        self.ax2_width = box2.width
        self.ax2_height = box2.height

        self.replot_inspection_figure()
        ##############################################################################

        ##############################################################################
        # place widgets

        self.setup_window([
            [[self.fits_figure, 0, 1, 11], [self.inspection_figure, 1, 4]],
            [[self.Label(text='Double-click on a point to see the frame on the left panel.'), 1, 4]],
            [],
            [[self.Label(text='SELECTED FRAME'), 1, 2],
             [self.Button(text='EXCLUDE', command=[self.exclude_selected, self.replot_inspection_figure]), 3, 1],
             [self.Button(text='INCLUDE', command=[self.include_selected, self.replot_inspection_figure]), 4, 1]],
            [[self.Label(text='ALL FRAMES IN SKY GRAPH'), 1, 2],
             [self.Button(text='EXCLUDE', command=[self.exclude_in_sky, self.replot_inspection_figure]), 3, 1],
             [self.Button(text='INCLUDE', command=[self.include_in_sky, self.replot_inspection_figure]), 4, 1]],
            [[self.Label(text='ALL FRAMES IN PSF GRAPH'), 1, 2],
             [self.Button(text='EXCLUDE', command=[self.exclude_in_psf, self.replot_inspection_figure]), 3, 1],
             [self.Button(text='INCLUDE', command=[self.include_in_psf, self.replot_inspection_figure]), 4, 1]],
            # [[self.Label(text='Sky Threshold'), 1], [self.sky_threshold, 2],
            #  [self.Label(text='PSF Threshold'), 3], [self.psf_threshold, 4]],
            [],
            [[self.Button(text='RETURN TO \n MAIN MENU', command=self.close,
                          bg='red', highlightbackground='red'), 1, 2],
             [self.Button(text='SAVE OPTIONS & \n PROCEED', command=self.save_and_proceed,
                          bg='green', highlightbackground='green'), 3, 2]],
            []
        ])
        ##############################################################################

        ##############################################################################
        # adjust figure sizes

        self.fits_figure.adjust_size()
        self.inspection_figure.canvas.get_tk_widget().config(
            width=0.9 * self.root.winfo_screenwidth() - self.fits_figure.canvas.get_tk_widget().winfo_reqwidth())
        self.inspection_figure.canvas.get_tk_widget().config(
            height=0.75 * self.fits_figure.canvas.get_tk_widget().winfo_reqheight())
        ##############################################################################

    def first_image_warning(self):

        self.showinfo('Inspection', 'IMPORTANT! The first frame must be of good quality, '
                                    'please exclude it if it has low S/N or startrails.')

    def select_frame(self, event=None):

        if isinstance(event, mpl_MouseEvent):
            if event.inaxes:
                if event.dblclick:
                    if event.y < 0.5 * self.inspection_figure.canvas.get_tk_widget().winfo_reqheight():

                        time_norm = self.ax2_width / (self.inspection_figure_ax2.get_xlim()[1] - self.inspection_figure_ax2.get_xlim()[0])
                        psf_norm = self.ax2_height / (self.inspection_figure_ax2.get_ylim()[1] - self.inspection_figure_ax2.get_ylim()[0])

                        self.fits_to_plot = np.argmin(np.sqrt(
                            ((event.xdata - np.array(self.time_array)) * time_norm) ** 2 + ((event.ydata - np.array(self.psf_array)) * psf_norm)**2))

                    else:

                        time_norm = self.ax1_width / (self.inspection_figure_ax1.get_xlim()[1] - self.inspection_figure_ax1.get_xlim()[0])
                        sky_norm = self.ax1_height / (self.inspection_figure_ax1.get_ylim()[1] - self.inspection_figure_ax1.get_ylim()[0])

                        self.fits_to_plot = np.argmin(np.sqrt(
                            ((event.xdata - np.array(self.time_array)) * time_norm) ** 2 + ((event.ydata - np.array(self.sky_mean_array)) * sky_norm) ** 2))

                    self.replot_fits()
                    self.replot_inspection_figure()

    def replot_fits(self):

        fits_data_to_plot, fits_header_to_plot = get_fits_data_and_header(self.science[self.fits_to_plot])
        self.fits_figure.load_fits(fits_data_to_plot, fits_header_to_plot, self.science[self.fits_to_plot],
                                   input_options=self.fits_figure.get_fov_options())

    def replot_inspection_figure(self):

        self.excluded1.set_xdata(self.time_array[np.where(self.skip_array)])
        self.excluded1.set_ydata(self.sky_mean_array[np.where(self.skip_array)])

        self.excluded2.set_xdata(self.time_array[np.where(self.skip_array)])
        self.excluded2.set_ydata(self.psf_array[np.where(self.skip_array)])

        self.arrow1.remove()
        if self.sky_mean_array[self.fits_to_plot] - min(self.sky_mean_array) > max(self.sky_mean_array) - self.sky_mean_array[self.fits_to_plot]:
            arrow1 = mpatches.Arrow(self.time_array[self.fits_to_plot],
                                    min(self.sky_mean_array),
                                    0, self.sky_mean_array[self.fits_to_plot] - min(self.sky_mean_array),
                                    width=0.02 * (max(self.time_array) - min(self.time_array)), fc='r')
        else:
            arrow1 = mpatches.Arrow(self.time_array[self.fits_to_plot],
                                    max(self.sky_mean_array),
                                    0, self.sky_mean_array[self.fits_to_plot] - max(self.sky_mean_array),
                                    width=0.02 * (max(self.time_array) - min(self.time_array)), fc='r')
        self.arrow1 = self.inspection_figure_ax1.add_patch(arrow1)

        self.arrow2.remove()
        if self.psf_array[self.fits_to_plot] - min(self.psf_array) > max(self.psf_array) - self.psf_array[self.fits_to_plot]:
            arrow2 = mpatches.Arrow(self.time_array[self.fits_to_plot],
                                    min(self.psf_array),
                                    0, self.psf_array[self.fits_to_plot] - min(self.psf_array),
                                    width=0.02 * (max(self.time_array) - min(self.time_array)), fc='r')
        else:
            arrow2 = mpatches.Arrow(self.time_array[self.fits_to_plot],
                                    max(self.psf_array),
                                    0, self.psf_array[self.fits_to_plot] - max(self.psf_array),
                                    width=0.02 * (max(self.time_array) - min(self.time_array)), fc='r')
        self.arrow2 = self.inspection_figure_ax2.add_patch(arrow2)

        self.inspection_figure.draw()

    def exclude_selected(self):

        self.skip_array[self.fits_to_plot] = True

    def exclude_in_sky(self):

        xlim1, xlim2 = self.inspection_figure_ax1.get_xlim()
        ylim1, ylim2 = self.inspection_figure_ax1.get_ylim()

        self.skip_array[np.where((self.time_array >= xlim1) *
                                 (self.time_array <= xlim2) *
                                 (self.sky_mean_array >= ylim1) *
                                 (self.sky_mean_array <= ylim2)
                                 )] = True

    def exclude_in_psf(self):

        xlim1, xlim2 = self.inspection_figure_ax2.get_xlim()
        ylim1, ylim2 = self.inspection_figure_ax2.get_ylim()

        self.skip_array[np.where((self.time_array >= xlim1) *
                                 (self.time_array <= xlim2) *
                                 (self.psf_array >= ylim1) *
                                 (self.psf_array <= ylim2)
                                 )] = True

    def include_selected(self):

        self.skip_array[self.fits_to_plot] = False

    def include_in_sky(self):

        xlim1, xlim2 = self.inspection_figure_ax1.get_xlim()
        ylim1, ylim2 = self.inspection_figure_ax1.get_ylim()

        self.skip_array[np.where((self.time_array >= xlim1) *
                                 (self.time_array <= xlim2) *
                                 (self.sky_mean_array >= ylim1) *
                                 (self.sky_mean_array <= ylim2)
                                 )] = False

    def include_in_psf(self):

        xlim1, xlim2 = self.inspection_figure_ax2.get_xlim()
        ylim1, ylim2 = self.inspection_figure_ax2.get_ylim()

        self.skip_array[np.where((self.time_array >= xlim1) *
                                 (self.time_array <= xlim2) *
                                 (self.psf_array >= ylim1) *
                                 (self.psf_array <= ylim2)
                                 )] = False

    def save(self):

        # self.log.set_param('sky_threshold', self.sky_threshold.get())
        # self.log.set_param('psf_threshold', self.psf_threshold.get())

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
