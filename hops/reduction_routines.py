
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


def reduction():

    print('Reduction...')

    if log.read_local_log('pipeline', 'reduction_complete'):
        if not askyesno('Overwrite files', 'Reduction has been completed, do you want to run again?'):
            return 0

    log.write_local_log('pipeline', False, 'reduction_complete')
    log.write_local_log('pipeline', False, 'alignment_complete')

    # get variables

    observation_files = log.read_local_log('pipeline', 'observation_files')
    reduction_directory = log.read_local_log('pipeline', 'reduction_directory')
    reduction_prefix = log.read_local_log('pipeline', 'reduction_prefix')
    exposure_time_key = log.read_local_log('pipeline_keywords', 'exposure_time_key')
    mean_key = log.read_local_log('pipeline_keywords', 'mean_key')
    std_key = log.read_local_log('pipeline_keywords', 'std_key')
    observation_date_key = log.read_local_log('pipeline_keywords', 'observation_date_key')
    observation_time_key = log.read_local_log('pipeline_keywords', 'observation_time_key')
    frame_low_std = log.read_local_log('windows', 'frame_low_std')
    frame_upper_std = log.read_local_log('windows', 'frame_upper_std')
    bias_files = log.read_local_log('reduction', 'bias_files')
    dark_files = log.read_local_log('reduction', 'dark_files')
    flat_files = log.read_local_log('reduction', 'flat_files')
    bin_fits = int(log.read_local_log('reduction', 'bin_fits'))
    bin_to = int(log.read_local_log('reduction', 'bin_to'))
    master_bias_method = log.read_local_log('reduction', 'master_bias_method')
    master_dark_method = log.read_local_log('reduction', 'master_dark_method')
    master_flat_method = log.read_local_log('reduction', 'master_flat_method')

    # check if reduction directory exists

    if os.path.isdir(reduction_directory):
        shutil.rmtree(reduction_directory)

    os.mkdir(reduction_directory)

    observation_files = find_fits_files(observation_files)

    def get_master_bias():

        bias_frames = []
        percent = 0
        lt0 = time.time()
        if len(bias_files) > 0:
            bb = find_fits_files(bias_files)
            for counter, bias_file in enumerate(bb):

                if show_progress.exit:
                    return None

                fits = pf.open(bias_file, memmap=False)
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

                bias_frames.append(np.ones_like(fits[0].data) * fits[0].data)

                new_percent = round(100 * (counter + 1) / float(len(bb)), 1)
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

                    progress_bar_1.update()
                    percent_label_1.update()

        if len(bias_frames) > 0:
            if master_bias_method == 'median':
                master_bias = np.median(bias_frames, 0)
            elif master_bias_method == 'mean':
                master_bias = np.mean(bias_frames, 0)
            else:
                master_bias = np.median(bias_frames, 0)
        else:
            master_bias = 0.0
            progress_bar_1['value'] = 0
            percent_label_1.configure(text='No bias frames')
            progress_bar_1.update()
            percent_label_1.update()

        print('Median Bias: ', round(np.median(master_bias), 3))

        return master_bias

    def get_master_dark(master_bias):

        dark_frames = []
        percent = 0
        lt0 = time.time()
        if len(str(dark_files)) > 0:
            dd = find_fits_files(dark_files)
            for counter, dark_file in enumerate(dd):

                if show_progress.exit:
                    return None

                fits = pf.open(dark_file, memmap=False)
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

                dark_frame = np.ones_like(fits[0].data) * fits[0].data
                dark_frames.append((dark_frame - master_bias) / fits[0].header[exposure_time_key])

                new_percent = round(100 * (counter + 1) / float(len(dd)), 1)
                if new_percent != percent:
                    lt1 = time.time()
                    rm_time = (100 - new_percent) * (lt1 - lt0) / new_percent
                    hours = rm_time / 3600.0
                    minutes = (hours - int(hours)) * 60
                    seconds = (minutes - int(minutes)) * 60

                    progress_bar_2['value'] = new_percent
                    percent_label_2.configure(text='{0} % ({1}h {2}m {3}s left)'.format(new_percent, int(hours),
                                                                                        int(minutes), int(seconds)))

                    percent = new_percent

                    progress_bar_2.update()
                    percent_label_2.update()

        if len(dark_frames) > 0:
            if master_dark_method == 'median':
                master_dark = np.median(dark_frames, 0)
            elif master_dark_method == 'mean':
                master_dark = np.mean(dark_frames, 0)
            else:
                master_dark = np.median(dark_frames, 0)
        else:
            master_dark = 0.0
            progress_bar_2['value'] = 0
            percent_label_2.configure(text='No dark frames')
            progress_bar_2.update()
            percent_label_2.update()

        print('Median Dark: ', round(np.median(master_dark), 3))

        return master_dark

    def get_master_flat(master_bias, master_dark):

        flat_frames = []
        percent = 0
        lt0 = time.time()
        if len(str(flat_files)) > 0:
            ff = find_fits_files(flat_files)
            for counter, flat_file in enumerate(ff):

                if show_progress.exit:
                    return None

                fits = pf.open(flat_file, memmap=False)
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

                flat_frame = np.ones_like(fits[0].data) * fits[0].data
                flat_frames.append(flat_frame - master_bias - fits[0].header[exposure_time_key] * master_dark)

                new_percent = round(100 * (counter + 1) / float(len(ff)), 1)
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

                    progress_bar_3.update()
                    percent_label_3.update()

        print('Median of each Flat: ', ' '.join([str(round(np.median(ff))) for ff in flat_frames]))
        if len(flat_frames) > 0:
            if master_flat_method == 'mean':
                flat_frames = [ff / np.mean(ff) for ff in flat_frames]
                master_flat = np.mean(flat_frames, 0)
            else:
                flat_frames = [ff / np.median(ff) for ff in flat_frames]
                master_flat = np.median(flat_frames, 0)
            print('Median Flat: ', round(np.median(master_flat), 3))
            master_flat = master_flat / np.median(master_flat)
            master_flat = np.where(master_flat == 0, 1, master_flat)
        else:
            master_flat = 1.0
            progress_bar_3['value'] = 0
            percent_label_3.configure(text='No flat frames')
            progress_bar_3.update()
            percent_label_3.update()

        return master_flat

    def reduce(master_bias, master_dark, master_flat):

        # correct each observation_files file

        percent = 0
        lt0 = time.time()

        testx = []
        testy = []
        testz = []

        for counter, science_file in enumerate(observation_files):

            if show_progress.exit:
                return None, None, None

            label_4.configure(text='Reducing data and calculating statistics: {0}'.format(os.path.split(science_file)[1]))
            label_4.update()
            # correct it with master bias_files, master dark_files and master flat_files

            fits_file = pf.open(science_file, memmap=False)

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
                fits = [fits_file[sci_id]]

            fits_file.close()

            data_frame = np.ones_like(fits[0].data) * fits[0].data
            data_frame = (data_frame - master_bias - fits[0].header[exposure_time_key] * master_dark) / master_flat

            if bin_fits > 1:
                data_frame = plc.bin_frame(data_frame, bin_fits)

            try:
                distribution = plc.one_d_distribution(data_frame.flatten()[::int(200000.0 / bin_to)],
                                                      gaussian_fit=True, mad_filter=5.0)
                mean = distribution[2]
                std = distribution[3]
            except:
                mean = np.median(data_frame)
                std = plc.mad(data_frame) * 1.5

            if observation_date_key == observation_time_key:
                observation_time = ' '.join(fits[0].header[observation_date_key].split('T'))
            else:
                observation_time = ' '.join([fits[0].header[observation_date_key].split('T')[0],
                                             fits[0].header[observation_time_key]])

            observation_time = plc.UTC(observation_time)

            testx.append(observation_time.jd)
            testy.append(mean / fits[0].header[exposure_time_key])
            testz.append(std)

            fits[0].header.set(mean_key, mean)
            fits[0].header.set(std_key, std)

            # write the new fits file
            # important to keep it like this for windows!
            time_in_file = observation_time.utc.isoformat()
            time_in_file = time_in_file.split('.')[0]
            time_in_file = time_in_file.replace('-', '_').replace(':', '_').replace('T', '_')

            hdu = pf.ImageHDU(header=fits[0].header,
                              data=np.nan_to_num(np.array(data_frame, dtype=np.float32),
                                                 nan=0.0, posinf=0.0, neginf=0.0))

            plc.save_fits(pf.HDUList([pf.PrimaryHDU(), hdu]), '{0}{1}{2}{3}_{4}'.format(
                reduction_directory, os.sep, reduction_prefix, time_in_file, science_file.split(os.sep)[-1]))

            if counter == 0:
                ax.cla()
                ax.imshow(data_frame[::2, ::2], origin='lower', cmap=cm.Greys_r,
                          vmin=fits[0].header[mean_key] + frame_low_std * fits[0].header[std_key],
                          vmax=fits[0].header[mean_key] + frame_upper_std * fits[0].header[std_key])
                ax.axis('off')

                canvas.draw()

            # counter

            new_percent = round(100 * (counter + 1) / float(len(observation_files)), 1)
            if new_percent != percent:
                lt1 = time.time()
                rm_time = (100 - new_percent) * (lt1 - lt0) / new_percent
                hours = rm_time / 3600.0
                minutes = (hours - int(hours)) * 60
                seconds = (minutes - int(minutes)) * 60

                progress_bar_4['value'] = new_percent
                percent_label_4.configure(text='{0} % ({1}h {2}m {3}s left)'.format(new_percent, int(hours),
                                                                                    int(minutes), int(seconds)))

                percent = new_percent

            show_progress.update()

        return testx, testy, testz

    def show_sky(testx, testy, testz):

        root = ProgressWindow('HOPS - Alignment', 0, 0, 5)

        testx = np.array(np.array(testx) - testx[0]) * 24.0 * 60.0

        reduction_trash_directory = log.read_local_log('pipeline', 'reduction_trash_directory')
        trash = log.read_local_log('pipeline', 'trash')
        if not trash:
            list_to_remove = []
        else:
            list_to_remove = list(np.int_(trash))

        f = Figure()
        ax = f.add_subplot(2, 1, 1)
        ax2 = f.add_subplot(2, 2, 3)
        ax3 = f.add_subplot(2, 2, 4)

        f.patch.set_facecolor('white')
        canvas = root.FigureCanvasTkAgg(f)
        canvas.get_tk_widget().pack()
        root.NavigationToolbar2Tk(canvas)

        ax.plot(testx, testy, 'ko', ms=3)
        for ii in list_to_remove:
            ax.plot(testx[ii], testy[ii], 'ro', ms=3)

        time_dt = np.median(np.array(testx[1:]) - np.array(testx[:-1]))
        arrow = mpatches.Arrow(testx[0], 0, 0, testy[0], width=time_dt, fc='r')
        ax.add_patch(arrow)
        ax.set_xlabel('Time (minutes from observation start)')
        ax.set_ylabel('Sky (counts / pixel / second)')
        ax.tick_params(top='on', bottom='off', labeltop='on', labelbottom='off')
        ax.xaxis.set_label_position('top')

        science = find_fits_files(os.path.join(reduction_directory, '*'))

        fits = pf.open(science[0], memmap=False)

        ax2.imshow(fits[1].data, origin='lower', cmap=cm.Greys_r,
                   vmin=fits[1].header[mean_key] + frame_low_std * fits[1].header[std_key],
                   vmax=fits[1].header[mean_key] + frame_upper_std * fits[1].header[std_key])
        ax2.axis('off')

        fits.close()

        ax3.text(-100105, -100100, 'Select faulty frames', va='center', ha='center')
        ax3.text(-100111, -100101, '>On the time-sky graph above\n'
                                   'double-click on a point to see\n'
                                   'the frame on the left panel.\n'
                                   '>To mark this point as faulty,\n'
                                   'use the right double-click.\n'
                                   '>To undo, use the right\n'
                                   'double-click again.', va='top')
        ax3.text(-100105, -100109, 'RUN ALIGNMENT', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 5},
                 va='center', ha='center')
        ax3.set_xlim(-100110, -100100)
        ax3.set_ylim(-100110, -100100)
        ax3.axis('off')

        def update_window_show(event):

            if event.inaxes is not None:

                if testx[0] < event.xdata < testx[-1]:

                    plot_hjd = np.argmin(np.abs(event.xdata - np.array(testx)))

                    if event.dblclick:

                        del ax.patches[0]
                        arrow2 = mpatches.Arrow(testx[plot_hjd], 0, 0, testy[plot_hjd], width=time_dt, fc='r')
                        ax.add_patch(arrow2)

                        ax2.cla()
                        fits2 = pf.open(science[plot_hjd], memmap=False)
                        ax2.imshow(fits2[1].data, origin='lower', cmap=cm.Greys_r,
                                   vmin=fits2[1].header[mean_key] + frame_low_std * fits2[1].header[std_key],
                                   vmax=fits2[1].header[mean_key] + frame_upper_std * fits2[1].header[std_key])
                        ax2.axis('off')

                        if event.button == 3:

                            pltxlim1, pltxlim2 = ax.get_xlim()
                            pltylim1, pltylim2 = ax.get_ylim()

                            if plot_hjd not in list_to_remove:
                                ax.plot(testx[plot_hjd], testy[plot_hjd], 'ro', ms=3)
                                list_to_remove.append(plot_hjd)
                            else:
                                ax.plot(testx[plot_hjd], testy[plot_hjd], 'ko', ms=3)
                                list_to_remove.remove(plot_hjd)

                            ax.set_xlim(pltxlim1, pltxlim2)
                            ax.set_ylim(pltylim1, pltylim2)

                    canvas.draw()

                else:
                    if -100110 < event.ydata < -100108:
                        if -100108 < event.xdata < -100102:
                            ax3.cla()
                            ax3.text(-100105, -100100, 'Select faulty frames', va='center', ha='center')
                            ax3.text(-100111, -100101, '>On the time-sky graph above\n'
                                                       'double-click on a point to see\n'
                                                       'the frame on the left panel.\n'
                                                       '>To mark this point as faulty,\n'
                                                       'use the right double-click.\n'
                                                       '>To undo, use the right\n'
                                                       'double-click again.', va='top')
                            ax3.text(-100105, -100109, 'RUN ALIGNMENT', color='w',
                                     bbox={'facecolor': 'blue', 'alpha': 0.5, 'pad': 5},
                                     va='center', ha='center')
                            ax3.set_xlim(-100110, -100100)
                            ax3.set_ylim(-100110, -100100)
                            ax3.axis('off')
                            canvas.draw()
                            time.sleep(0.5)
                            root.exit = True

        f.canvas.callbacks.connect('button_press_event', update_window_show)
        # f.canvas.callbacks.connect('button_press_event', run_alignment)

        root.show()

        while not root.exit:
            root.update()

        root.close()
        if not os.path.isdir(os.path.join(reduction_directory, reduction_trash_directory)):
            os.mkdir(os.path.join(reduction_directory, reduction_trash_directory))

        for iii in list_to_remove:
            shutil.move(science[iii], os.path.join(reduction_directory, reduction_trash_directory))

        log.write_local_log('pipeline', list(map(str, list_to_remove)), 'trash')

    def run():

        master_bias = get_master_bias()
        if not show_progress.exit:
            master_dark = get_master_dark(master_bias)
            if not show_progress.exit:
                master_flat = get_master_flat(master_bias, master_dark)
                if not show_progress.exit:
                    x, y, z = reduce(master_bias, master_dark, master_flat)
                    if not show_progress.exit:
                        show_sky(x, y, z)
                        if not show_progress.exit:
                            log.write_local_log('pipeline', True, 'reduction_complete')

        show_progress.close()

    # progress window

    show_progress = ProgressWindow('HOPS - Reduction', 0, 0, 5)

    f = Figure()
    f.patch.set_facecolor('white')
    ax = f.add_subplot(111)
    f.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    ax.axis('off')
    canvas = show_progress.FigureCanvasTkAgg(f)
    canvas.get_tk_widget().pack()

    fits_file = pf.open(observation_files[0], memmap=False)
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
        fits = [fits_file[sci_id]]
    mean = np.median(fits[0].data)
    std = plc.mad(fits[0].data) * 1.5
    ax.imshow(fits[0].data, origin='lower', cmap=cm.Greys_r,
              vmin=mean + frame_low_std * std,
              vmax=mean + frame_upper_std * std)
    fits_file.close()

    frame1 = show_progress.Frame()
    frame1.pack()

    label_1 = Label(frame1, text="Creating master bias")
    progress_bar_1 = ttk.Progressbar(frame1, orient=HORIZONTAL,
                                     length=300, maximum=100, mode='determinate', value=0)
    percent_label_1 = Label(frame1, text='0.0 %')

    label_2 = Label(frame1, text="Creating master dark")
    progress_bar_2 = ttk.Progressbar(frame1, orient=HORIZONTAL,
                                     length=300, maximum=100, mode='determinate', value=0)
    percent_label_2 = Label(frame1, text='0.0 %')

    label_3 = Label(frame1, text="Creating master flat")
    progress_bar_3 = ttk.Progressbar(frame1, orient=HORIZONTAL,
                                     length=300, maximum=100, mode='determinate', value=0)
    percent_label_3 = Label(frame1, text='0.0 %')

    label_4 = Label(frame1, text="Reducing data and calculating statistics")
    progress_bar_4 = ttk.Progressbar(frame1, orient=HORIZONTAL,
                                     length=300, maximum=100, mode='determinate', value=0)
    percent_label_4 = Label(frame1, text='0.0 %')

    setup_window(frame1, [
        [[label_1, 0, 2]],
        [[progress_bar_1, 0, 1, 1, (20, 0)], [percent_label_1, 1]],
        [[label_2, 0, 2]],
        [[progress_bar_2, 0, 1, 1, (20, 0)], [percent_label_2, 1]],
        [[label_3, 0, 2]],
        [[progress_bar_3, 0, 1, 1, (20, 0)], [percent_label_3, 1]],
        [[label_4, 0, 2]],
        [[progress_bar_4, 0, 1, 1, (20, 0)], [percent_label_4, 1]],
        []
    ], main_font='Courier')

    canvas.draw()
    show_progress.after(200, run)
    show_progress.loop()
