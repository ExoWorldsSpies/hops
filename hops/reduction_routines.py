from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .hops_basics import *


def initialise_window(window, window_name=None, exit_command=None):

    if not window_name:
        window_name = read_log('windows', 'software_window')

    if not exit_command:
        def exit_command():
            os._exit(-1)

    window.wm_title(window_name)
    window.protocol('WM_DELETE_WINDOW', exit_command)

    window.withdraw()


def setup_window(window, objects):

    main_font = tuple(read_log('windows', 'main_font'))
    title_font = tuple(read_log('windows', 'title_font'))
    button_font = tuple(read_log('windows', 'button_font'))
    entries_bd = read_log('windows', 'entries_bd')

    for row in range(len(objects)):
        if len(objects[row]) == 0:
            label_empty = Label(window, text='')
            label_empty.grid(row=row, column=100)
        else:
            for obj in objects[row]:

                if obj[0].winfo_class() == 'Button':
                    obj[0].configure(font=button_font)
                elif obj[0].winfo_class() == 'Entry':
                    obj[0].configure(bd=entries_bd, font=main_font)
                elif obj[0].winfo_class() in ['Label', 'Radiobutton']:
                    if len(obj) == 5:
                        if obj[4] == 'title':
                            obj[0].configure(font=title_font)
                        else:
                            obj[0].configure(font=main_font)
                    else:
                        obj[0].configure(font=main_font)

                if len(obj) >= 4:
                    obj[0].grid(row=row, column=obj[1], columnspan=obj[2], rowspan=obj[3])
                elif len(obj) == 3:
                    obj[0].grid(row=row, column=obj[1], columnspan=obj[2])
                else:
                    obj[0].grid(row=row, column=obj[1])


def finalise_window(window, center=True, topmost=False):

    window.update_idletasks()

    if center:
        x = (window.winfo_screenwidth() - window.winfo_reqwidth()) / 2
        y = (window.winfo_screenheight() - window.winfo_reqheight()) / 2
        window.geometry('+%d+%d' % (x, y))

    else:
        window.geometry('+%d+%d' % (0, 0))

    window.update_idletasks()

    window.lift()
    window.wm_attributes("-topmost", 1)
    # if not topmost:
    window.after_idle(window.attributes, '-topmost', 0)

    window.deiconify()


def reduction():

    if read_local_log('pipeline', 'reduction_complete'):
        if not askyesno('Overwrite files', 'Reduction has been completed, do you want to run again?'):
            return 0

    write_local_log('pipeline', False, 'reduction_complete')
    write_local_log('pipeline', False, 'alignment_complete')

    # get variables

    observation_files = read_local_log('pipeline', 'observation_files')
    reduction_directory = read_local_log('pipeline', 'reduction_directory')
    reduction_prefix = read_local_log('pipeline', 'reduction_prefix')
    exposure_time_key = read_local_log('pipeline_keywords', 'exposure_time_key')
    mean_key = read_local_log('pipeline_keywords', 'mean_key')
    std_key = read_local_log('pipeline_keywords', 'std_key')
    observation_date_key = read_local_log('pipeline_keywords', 'observation_date_key')
    observation_time_key = read_local_log('pipeline_keywords', 'observation_time_key')
    frame_low_std = read_local_log('windows', 'frame_low_std')
    frame_upper_std = read_local_log('windows', 'frame_upper_std')
    bias_files = read_local_log('reduction', 'bias_files')
    dark_files = read_local_log('reduction', 'dark_files')
    flat_files = read_local_log('reduction', 'flat_files')
    bin_fits = int(read_local_log('reduction', 'bin_fits'))
    bin_to = int(read_local_log('reduction', 'bin_to'))
    master_bias_method = read_local_log('reduction', 'master_bias_method')
    master_dark_method = read_local_log('reduction', 'master_dark_method')
    master_flat_method = read_local_log('reduction', 'master_flat_method')
    target_ra_dec = read_local_log('photometry', 'target_ra_dec')

    # check if reduction directory exists

    if os.path.isdir(reduction_directory):
        shutil.rmtree(reduction_directory)

    os.mkdir(reduction_directory)

    # create master bias

    bias_frames = []
    if len(bias_files) > 0:
        for bias_file in find_fits_files(bias_files):
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

            try:
                bias_frames.append(float(fits[0].header['BZERO']) +
                                   float(fits[0].header['BSCALE']) * np.ones_like(fits[0].data) * fits[0].data)
            except KeyError:
                bias_frames.append(1.0 * np.ones_like(fits[0].data) * fits[0].data)

    if len(bias_frames) > 0:
        if master_bias_method == 'median':
            master_bias = np.median(bias_frames, 0)
        elif master_bias_method == 'mean':
            master_bias = np.mean(bias_frames, 0)
        else:
            master_bias = np.median(bias_frames, 0)
        master_bias_check = False
    else:
        master_bias = 0.0
        master_bias_check = True

    print(np.median(master_bias), master_bias_check)

    # create master dark

    dark_frames = []
    if len(str(dark_files)) > 0:
        for dark_file in find_fits_files(dark_files):
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

            try:
                if master_bias_check:
                    dark_frame = float(fits[0].header['BSCALE']) * np.ones_like(fits[0].data) * fits[0].data
                else:
                    dark_frame = float(fits[0].header['BZERO']) + \
                                 float(fits[0].header['BSCALE']) * np.ones_like(fits[0].data) * fits[0].data
            except KeyError:
                dark_frame = 1.0 * np.ones_like(fits[0].data) * fits[0].data
            dark_frames.append((dark_frame - master_bias) / fits[0].header[exposure_time_key])

    if len(dark_frames) > 0:
        if master_dark_method == 'median':
            master_dark = np.median(dark_frames, 0)
        elif master_dark_method == 'mean':
            master_dark = np.mean(dark_frames, 0)
        else:
            master_dark = np.median(dark_frames, 0)
    else:
        master_dark = 0.0

    print(np.median(master_dark))

    # create master flat

    flat_frames = []
    if len(str(flat_files)) > 0:
        for flat_file in find_fits_files(flat_files):
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

            try:
                if master_bias_check:
                    flat_frame = float(fits[0].header['BSCALE']) * np.ones_like(fits[0].data) * fits[0].data
                else:
                    flat_frame = float(fits[0].header['BZERO']) + \
                                 float(fits[0].header['BSCALE']) * np.ones_like(fits[0].data) * fits[0].data
            except KeyError:
                flat_frame = 1.0 * np.ones_like(fits[0].data) * fits[0].data
            flat_frames.append(flat_frame - master_bias - fits[0].header[exposure_time_key] * master_dark)

    if len(flat_frames) > 0:
        if master_flat_method == 'median':
            flat_frames = [ff / np.median(ff) for ff in flat_frames]
            master_flat = np.median(flat_frames, 0)
        elif master_flat_method == 'mean':
            master_flat = np.mean(flat_frames, 0)
        else:
            flat_frames = [ff / np.median(ff) for ff in flat_frames]
            master_flat = np.median(flat_frames, 0)
        master_flat = master_flat / np.median(master_flat)
    else:
        master_flat = 1.0

    print(np.median(master_flat))

    # setup counter window

    root = Tk()

    exit_var = BooleanVar(value=False)

    def break_and_exit():
        exit_var.set(True)

    initialise_window(root, exit_command=break_and_exit)

    f = Figure()
    f.patch.set_facecolor('white')
    ax = f.add_subplot(111)
    ax.axis('off')
    canvas = FigureCanvasTkAgg(f, root)
    canvas.get_tk_widget().pack()

    frame1 = Frame(root)
    frame1.pack()

    label1 = Label(frame1, text='REDUCTION')
    label2 = Label(frame1, text='FILE:')
    label3 = Label(frame1, text=' ')
    label4 = Label(frame1, text='COMPLETE:')
    label5 = Label(frame1, text=' ')
    label6 = Label(frame1, text='TIME LEFT:')
    label7 = Label(frame1, text=' ')

    setup_window(frame1, [
        [[label1, 1, 2]],
        [[label2, 1], [label3, 2]],
        [[label4, 1], [label5, 2]],
        [[label6, 1], [label7, 2]],
    ])

    # correct each observation_files file

    observation_files = find_fits_files(observation_files)
    percent = 0
    lt0 = time.time()

    testx = []
    testy = []
    testz = []

    for counter, science_file in enumerate(observation_files):

        # correct it with master bias_files, master dark_files and master flat_files

        fits = pf.open(science_file, memmap=False)

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

        try:
            if master_bias_check:
                data_frame = float(fits[0].header['BSCALE']) * np.ones_like(fits[0].data) * fits[0].data
            else:
                data_frame = float(fits[0].header['BZERO']) + float(fits[0].header['BSCALE']) * np.ones_like(fits[0].data) * fits[0].data
        except KeyError:
            data_frame = 1.0 * np.ones_like(fits[0].data) * fits[0].data

        data_frame = (data_frame - master_bias - fits[0].header[exposure_time_key] * master_dark) / master_flat

        if bin_fits > 1:
            data_frame = bin_frame(data_frame, bin_fits)

        fits[0].header.set('BZERO', 0.0)
        fits[0].header.set('BSCALE', 1.0)

        norm, floor, mean, std = fit_distribution1d_gaussian(data_frame, binning=bin_to/100000.0)

        if np.isnan(norm):
            mean = np.mean(data_frame)
            std = np.std(data_frame)

        if observation_date_key == observation_time_key:
            local_time = ' '.join(fits[0].header[observation_date_key].split('T'))
            if counter == 1:
                write_local_log('fitting', fits[0].header[observation_date_key].split('T')[0], 'date')
        else:
            local_time = ' '.join([fits[0].header[observation_date_key].split('T')[0],
                                   fits[0].header[observation_time_key]])
            if counter == 1:
                write_local_log('fitting', fits[0].header[observation_date_key].split('T')[0], 'date')

        ra_target, dec_target = ra_dec_string_to_deg(target_ra_dec)

        heliocentric_julian_date = plc.ut_to_hjd(ra_target, dec_target, local_time)

        testx.append(heliocentric_julian_date)
        testy.append(mean / fits[0].header[exposure_time_key])
        testz.append(std)

        fits[0].header.set(mean_key, mean)
        fits[0].header.set(std_key, std)

        # write the new fits file

        if observation_date_key == observation_time_key:
                local_time = fits[0].header[observation_date_key]
                local_time = '{0}_'.format(local_time.replace('-', '_').replace('T', '_').replace(':', '_'))
        else:
                local_time = '{0}_{1}_'.format(fits[0].header[observation_date_key].split('T')[0].replace('-', '_'),
                                               fits[0].header[observation_time_key].replace(':', '_'))

        hdu = pf.ImageHDU(header=fits[0].header, data=np.array(np.round(data_frame, 0), dtype=np.int32))
        hdu.writeto('{0}{1}{2}{3}{4}'.format(reduction_directory,
                                             os.sep, reduction_prefix, local_time, science_file.split(os.sep)[-1]))

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
            label3.configure(text='     {0}     '.format(science_file.split(os.sep)[-1]))
            label5.configure(text='     {0}%    '.format(new_percent))
            label7.configure(text='     %dh %02dm %02ds     ' % (int(hours), int(minutes), int(seconds)))
            percent = new_percent

        if counter == 0:
            finalise_window(root, topmost=True)

        root.update()

        if exit_var.get():
            break

        if counter + 1 == len(observation_files):
            write_local_log('pipeline', True, 'reduction_complete')

    testx = np.array(np.array(testx) - testx[0]) * 24.0 * 60.0

    root.destroy()

    # TODO discard all in frame

    if not exit_var.get():

        reduction_trash_directory = read_local_log('pipeline', 'reduction_trash_directory')
        trash = read_local_log('pipeline', 'trash')
        if not trash:
            list_to_remove = []
        else:
            list_to_remove = np.int_(trash)

        root = Tk()
        f = Figure()
        ax = f.add_subplot(2, 1, 1)
        ax2 = f.add_subplot(2, 2, 3)
        ax3 = f.add_subplot(2, 2, 4)

        exit_var_2 = BooleanVar(value=False)

        def break_and_exit():
            exit_var_2.set(True)

        initialise_window(root, exit_command=break_and_exit)

        f.patch.set_facecolor('white')
        canvas = FigureCanvasTkAgg(f, root)
        canvas.get_tk_widget().pack()
        NavigationToolbar2TkAgg(canvas, root)

        ax.plot(testx, testy, 'ko', ms=3)
        for ii in list_to_remove:
            ax.plot(testx[ii], testy[ii], 'ro', ms=3)

        time_dt = np.median(np.array(testx[1:]) - np.array(testx[:-1]))
        arrow = mpatches.Arrow(testx[0], 0, 0, testy[0], width=time_dt, fc='r')
        ax.add_patch(arrow)
        ax.set_xlabel('time (m)')
        ax.set_ylabel('sky (counts / pixel / second)')
        ax.tick_params(top='on', bottom='off', labeltop='on', labelbottom='off')
        ax.xaxis.set_label_position('top')

        science = find_fits_files(os.path.join(reduction_directory, '*'))

        fits = pf.open(science[0], memmap=False)

        ax2.imshow(fits[1].data, origin='lower', cmap=cm.Greys_r,
                   vmin=fits[1].header[mean_key] + frame_low_std * fits[1].header[std_key],
                   vmax=fits[1].header[mean_key] + frame_upper_std * fits[1].header[std_key])
        ax2.axis('off')

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

        def run_alignment(event):
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
                    exit_var_2.set(True)

        f.canvas.callbacks.connect('button_press_event', update_window_show)
        f.canvas.callbacks.connect('button_press_event', run_alignment)

        finalise_window(root, topmost=True)

        while not exit_var_2.get():
            root.update()

        root.destroy()
        if not os.path.isdir(os.path.join(reduction_directory, reduction_trash_directory)):
            os.mkdir(os.path.join(reduction_directory, reduction_trash_directory))

        for iii in list_to_remove:
            shutil.move(science[iii], os.path.join(reduction_directory, reduction_trash_directory))

        write_local_log('pipeline', list(map(str, list_to_remove)), 'trash')
