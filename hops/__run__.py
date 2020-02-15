from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .hops_basics import *
from .reduction_routines import reduction as rdr_reduction
from .alignment_routines import alignment as alr_alignment
from .photometry_routines import photometry as phr_photometry
from .fitting_routines import fitting as ftr_fitting
from .observing_planner import run_observing_planner

class AddOnWindow:

    def __init__(self, name, sizex, sizey, position=5):

        self.root = Tk()
        self.root.wm_title(name)

        self.root.protocol('WM_DELETE_WINDOW', self.root.withdraw)
        if sizex and sizey:
            self.root.geometry('{0}x{1}'.format(int(self.root.winfo_screenwidth() / sizex),
                                                int(self.root.winfo_screenheight() / sizey)))

        self.root.withdraw()
        self.finalised = False
        self.position = position

    def mainloop(self):

        self.root.mainloop()

    def finalise(self):

        self.root.update_idletasks()

        if self.position == 1:
            x = 0
            y = 0

        elif self.position == 2:
            x = (self.root.winfo_screenwidth() - self.root.winfo_reqwidth()) / 2
            y = 0

        elif self.position == 3:
            x = self.root.winfo_screenwidth() - self.root.winfo_reqwidth()
            y = 0

        elif self.position == 4:
            x = 0
            y = (self.root.winfo_screenheight() - self.root.winfo_reqheight()) / 2

        elif self.position == 5:
            x = (self.root.winfo_screenwidth() - self.root.winfo_reqwidth()) / 2
            y = (self.root.winfo_screenheight() - self.root.winfo_reqheight()) / 2

        elif self.position == 6:
            x = self.root.winfo_screenwidth() - self.root.winfo_reqwidth()
            y = (self.root.winfo_screenheight() - self.root.winfo_reqheight()) / 2

        elif self.position == 7:
            x = 0
            y = self.root.winfo_screenheight() - self.root.winfo_reqheight()

        elif self.position == 8:
            x = (self.root.winfo_screenwidth() - self.root.winfo_reqwidth()) / 2
            y = self.root.winfo_screenheight() - self.root.winfo_reqheight()

        elif self.position == 9:
            x = self.root.winfo_screenwidth() - self.root.winfo_reqwidth()
            y = self.root.winfo_screenheight() - self.root.winfo_reqheight()

        else:
            x = 0
            y = 0

        self.root.geometry('+%d+%d' % (int(x), int(y)))

        self.root.update_idletasks()

        self.root.lift()
        self.root.wm_attributes("-topmost", 1)
        self.root.after_idle(self.root.attributes, '-topmost', 0)

    def show(self):

        if not self.finalised:
            self.finalise()
            self.finalised = True

        self.root.deiconify()

    def hide(self):

        self.root.withdraw()

    def close(self):

        self.root.destroy()


def initialise_window(window, window_name, run, other_exit_command=None):

    def exit_command():

        if other_exit_command:
            other_exit_command()

        if run:
            run.exit = True
            os._exit(-1)

    window.wm_title(window_name)
    window.protocol('WM_DELETE_WINDOW', exit_command)

    window.withdraw()


def setup_window(window, objects, title_font=None, main_font=None, button_font=None, entries_bd=3, buttons_bd=5):
    screenheigth = window.winfo_screenheight()

    if button_font is None:
        button_font = ['times', int(screenheigth/55), 'bold']

    if main_font is None:
        main_font = ['times', int(screenheigth/60)]

    if title_font is None:
        title_font = ['times', int(screenheigth/40), 'bold']

    for row in range(len(objects)):
        if len(objects[row]) == 0:
            label_empty = Label(window, text='')
            label_empty.grid(row=row, column=100)
        else:
            for obj in objects[row]:

                if obj[0].winfo_class() == 'Button':
                    obj[0].config(borderwidth=buttons_bd, font=button_font, padx=3, pady=3)
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


def finalise_window(window, position=5):

    window.update_idletasks()

    if position == 1:
        x = 0
        y = 0

    elif position == 2:
        x = (window.winfo_screenwidth() - window.winfo_reqwidth()) / 2
        y = 0

    elif position == 3:
        x = window.winfo_screenwidth() - window.winfo_reqwidth()
        y = 0

    elif position == 4:
        x = 0
        y = (window.winfo_screenheight() - window.winfo_reqheight()) / 2

    elif position == 5:
        x = (window.winfo_screenwidth() - window.winfo_reqwidth()) / 2
        y = (window.winfo_screenheight() - window.winfo_reqheight()) / 2

    elif position == 6:
        x = window.winfo_screenwidth() - window.winfo_reqwidth()
        y = (window.winfo_screenheight() - window.winfo_reqheight()) / 2

    elif position == 7:
        x = 0
        y = window.winfo_screenheight() - window.winfo_reqheight()

    elif position == 8:
        x = (window.winfo_screenwidth() - window.winfo_reqwidth()) / 2
        y = window.winfo_screenheight() - window.winfo_reqheight()

    elif position == 9:
        x = window.winfo_screenwidth() - window.winfo_reqwidth()
        y = window.winfo_screenheight() - window.winfo_reqheight()

    else:
        x = 0
        y = 0

    window.geometry('+%d+%d' % (int(x), int(y)))

    window.update_idletasks()

    window.lift()
    window.wm_attributes("-topmost", 1)
    window.after_idle(window.attributes, '-topmost', 0)

    window.deiconify()


def reduction_alignment_window(run):

    # #########
    # create and initialise the window
    # #########

    root = Tk()
    initialise_window(root, 'Reduction & Alignment', run)

    # get variables from log and set as tk variables those to be modified

    directory = StringVar(root, value=read_local_log_user('directory'))
    directory_short = StringVar(root, value=read_local_log_user('directory_short'))
    observation_files = StringVar(root, value=read_local_log_profile('observation_files'))
    bias_files = StringVar(root, value=read_local_log_profile('bias_files'))
    dark_files = StringVar(root, value=read_local_log_profile('dark_files'))
    flat_files = StringVar(root, value=read_local_log_profile('flat_files'))
    bin_fits = StringVar(root, value=read_local_log_profile('bin_fits'))
    exposure_time_key = StringVar(root, value=read_local_log_profile('exposure_time_key'))
    observation_date_key = StringVar(root, value=read_local_log_profile('observation_date_key'))
    observation_time_key = StringVar(root, value=read_local_log_profile('observation_time_key'))
    target_ra_dec = StringVar(root, value=read_log('photometry', 'target_ra_dec'))
    auto_target_ra_dec = StringVar(root, value=read_log('photometry', 'auto_target_ra_dec'))
    use_auto_target_ra_dec = BooleanVar(root, value=read_log('photometry', 'use_auto_target_ra_dec'))
    mid_exposure = BooleanVar(root, value=read_log('photometry', 'mid_exposure'))

    # set progress variables, useful for updating the window

    update_directory = BooleanVar(root, value=False)
    running = BooleanVar(root, value=False)

    # create widgets

    observing_planner_button = Button(root, text='OBSERVATION\nPLANNER')
    my_profile_button = Button(root, text='MY PROFILE')

    directory_label = Label(root, text='Directory')
    directory_entry = Button(root, textvariable=directory_short)

    observation_files_label = Label(root, text='Name identifier for observation files')
    observation_files_entry = Entry(root, textvariable=observation_files)
    observation_files_test = Label(root, text=' ')

    bias_files_label = Label(root, text='Name identifier for bias files')
    bias_files_entry = Entry(root, textvariable=bias_files)
    bias_files_test = Label(root, text=' ')

    dark_files_label = Label(root, text='Name identifier for dark files')
    dark_files_entry = Entry(root, textvariable=dark_files)
    dark_files_test = Label(root, text=' ')

    flat_files_label = Label(root, text='Name identifier for flat files')
    flat_files_entry = Entry(root, textvariable=flat_files)
    flat_files_test = Label(root, text=' ')

    bin_fits_label = Label(root, text='Bin fits files (reduced only)')
    bin_fits_entry = Entry(root, textvariable=bin_fits, validate='key')
    bin_fits_entry['validatecommand'] = (bin_fits_entry.register(test_int_positive_non_zero_input), '%P', '%d')

    show_files_button = Button(root, text='Show files')

    exposure_time_key_label = Label(root, text='Exposure time header keyword')
    exposure_time_key_entry = Entry(root, textvariable=exposure_time_key)
    exposure_time_key_test = Label(root, text=' ')

    observation_date_key_label = Label(root, text='Observation date header keyword')
    observation_date_key_entry = Entry(root, textvariable=observation_date_key)
    observation_date_key_test = Label(root, text=' ')

    observation_time_key_label = Label(root, text='Observation time header keyword')
    observation_time_key_entry = Entry(root, textvariable=observation_time_key)
    observation_time_key_test = Label(root, text=' ')

    auto_target_ra_dec_label = Label(root, text='Detected target RA DEC')
    auto_target_ra_dec_entry = Label(root, text=auto_target_ra_dec.get())
    use_auto_target_ra_dec_entry = Checkbutton(root, text='Use detected values', variable=use_auto_target_ra_dec)

    target_ra_dec_label = Label(root, text='Manual target RA DEC\n(hh:mm:ss +/-dd:mm:ss)')
    target_ra_dec_entry = Entry(root, textvariable=target_ra_dec)
    target_ra_dec_test = Label(root, text=' ')

    mid_exposure_entry = Checkbutton(root, text='If the time stamp in your fits files refers to the mid-exposure\n'
                                                'instead of the exposure start, please tick this box.',
                                     variable=mid_exposure)

    show_header_button = Button(root, text='Show header')

    run_reduction_alignment_button = Button(root, text='RUN REDUCTION & ALIGNMENT')

    # define additional windows

    show_content_window = AddOnWindow('Files list', 3, 3, 1)
    scrollbar = Scrollbar(show_content_window.root)
    scrollbar.pack(side=RIGHT, fill=Y)
    files_list = Listbox(show_content_window.root, yscrollcommand=scrollbar.set, font='Courier')
    files_list.pack(side=LEFT, fill=BOTH, expand=True)
    scrollbar.config(command=files_list.yview)

    show_header_window = AddOnWindow('Header keywords list', 3, 3, 7)
    scrollbar = Scrollbar(show_header_window.root)
    scrollbar.pack(side=RIGHT, fill=Y)
    header_list = Listbox(show_header_window.root, yscrollcommand=scrollbar.set, font='Courier')
    header_list.pack(side=LEFT, fill=BOTH, expand=True)
    scrollbar.config(command=header_list.yview)

    my_profile_window = AddOnWindow('My Profile', 2, 1.1, 1)

    core_headers = {ff.split(':')[0]: read_log_profile(ff.split(':')[0])
                    for ff in open(log_profile_file, 'r').readlines()}

    local_headers = yaml.load(open(local_log_profile_file, 'r'), Loader=yaml.SafeLoader)

    variables = {}
    labels = {}
    entries = {}
    for row, header in enumerate(core_headers):
        if header in local_headers:
            variables[header] = StringVar(my_profile_window.root, value=local_headers[header])
            labels[header] = Label(my_profile_window.root, text=header)
            entries[header] = Entry(my_profile_window.root, textvariable=variables[header])
        else:
            variables[header] = StringVar(my_profile_window.root, value=core_headers[header])
            labels[header] = Label(my_profile_window.root, text=header)
            entries[header] = Entry(my_profile_window.root, textvariable=variables[header])

    def update_headers():
        new_local_profile = {}
        for header2 in variables:
            new_local_profile[header2] = variables[header2].get()
        yaml.dump(new_local_profile, open(local_log_profile_file, 'w'), default_flow_style=False)
        update_window(None)

    update_headers_button = Button(my_profile_window.root, text='UPDATE')
    update_headers_button['command'] = update_headers

    stucture = [[], [[update_headers_button, 2]], []]
    for header in list(core_headers.keys())[:int(len(list(core_headers.keys()))/2)]:
        stucture.append([[labels[header], 1], [entries[header], 2]])

    stucture.append([])

    for jj, header in enumerate(list(core_headers.keys())[int(len(list(core_headers.keys()))/2):]):
        stucture[3+jj].append([labels[header], 3])
        stucture[3+jj].append([entries[header], 4])

    stucture.append([])

    setup_window(my_profile_window.root, stucture)

    # define the function that updates the window

    def update_window(*entry):

        if not entry:
            pass

        if running.get():

            directory_entry['state'] = DISABLED
            observation_files_entry['state'] = DISABLED
            bias_files_entry['state'] = DISABLED
            dark_files_entry['state'] = DISABLED
            flat_files_entry['state'] = DISABLED
            bin_fits_entry['state'] = DISABLED
            show_files_button['state'] = DISABLED
            exposure_time_key_entry['state'] = DISABLED
            observation_date_key_entry['state'] = DISABLED
            observation_time_key_entry['state'] = DISABLED
            target_ra_dec_entry['state'] = DISABLED
            use_auto_target_ra_dec_entry['state'] = DISABLED
            mid_exposure_entry['state'] = DISABLED
            show_header_button['state'] = DISABLED
            run_reduction_alignment_button['state'] = DISABLED
            observing_planner_button['state'] = DISABLED
            my_profile_button['state'] = DISABLED

        elif not os.path.isdir(directory.get()):

            directory_entry['state'] = NORMAL
            observation_files_entry['state'] = DISABLED
            bias_files_entry['state'] = DISABLED
            dark_files_entry['state'] = DISABLED
            flat_files_entry['state'] = DISABLED
            bin_fits_entry['state'] = DISABLED
            show_files_button['state'] = DISABLED
            exposure_time_key_entry['state'] = DISABLED
            observation_date_key_entry['state'] = DISABLED
            observation_time_key_entry['state'] = DISABLED
            target_ra_dec_entry['state'] = DISABLED
            use_auto_target_ra_dec_entry['state'] = DISABLED
            mid_exposure_entry['state'] = DISABLED
            show_header_button['state'] = DISABLED
            run_reduction_alignment_button['state'] = DISABLED
            observing_planner_button['state'] = NORMAL
            my_profile_button['state'] = NORMAL

        else:

            directory_entry['state'] = NORMAL

            if update_directory.get():

                os.chdir(directory.get())
                initiate_local_log_file()

                observation_files.set(read_local_log('pipeline', 'observation_files'))
                bias_files.set(read_local_log('reduction', 'bias_files'))
                dark_files.set(read_local_log('reduction', 'dark_files'))
                flat_files.set(read_local_log('reduction', 'flat_files'))
                bin_fits.set(read_local_log('reduction', 'bin_fits'))
                exposure_time_key.set(read_local_log('pipeline_keywords', 'exposure_time_key'))
                observation_date_key.set(read_local_log('pipeline_keywords', 'observation_date_key'))
                observation_time_key.set(read_local_log('pipeline_keywords', 'observation_time_key'))
                target_ra_dec.set(read_local_log('photometry', 'target_ra_dec'))
                auto_target_ra_dec.set(read_local_log('photometry', 'auto_target_ra_dec'))
                auto_target_ra_dec_entry.configure(text=auto_target_ra_dec.get())
                use_auto_target_ra_dec.set(read_local_log('photometry', 'use_auto_target_ra_dec'))
                mid_exposure.set(read_local_log('photometry', 'mid_exposure'))

            observation_files_entry['state'] = NORMAL
            bias_files_entry['state'] = NORMAL
            dark_files_entry['state'] = NORMAL
            flat_files_entry['state'] = NORMAL
            bin_fits_entry['state'] = NORMAL
            show_files_button['state'] = NORMAL
            observing_planner_button['state'] = NORMAL
            my_profile_button['state'] = NORMAL

            files_list.delete(0, END)
            files_list.insert(END, '  List of files in your directory:')
            files_list.insert(END, '  ')

            xx = find_fits_files('*')

            for ii in xx:
                files_list.insert(END, '  {0}'.format(str(ii).split(os.sep)[-1]))

            check_science = test_file_number(observation_files_entry.get())
            observation_files_test.configure(text=check_science[1])

            check_bias = test_file_number(bias_files_entry.get())
            bias_files_test.configure(text=check_bias[1])

            check_dark = test_file_number(dark_files_entry.get())
            dark_files_test.configure(text=check_dark[1])

            check_flat = test_file_number(flat_files_entry.get())
            flat_files_test.configure(text=check_flat[1])

            if not check_science[0]:

                exposure_time_key_entry['state'] = DISABLED
                observation_date_key_entry['state'] = DISABLED
                observation_time_key_entry['state'] = DISABLED
                target_ra_dec_entry['state'] = DISABLED
                use_auto_target_ra_dec_entry['state'] = DISABLED
                mid_exposure_entry['state'] = DISABLED
                show_header_button['state'] = DISABLED
                run_reduction_alignment_button['state'] = DISABLED

            else:

                target_ra_dec_entry['state'] = NORMAL
                use_auto_target_ra_dec_entry['state'] = NORMAL
                mid_exposure_entry['state'] = NORMAL
                exposure_time_key_entry['state'] = NORMAL
                observation_date_key_entry['state'] = NORMAL
                observation_time_key_entry['state'] = NORMAL
                show_header_button['state'] = NORMAL

                fits = pf.open(find_fits_files(observation_files_entry.get())[0], memmap=False)

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

                science_header = fits[0].header

                header_list.delete(0, END)
                header_list.insert(END, '  Keywords:      Values:')
                header_list.insert(END, '  ')

                for ii in science_header:
                    if ii != '':
                        header_list.insert(END, '  {0}{1}{2}'.format(str(ii[:10]), ' ' * (15 - len(str(ii[:10]))),
                                                                     str(science_header[ii])))

                check_ra = [None, None]
                for key in read_local_log_profile('target_ra_key').split(','):
                    check_ra = test_fits_keyword(observation_files_entry.get(), key)
                    if check_ra[0]:
                        break

                check_dec = [None, None]
                for key in read_local_log_profile('target_dec_key').split(','):
                    check_dec = test_fits_keyword(observation_files_entry.get(), key)
                    if check_dec[0]:
                        break

                if check_ra[0] and check_dec[0]:
                    try:
                        if isinstance(check_ra[2], str):
                            target = plc.Target(plc.Hours(check_ra[2].replace(',', '.')), plc.Degrees(check_dec[2].replace(',', '.')))
                        elif isinstance(check_ra[2], float):
                            target = plc.Target(plc.Degrees(check_ra[2]), plc.Degrees(check_dec[2]))
                        auto_target_ra_dec_entry.configure(text=target.coord)
                        use_auto_target_ra_dec_entry['state'] = NORMAL
                    except:
                        auto_target_ra_dec_entry.configure(text='None detected')
                        use_auto_target_ra_dec.set(0)
                        use_auto_target_ra_dec_entry['state'] = DISABLED
                else:
                    auto_target_ra_dec_entry.configure(text='None detected')
                    use_auto_target_ra_dec.set(0)
                    use_auto_target_ra_dec_entry['state'] = DISABLED

                if use_auto_target_ra_dec.get():
                    target_ra_dec.set(target.coord)
                    target_ra_dec_entry['state'] = DISABLED
                else:
                    target_ra_dec_entry['state'] = NORMAL

                check_ra_dec = test_coordinates(target_ra_dec_entry.get())
                target_ra_dec_test.configure(text=check_ra_dec[1])

                check_exposure_time = test_fits_keyword(observation_files_entry.get(),
                                                        exposure_time_key_entry.get())
                exposure_time_key_test.configure(text=check_exposure_time[1])

                check_observation_date = test_fits_keyword(observation_files_entry.get(),
                                                           observation_date_key_entry.get())
                observation_date_key_test.configure(text=check_observation_date[1])

                if check_observation_date[0]:
                    if len(check_observation_date[2].split('T')) == 2:
                        observation_time_key.set(observation_date_key_entry.get())
                        observation_time_key_entry['state'] = DISABLED

                check_observation_time = test_fits_keyword(observation_files_entry.get(),
                                                           observation_time_key_entry.get())
                observation_time_key_test.configure(text=check_observation_time[1])

                if (check_ra_dec[0] and check_exposure_time[0] and
                        check_observation_date[0] and check_observation_time[0]):

                    run_reduction_alignment_button['state'] = NORMAL

                else:

                    run_reduction_alignment_button['state'] = DISABLED

    update_directory = BooleanVar(root, value=True)
    update_window(None)
    update_directory = BooleanVar(root, value=False)

    # define actions for the different buttons, including calls to the function that updates the window

    def choose_directory():

        new_directory = tkFileDialog.askdirectory()

        if len(new_directory) > 0:
            directory.set(new_directory)
            directory_short.set('/'.join(new_directory.split('/')[-2:]))
            write_local_log_user('directory', directory.get())
            write_local_log_user('directory_short', directory_short.get())

            update_directory.set(True)
            update_window('a')
            update_directory.set(False)

    def run_reduction_alignment():

        running.set(True)
        update_window(None)

        write_local_log('pipeline', directory.get(), 'directory')
        write_local_log_user('directory', directory.get())
        write_local_log('pipeline', directory_short.get(), 'directory_short')
        write_local_log_user('directory_short', directory_short.get())
        write_local_log('pipeline', observation_files.get(), 'observation_files')
        write_local_log('reduction', bias_files.get(), 'bias_files')
        write_local_log('reduction', dark_files.get(), 'dark_files')
        write_local_log('reduction', flat_files.get(), 'flat_files')
        try:
            bins_test = int(bin_fits.get())
        except:
            bins_test = 1
        write_local_log('reduction', bins_test, 'bin_fits')
        write_local_log('photometry', target_ra_dec.get(), 'target_ra_dec')
        write_local_log('photometry', use_auto_target_ra_dec.get(), 'use_auto_target_ra_dec')
        write_local_log('photometry', mid_exposure.get(), 'mid_exposure')
        write_local_log('pipeline_keywords', exposure_time_key.get(), 'exposure_time_key')
        write_local_log('pipeline_keywords', observation_date_key.get(), 'observation_date_key')
        write_local_log('pipeline_keywords', observation_time_key.get(), 'observation_time_key')

        rdr_reduction()

        if read_local_log('pipeline', 'reduction_complete'):
            alr_alignment()

        if read_local_log('pipeline', 'reduction_complete') and read_local_log('pipeline', 'alignment_complete'):
            run.run_from_reduction = False
            run.run_from_photometry = True
            run.run_from_fitting = False
            root.destroy()
            show_content_window.close()
            show_header_window.close()
            my_profile_window.close()
        else:
            running.set(False)
            update_window(None)

    # connect actions to widgets

    observing_planner_button['command'] = run_observing_planner
    my_profile_button['command'] = my_profile_window.show
    directory_entry['command'] = choose_directory
    observation_files_entry.bind(sequence='<KeyRelease>', func=update_window)
    bias_files_entry.bind(sequence='<KeyRelease>', func=update_window)
    dark_files_entry.bind(sequence='<KeyRelease>', func=update_window)
    flat_files_entry.bind(sequence='<KeyRelease>', func=update_window)
    bin_fits_entry.bind(sequence='<KeyRelease>', func=update_window)
    show_files_button['command'] = show_content_window.show
    exposure_time_key_entry.bind(sequence='<KeyRelease>', func=update_window)
    observation_date_key_entry.bind(sequence='<KeyRelease>', func=update_window)
    observation_time_key_entry.bind(sequence='<KeyRelease>', func=update_window)
    use_auto_target_ra_dec_entry['command'] = update_window
    target_ra_dec_entry.bind(sequence='<KeyRelease>', func=update_window)
    show_header_button['command'] = show_header_window.show
    run_reduction_alignment_button['command'] = run_reduction_alignment

    Btn = Button(root, text="HOPS UPDATES &\nUSER MANUAL", command=openweb)

    # setup window

    photo = PhotoImage(file=holomon_logo)
    logo_label = Label(root, image=photo)
    window_label = Label(root, text='Reduction & Alignment')
    created_by_label = Label(root, text=read_log('windows', 'created_by').replace(',', '\n'))

    setup_window(root, [
        [[logo_label, 0, 1, 8]],
        [],
        [[window_label, 1, 3, 1, 'title']],
        [],
        [[directory_label, 1], [directory_entry, 2, 2]],
        [[observation_files_label, 1], [observation_files_entry, 2], [observation_files_test, 3]],
        [[bias_files_label, 1], [bias_files_entry, 2], [bias_files_test, 3]],
        [[dark_files_label, 1], [dark_files_entry, 2], [dark_files_test, 3]],
        [[created_by_label, 0, 1, 3], [flat_files_label, 1], [flat_files_entry, 2], [flat_files_test, 3]],
        [[bin_fits_label, 1], [bin_fits_entry, 2]],
        [[show_files_button, 2]],
        [[my_profile_button, 0]],
        [[auto_target_ra_dec_label, 1], [auto_target_ra_dec_entry, 2], [use_auto_target_ra_dec_entry, 3]],
        [[observing_planner_button, 0], [target_ra_dec_label, 1], [target_ra_dec_entry, 2], [target_ra_dec_test, 3]],
        [[exposure_time_key_label, 1], [exposure_time_key_entry, 2], [exposure_time_key_test, 3]],
        [[Btn, 0, 1, 2], [observation_date_key_label, 1], [observation_date_key_entry, 2], [observation_date_key_test, 3]],
        [[observation_time_key_label, 1], [observation_time_key_entry, 2], [observation_time_key_test, 3]],
        [[mid_exposure_entry, 1, 3]],
        [[show_header_button, 2]],
        [],
        [[run_reduction_alignment_button, 1, 3]],
        [],
    ])

    # finalise and show window

    finalise_window(root, position=5)
    try:
        location = os.path.abspath(os.path.dirname(__file__))

        current_version = '0.0.0'
        for i in open(os.path.join(location, '__init__.py')):
            if len(i.split('__version__')) > 1:
                current_version = i.split()[-1][1:-1]

        c1 = int(current_version.split('.')[0]) * 100 * 100 * 100
        c2 = int(current_version.split('.')[1]) * 100 * 100
        c3 = int(current_version.split('.')[2]) * 100

        version = '0.0.0'
        message = ''
        for i in urlopen('https://raw.githubusercontent.com/ExoWorldsSpies/hops/master/hops/__init__.py').readlines():
            if len(str(i).split('__version__')) > 1:
                version = str(i).split()[-1][1:-4]
            if len(str(i).split('__message__')) > 1:
                message = str(i).split('__message__ = ')[-1][1:-4]

        v1 = int(version.split('.')[0]) * 100 * 100 * 100
        v2 = int(version.split('.')[1]) * 100 * 100
        v3 = int(version.split('.')[2]) * 100

        if v1 + v2 + v3 > c1 + c2 + c3:
            showinfo('Update available', 'There is a newer version ({0}) of the code available!\n\n{1}\n\nDownload and install it from:'
                                         '\nhttps://www.exoworldsspies.com/en/software'.format(version, message))
    except:
        pass
    show_content_window.mainloop()
    show_header_window.mainloop()
    my_profile_window.mainloop()
    root.mainloop()


def photometry_window(run):

    # #########
    # create and initialise window
    # #########

    root = Tk()
    initialise_window(root, 'Photometry', run)

    reduction_directory = read_local_log('pipeline', 'reduction_directory')
    light_curve_aperture_file = read_local_log('pipeline', 'light_curve_aperture_file')
    photometry_directory = read_local_log('pipeline', 'photometry_directory')
    fov_figure = read_local_log('pipeline', 'fov_figure')
    mean_key = read_local_log('pipeline_keywords', 'mean_key')
    std_key = read_local_log('pipeline_keywords', 'std_key')
    align_x0_key = read_local_log('pipeline_keywords', 'align_x0_key')
    align_y0_key = read_local_log('pipeline_keywords', 'align_y0_key')
    frame_low_std = read_local_log('windows', 'frame_low_std')
    frame_upper_std = read_local_log('windows', 'frame_upper_std')
    bin_fits = int(read_local_log('reduction', 'bin_fits'))
    burn_limit = int(read_local_log('alignment', 'burn_limit')) * bin_fits * bin_fits
    star_std = read_local_log('alignment', 'star_std')
    search_window_std = read_local_log('alignment', 'search_window_std')
    max_comparisons = read_local_log('photometry', 'max_comparisons')
    max_targets = max_comparisons + 1
    target_ra_dec = read_local_log('photometry', 'target_ra_dec')

    targets_x_position = [DoubleVar(root, value=read_local_log('photometry', 'target_x_position'))]
    for comparison in range(max_comparisons):
        targets_x_position.append(
            DoubleVar(root, value=read_local_log('photometry', 'comparison_{0}_x_position'.format(comparison + 1))))
        if not targets_x_position[-1].get():
            targets_x_position[-1].set(0)

    targets_y_position = [DoubleVar(root, value=read_local_log('photometry', 'target_y_position'))]
    for comparison in range(max_comparisons):
        targets_y_position.append(
            DoubleVar(root, value=read_local_log('photometry', 'comparison_{0}_y_position'.format(comparison + 1))))
        if not targets_y_position[-1].get():
            targets_y_position[-1].set(0)

    targets_peak_counts = [DoubleVar(root, value=0)]
    for comparison in range(max_comparisons):
        targets_peak_counts.append(DoubleVar(root, value=0))

    targets_aperture = [IntVar(root, value=read_local_log('photometry', 'target_aperture'))]
    for comparison in range(max_comparisons):
        targets_aperture.append(
            IntVar(root, value=read_local_log('photometry', 'comparison_{0}_aperture'.format(comparison + 1))))
        if not targets_aperture[-1].get():
            targets_aperture[-1].set(0)

    # set progress variables, useful for updating the tk windows

    targets_indication = IntVar(root, value=0)
    click_test_x = DoubleVar(root, value=-100)
    click_test_y = DoubleVar(root, value=-100)
    click_off_axis = BooleanVar(root, value=False)
    running = BooleanVar(root, value=False)

    # create widgets

    help_label =Label(root, text="Remember, the best comparison stars need to be:\n"
                                 "a) close to your target, b) of similar magnitude to the target,\n"
                                 "c) of similar colour to the target, d) photometrically stable, i.e. not variables!")

    position_x_label = Label(root, text='X')

    position_y_label = Label(root, text='Y')

    peak_counts_label = Label(root, text='Peak')

    box_semi_length_label = Label(root, text='Box semi-length')

    targets_indication_entry = [Radiobutton(root, text='      Target           ', variable=targets_indication, value=0)]
    for comparison in range(max_comparisons):
        targets_indication_entry.append(Radiobutton(root, text='Comparison {0}     '.format(comparison + 1),
                                                    variable=targets_indication, value=comparison + 1))

    targets_x_position_label = [Label(root, textvar=targets_x_position[0])]
    for comparison in range(max_comparisons):
        targets_x_position_label.append(Label(root, textvar=targets_x_position[comparison + 1]))

    targets_y_position_label = [Label(root, textvar=targets_y_position[0])]
    for comparison in range(max_comparisons):
        targets_y_position_label.append(Label(root, textvar=targets_y_position[comparison + 1]))

    targets_peak_counts_label = [Label(root, textvar=targets_peak_counts[0])]
    for comparison in range(max_comparisons):
        targets_peak_counts_label.append(Label(root, textvar=targets_peak_counts[comparison + 1]))

    targets_aperture_entry = [Entry(root, textvar=targets_aperture[0], validate='key')]
    for comparison in range(max_comparisons):
        targets_aperture_entry.append(Entry(root, textvar=targets_aperture[comparison + 1], validate='key'))

    for target in range(max_targets):
        targets_aperture_entry[target]['validatecommand'] = \
            (targets_aperture_entry[target].register(test_int_positive_non_zero_input), '%P', '%d')

    show_fov_button = Button(root, text='Show FOV')

    flip_fov_button = Button(root, text='Flip FOV')

    mirror_fov_button = Button(root, text='Mirror FOV')

    photometry_button = Button(root, text='RUN PHOTOMETRY')

    proceed_to_fitting_button = Button(root, text='PROCEED TO FITTING')

    return_to_reduction_button = Button(root, text='RETURN TO REDUCTION')

    # define additional windows

    show_fov_window = AddOnWindow('FOV', None, None, 1)
    test = find_fits_files(os.path.join(reduction_directory, '*'))

    fits = pf.open(test[0], memmap=False)
    f = Figure()
    f.patch.set_facecolor('white')
    ax = f.add_subplot(111)
    canvas = FigureCanvasTkAgg(f, show_fov_window.root)
    canvas.get_tk_widget().pack()
    NavigationToolbar2TkAgg(canvas, show_fov_window.root)

    ax.imshow(fits[1].data, origin='lower', cmap=cm.Greys_r,
              vmin=fits[1].header[mean_key] + frame_low_std * fits[1].header[std_key],
              vmax=fits[1].header[mean_key] + frame_upper_std * fits[1].header[std_key])

    targets_box = [mpatches.Rectangle((targets_x_position[0].get() - targets_aperture[0].get(),
                                       targets_y_position[0].get() - targets_aperture[0].get()),
                                      2 * targets_aperture[0].get() + 1, 2 * targets_aperture[0].get() + 1,
                                      ec='r', fill=False)]
    for comparison in range(max_comparisons):
        targets_box.append(mpatches.Rectangle((targets_x_position[comparison + 1].get() -
                                               targets_aperture[comparison + 1].get(),
                                               targets_y_position[comparison + 1].get() -
                                               targets_aperture[comparison + 1].get()),
                                              2 * targets_aperture[comparison + 1].get() + 1,
                                              2 * targets_aperture[comparison + 1].get() + 1,
                                              ec='#07fefc', fill=False))

    for target in range(max_targets):
        ax.add_patch(targets_box[target])

    targets_text = [ax.text(targets_x_position[0].get(),
                            targets_y_position[0].get() - targets_aperture[0].get() - 1, 'T',
                            color='r', fontsize=20, va='top')]

    for comparison in range(max_comparisons):
        targets_text.append(ax.text(targets_x_position[comparison + 1].get()
                                    + targets_aperture[comparison + 1].get() + 1,
                                    targets_y_position[comparison + 1].get()
                                    - targets_aperture[comparison + 1].get() - 1,
                                    'C{0}'.format(comparison + 1), color='#07fefc', fontsize=20, va='top'))

    # define the function that updates the window

    def update_window(event):

        if running.get():

            for i_target in range(max_targets):
                targets_indication_entry[i_target]['state'] = DISABLED
                targets_aperture_entry[i_target]['state'] = DISABLED

            show_fov_button['state'] = DISABLED
            flip_fov_button['state'] = DISABLED
            mirror_fov_button['state'] = DISABLED
            photometry_button['state'] = DISABLED
            proceed_to_fitting_button['state'] = DISABLED
            return_to_reduction_button['state'] = DISABLED

        else:

            show_fov_button['state'] = NORMAL
            flip_fov_button['state'] = NORMAL
            mirror_fov_button['state'] = NORMAL
            return_to_reduction_button['state'] = NORMAL

            for i_target in range(max_targets):
                targets_indication_entry[i_target]['state'] = NORMAL

            try:
                if event.inaxes is not None:

                    click_off_axis.set(False)

                    if (event.xdata, event.ydata) == (click_test_x.get(), click_test_y.get()):

                        star = plc.find_single_star(
                            fits[1].data, event.xdata, event.ydata,
                            mean=fits[1].header[mean_key], std=fits[1].header[std_key],
                            burn_limit=burn_limit * 7.0 / 8.0, star_std=star_std)

                        if star:

                            targets_x_position[targets_indication.get()].set(round(star[0], 1))
                            targets_y_position[targets_indication.get()].set(round(star[1], 1))
                            targets_aperture[targets_indication.get()].set(abs(int(3 * star[4])))

                        click_test_x.set(-100)
                        click_test_y.set(-100)

                    else:
                        click_test_x.set(event.xdata)
                        click_test_y.set(event.ydata)

                else:

                    click_test_x.set(-100)
                    click_test_y.set(-100)

                    if click_off_axis.get():

                        targets_x_position[targets_indication.get()].set(0)
                        targets_y_position[targets_indication.get()].set(0)
                        targets_aperture[targets_indication.get()].set(0)
                        targets_peak_counts[targets_indication.get()].set(0)

                        click_off_axis.set(False)

                    else:
                        click_off_axis.set(True)

            except AttributeError:
                pass

            try:

                for i_target in range(max_targets):

                    if 0 in [targets_x_position[i_target].get(), targets_y_position[i_target].get()]:

                        targets_box[i_target].set_xy((-10000, -10000))

                        targets_text[i_target].set_x(-10000)
                        targets_text[i_target].set_y(-10000)

                        targets_aperture_entry[i_target]['state'] = DISABLED

                    else:

                        targets_box[i_target].set_xy((targets_x_position[i_target].get() -
                                                      targets_aperture[i_target].get() - 0.5,
                                                      targets_y_position[i_target].get() -
                                                      targets_aperture[i_target].get() - 0.5))

                        targets_box[i_target].set_width(2 * targets_aperture[i_target].get() + 1)
                        targets_box[i_target].set_height(2 * targets_aperture[i_target].get() + 1)

                        targets_text[i_target].set_x(targets_x_position[i_target].get() +
                                                     targets_aperture[i_target].get() + 1)
                        targets_text[i_target].set_y(targets_y_position[i_target].get() -
                                                     targets_aperture[i_target].get() - 1)

                        y1 = int(targets_y_position[i_target].get() - targets_aperture[i_target].get())
                        y2 = y1 + int(2 * targets_aperture[i_target].get()) + 2
                        x1 = int(targets_x_position[i_target].get() - targets_aperture[i_target].get())
                        x2 = x1 + int(2 * targets_aperture[i_target].get()) + 2
                        targets_peak_counts[i_target].set(int(np.max(fits[1].data[y1:y2, x1:x2])))

                        targets_aperture_entry[i_target]['state'] = NORMAL

            except ValueError:
                photometry_button['state'] = DISABLED

            if 0 in [targets_x_position[0].get(), targets_y_position[0].get(),
                     targets_x_position[1].get(), targets_y_position[1].get()]:
                photometry_button['state'] = DISABLED

            else:
                photometry_button['state'] = NORMAL

            if (read_local_log('pipeline', 'photometry_complete')
               and len(glob.glob(os.path.join('{0}*'.format(photometry_directory), light_curve_aperture_file))) > 0):
                proceed_to_fitting_button['state'] = NORMAL

            else:
                proceed_to_fitting_button['state'] = DISABLED

        canvas.draw()

    update_window(None)

    # define actions for the different buttons, including calls to the function that updates the window

    def photometry():

        running.set(True)
        update_window(None)

        write_local_log('photometry', targets_x_position[0].get(), 'target_x_position')
        write_local_log('photometry', targets_y_position[0].get(), 'target_y_position')
        write_local_log('photometry', targets_aperture[0].get(), 'target_aperture')
        target_polar = plc.cartesian_to_polar(targets_x_position[0].get(), targets_y_position[0].get(),
                                              fits[1].header[align_x0_key], fits[1].header[align_y0_key])
        write_local_log('photometry', float(target_polar[0]), 'target_r_position')
        write_local_log('photometry', float(target_polar[1]), 'target_u_position')

        for i_comparison in range(max_comparisons):
            write_local_log('photometry', targets_x_position[i_comparison + 1].get(),
                            'comparison_{0}_x_position'.format(i_comparison + 1))
            write_local_log('photometry', targets_y_position[i_comparison + 1].get(),
                            'comparison_{0}_y_position'.format(i_comparison + 1))
            write_local_log('photometry', targets_aperture[i_comparison + 1].get(),
                            'comparison_{0}_aperture'.format(i_comparison + 1))

            if 0 not in [targets_x_position[i_comparison + 1].get(), targets_y_position[i_comparison + 1].get()]:

                target_polar = plc.cartesian_to_polar(targets_x_position[i_comparison + 1].get(),
                                                      targets_y_position[i_comparison + 1].get(),
                                                      fits[1].header[align_x0_key], fits[1].header[align_y0_key])

            else:

                target_polar = [0, 0]

            write_local_log('photometry', float(target_polar[0]), 'comparison_{0}_r_position'.format(i_comparison + 1))
            write_local_log('photometry', float(target_polar[1]), 'comparison_{0}_u_position'.format(i_comparison + 1))

        f.savefig(fov_figure, dpi=200)
        phr_photometry()

        running.set(False)
        update_window(None)

    def flip_fov():
        ax.set_ylim(ax.get_ylim()[1], ax.get_ylim()[0])
        canvas.draw()

    def mirror_fov():
        ax.set_xlim(ax.get_xlim()[1], ax.get_xlim()[0])
        canvas.draw()

    def proceed_to_fitting():
        run.run_from_reduction = False
        run.run_from_photometry = False
        run.run_from_fitting = True
        show_fov_window.close()
        root.destroy()

    def return_to_reduction():
        run.run_from_reduction = True
        run.run_from_photometry = False
        run.run_from_fitting = False
        show_fov_window.close()
        root.destroy()

    # connect actions to widgets

    f.canvas.callbacks.connect('button_press_event', update_window)
    show_fov_button['command'] = show_fov_window.show
    flip_fov_button['command'] = flip_fov
    mirror_fov_button['command'] = mirror_fov
    photometry_button['command'] = photometry
    proceed_to_fitting_button['command'] = proceed_to_fitting
    return_to_reduction_button['command'] = return_to_reduction

    for target in range(max_targets):
        targets_aperture_entry[target].bind(sequence='<KeyRelease>', func=update_window)

    # setup window

    Btn = Button(root, text="HOPS UPDATES &\nUSER MANUAL", command=openweb)

    Btn2 = Button(root, text="CHECK SIMBAD", command=openweb_simbad(target_ra_dec))

    photo = PhotoImage(file=holomon_logo)
    logo_label = Label(root, image=photo)
    window_label = Label(root, text='Photometry')
    created_by_label = Label(root, text=read_log('windows', 'created_by').replace(',', '\n'))

    setup_list = [
        [[logo_label, 0, 1, 8]],
        [],
        [[window_label, 1, 6, 1, 'title']],
        [[help_label, 1, 6]],
        [[Btn2, 1, 6]],
        [],
        [[position_x_label, 2], [position_y_label, 3], [peak_counts_label, 4], [box_semi_length_label, 5, 2]],
    ]

    for target in range(max_targets):

        if target == 0:
            setup_list.append([[created_by_label, 0, 1, 3],
                               [targets_indication_entry[target], 1], [targets_x_position_label[target], 2],
                               [targets_y_position_label[target], 3], [targets_peak_counts_label[target], 4],
                               [targets_aperture_entry[target], 5, 2]])
        elif target == 2:
            setup_list.append([[Btn, 0, 1, 3],
                               [targets_indication_entry[target], 1], [targets_x_position_label[target], 2],
                               [targets_y_position_label[target], 3], [targets_peak_counts_label[target], 4],
                               [targets_aperture_entry[target], 5, 2]])
        else:
            setup_list.append([[targets_indication_entry[target], 1], [targets_x_position_label[target], 2],
                               [targets_y_position_label[target], 3], [targets_peak_counts_label[target], 4],
                               [targets_aperture_entry[target], 5, 2]])

    setup_list.append([[show_fov_button, 5, 2]])
    setup_list.append([[flip_fov_button, 5], [mirror_fov_button, 6]])
    setup_list.append([])
    setup_list.append([[photometry_button, 1, 6]])
    setup_list.append([[proceed_to_fitting_button, 1, 3], [return_to_reduction_button, 4, 3]])
    setup_list.append([])

    setup_window(root, setup_list)

    # finalise and show window

    finalise_window(root, position=5)
    show_fov_window.mainloop()
    root.mainloop()


def fitting_window(run):

    # #########
    # create and initialise window
    # #########

    root = Tk()
    initialise_window(root, 'Fitting', run)

    # get variables from log and set as tk variables those to be modified

    catalogue = plc.oec_catalogue()

    observation_files = StringVar(root, value=read_local_log('pipeline', 'observation_files'))
    light_curve_file = StringVar(value=read_local_log('fitting', 'light_curve_file'))
    light_curve_aperture_file = read_local_log('pipeline', 'light_curve_aperture_file')
    light_curve_gauss_file = read_local_log('pipeline', 'light_curve_gauss_file')
    photometry_directory = read_local_log('pipeline', 'photometry_directory')
    light_curve_file.set(glob.glob(os.path.join('{0}*'.format(photometry_directory), light_curve_aperture_file))[-1])
    planet_search = StringVar(value=read_local_log('fitting', 'planet_search'))
    planet = StringVar(value=read_local_log('fitting', 'planet'))
    binning = IntVar(value=read_local_log('fitting', 'binning'))
    scatter = DoubleVar(value=read_local_log('fitting', 'scatter'))
    iterations = IntVar(value=read_local_log('fitting', 'iterations'))
    burn = IntVar(value=read_local_log('fitting', 'burn'))
    metallicity = DoubleVar(value=read_local_log('fitting', 'metallicity'))
    temperature = DoubleVar(value=read_local_log('fitting', 'temperature'))
    logg = DoubleVar(value=read_local_log('fitting', 'logg'))
    period = DoubleVar(value=read_local_log('fitting', 'period'))
    mid_time = DoubleVar(value=read_local_log('fitting', 'mid_time'))
    rp_over_rs = DoubleVar(value=read_local_log('fitting', 'rp_over_rs'))
    sma_over_rs = DoubleVar(value=read_local_log('fitting', 'sma_over_rs'))
    inclination = DoubleVar(value=read_local_log('fitting', 'inclination'))
    eccentricity = DoubleVar(value=read_local_log('fitting', 'eccentricity'))
    periastron = DoubleVar(value=read_local_log('fitting', 'periastron'))
    target_ra_dec = StringVar(root, value=read_local_log('fitting', 'target_ra_dec'))
    observer = StringVar(value=read_local_log('fitting', 'observer'))
    observatory = StringVar(value=read_local_log('fitting', 'observatory'))
    telescope = StringVar(value=read_local_log('fitting', 'telescope'))
    camera = StringVar(value=read_local_log('fitting', 'camera'))
    phot_filter = StringVar(value=read_local_log('fitting', 'phot_filter'))

    # set progress variables, useful for updating the window

    update_preview = BooleanVar(root, value=True)
    update_planet = BooleanVar(root, value=False)
    running = BooleanVar(root, value=False)

    # create widgets
    combostyle = ttk.Style()
    combostyle.theme_create('combostyle', parent='alt',
                            settings={'TCombobox': {'configure':
                                                    {'selectbackground': 'white',
                                                     'fieldbackground': 'white',
                                                     'background': 'white'}}})
    combostyle.theme_use('combostyle')

    light_curve_file_label = Label(root, text='Light-curve file')
    light_curve_file_entry = ttk.Combobox(root, textvariable=light_curve_file, state='readonly', width=55)

    binning_label = Label(root, text='Binning')
    binning_entry = Entry(root, textvariable=binning, validate='key')
    binning_entry['validatecommand'] = (binning_entry.register(test_int_positive_non_zero_input), '%P', '%d')

    scatter_label = Label(root, text='Scatter limit')
    scatter_entry = Entry(root, textvariable=scatter, validate='key')
    scatter_entry['validatecommand'] = (scatter_entry.register(test_float_positive_input), '%P', '%d')

    iterations_label = Label(root, text='Iterations')
    iterations_entry = Entry(root, textvariable=iterations, validate='key')
    iterations_entry['validatecommand'] = (iterations_entry.register(test_int_positive_non_zero_input), '%P', '%d')

    burn_label = Label(root, text='Burned iterations')
    burn_entry = Entry(root, textvariable=burn, validate='key')
    burn_entry['validatecommand'] = (burn_entry.register(test_int_positive_non_zero_input), '%P', '%d')

    metallicity_label = Label(root, text='M* [Fe/H, dex]')
    metallicity_entry = Entry(root, textvariable=metallicity, validate='key')
    metallicity_entry['validatecommand'] = (metallicity_entry.register(test_float_input), '%P', '%d')

    temperature_label = Label(root, text='T* [K]')
    temperature_entry = Entry(root, textvariable=temperature, validate='key')
    temperature_entry['validatecommand'] = (temperature_entry.register(test_float_positive_input), '%P', '%d')

    logg_label = Label(root, text='log(g*) [cm/s^2]')
    logg_entry = Entry(root, textvariable=logg, validate='key')
    logg_entry['validatecommand'] = (logg_entry.register(test_float_positive_input), '%P', '%d')

    period_label = Label(root, text='Period [days]')
    period_entry = Entry(root, textvariable=period, validate='key')
    period_entry['validatecommand'] = (period_entry.register(test_float_positive_input), '%P', '%d')

    mid_time_label = Label(root, text='Mid-time [days, BJD_TDB]')
    mid_time_entry = Entry(root, textvariable=mid_time, validate='key')
    mid_time_entry['validatecommand'] = (mid_time_entry.register(test_float_positive_input), '%P', '%d')

    rp_over_rs_label = Label(root, text='Rp/Rs')
    rp_over_rs_entry = Entry(root, textvariable=rp_over_rs, validate='key')
    rp_over_rs_entry['validatecommand'] = (rp_over_rs_entry.register(test_float_positive_input), '%P', '%d')

    sma_over_rs_label = Label(root, text='a/Rs')
    sma_over_rs_entry = Entry(root, textvariable=sma_over_rs, validate='key')
    sma_over_rs_entry['validatecommand'] = (sma_over_rs_entry.register(test_float_positive_input), '%P', '%d')

    inclination_label = Label(root, text='Inclination [deg]')
    inclination_entry = Entry(root, textvariable=inclination, validate='key')
    inclination_entry['validatecommand'] = (inclination_entry.register(test_float_positive_input), '%P', '%d')

    eccentricity_label = Label(root, text='Eccentricity')
    eccentricity_entry = Entry(root, textvariable=eccentricity, validate='key')
    eccentricity_entry['validatecommand'] = (eccentricity_entry.register(test_float_positive_input), '%P', '%d')

    periastron_label = Label(root, text='Periastron [deg]')
    periastron_entry = Entry(root, textvariable=periastron, validate='key')
    periastron_entry['validatecommand'] = (periastron_entry.register(test_float_positive_input), '%P', '%d')

    target_ra_dec_label = Label(root, text='Planet RA DEC\n(hh:mm:ss +/-dd:mm:ss)')
    target_ra_dec_entry = Entry(root, textvariable=target_ra_dec)
    target_ra_dec_test = Label(root, text=' ')

    observer_label = Label(root, text='Observer')
    observer_entry = Entry(root, textvariable=observer)

    observatory_label = Label(root, text='Observatory')
    observatory_entry = Entry(root, textvariable=observatory)

    telescope_label = Label(root, text='Telescope')
    telescope_entry = Entry(root, textvariable=telescope)

    camera_label = Label(root, text='Camera')
    camera_entry = Entry(root, textvariable=camera)

    planet_label = Label(root, text='Planet')
    planet_entry = ttk.Combobox(root, textvariable=planet, state='readonly', width=17)
    planet_search_entry = Entry(root, textvariable=planet_search)

    phot_filter_label = Label(root, text='Filter')
    phot_filter_entry = ttk.Combobox(root, textvariable=phot_filter, state='readonly', width=17)
    phot_filter_entry['values'] = tuple([ff for ff in filter_map])

    show_preview_button = Button(root, text='Show Preview')

    return_to_photometry_button = Button(root, text='RETURN TO PHOTOMETRY')

    return_to_reduction_button = Button(root, text='RETURN TO REDUCTION')

    fitting_button = Button(root, text='RUN FITTING', bg='green')

    my_profile_button = Button(root, text='MY PROFILE')

    # define additional windows

    show_preview_window = AddOnWindow('Preview', None, None, 1)
    f = Figure()
    f.patch.set_facecolor('white')
    ax1 = f.add_subplot(211)
    ax2 = f.add_subplot(212)
    canvas = FigureCanvasTkAgg(f, show_preview_window.root)
    canvas.get_tk_widget().pack()
    NavigationToolbar2TkAgg(canvas, show_preview_window.root)

    my_profile_window = AddOnWindow('My Profile', 2, 1.1, 1)

    core_headers = {ff.split(':')[0]: read_log_profile(ff.split(':')[0])
                    for ff in open(log_profile_file, 'r').readlines()}

    local_headers = yaml.load(open(local_log_profile_file, 'r'), Loader=yaml.SafeLoader)

    variables = {}
    labels = {}
    entries = {}
    for row, header in enumerate(core_headers):
        if header in local_headers:
            variables[header] = StringVar(my_profile_window.root, value=local_headers[header])
            labels[header] = Label(my_profile_window.root, text=header)
            entries[header] = Entry(my_profile_window.root, textvariable=variables[header])
        else:
            variables[header] = StringVar(my_profile_window.root, value=core_headers[header])
            labels[header] = Label(my_profile_window.root, text=header)
            entries[header] = Entry(my_profile_window.root, textvariable=variables[header])

    def update_headers():
        new_local_profile = {}
        for header2 in variables:
            new_local_profile[header2] = variables[header2].get()
        yaml.dump(new_local_profile, open(local_log_profile_file, 'w'), default_flow_style=False)
        update_window(None)

    update_headers_button = Button(my_profile_window.root, text='UPDATE')
    update_headers_button['command'] = update_headers

    stucture = [[], [[update_headers_button, 2]], []]
    for header in list(core_headers.keys())[:int(len(list(core_headers.keys())) / 2)]:
        stucture.append([[labels[header], 1], [entries[header], 2]])

    stucture.append([])

    for jj, header in enumerate(list(core_headers.keys())[int(len(list(core_headers.keys())) / 2):]):
        stucture[3 + jj].append([labels[header], 3])
        stucture[3 + jj].append([entries[header], 4])

    stucture.append([])

    setup_window(my_profile_window.root, stucture)

    # define the function that updates the window

    def update_window(entry):

        if not entry:
            pass

        planet_entry.selection_clear()
        phot_filter_entry.selection_clear()
        light_curve_file_entry.selection_clear()

        all_files = (glob.glob(os.path.join('{0}*'.format(photometry_directory), light_curve_aperture_file)) +
                     glob.glob(os.path.join('{0}*'.format(photometry_directory), light_curve_gauss_file)))
        all_files.sort()
        light_curve_file_entry['values'] = tuple(all_files)

        if planet.get() == 'Choose Planet':

            try:
                catalogue_planets = []

                ra_target, dec_target = read_local_log('photometry', 'target_ra_dec').split(' ')
                target = plc.Target(plc.Hours(ra_target), plc.Degrees(dec_target))

                for catalogue_planet in catalogue.planets:
                    if not np.isnan(catalogue_planet.system.dec):
                        test_target = plc.Target(plc.Degrees(catalogue_planet.system.ra.deg),
                                                 plc.Degrees(catalogue_planet.system.dec.deg))
                        catalogue_planets.append([test_target.distance_from_target(target).rad,
                                                  catalogue_planet.name])
                catalogue_planets.sort()

                planet.set(catalogue_planets[0][1])
                planet_search.set(catalogue_planets[0][1])

                parameters = plc.find_oec_parameters(planet.get(), catalogue=catalogue)
                coordinates = plc.find_oec_coordinates(planet.get(), catalogue=catalogue, output='str')
                logg.set(parameters[1])
                temperature.set(parameters[2])
                metallicity.set(parameters[3])
                rp_over_rs.set(parameters[4])
                period.set(parameters[6])
                sma_over_rs.set(parameters[7])
                eccentricity.set(parameters[8])
                inclination.set(parameters[9])
                periastron.set(parameters[10])
                mid_time.set(parameters[11])
                target_ra_dec.set(coordinates)

            except:
                pass

        if running.get():

            light_curve_file_entry['state'] = DISABLED
            binning_entry['state'] = DISABLED
            scatter_entry['state'] = DISABLED
            iterations_entry['state'] = DISABLED
            burn_entry['state'] = DISABLED
            metallicity_entry['state'] = DISABLED
            temperature_entry['state'] = DISABLED
            logg_entry['state'] = DISABLED
            phot_filter_entry['state'] = DISABLED
            period_entry['state'] = DISABLED
            mid_time_entry['state'] = DISABLED
            rp_over_rs_entry['state'] = DISABLED
            sma_over_rs_entry['state'] = DISABLED
            inclination_entry['state'] = DISABLED
            eccentricity_entry['state'] = DISABLED
            periastron_entry['state'] = DISABLED
            target_ra_dec_entry['state'] = DISABLED
            observer_entry['state'] = DISABLED
            observatory_entry['state'] = DISABLED
            telescope_entry['state'] = DISABLED
            camera_entry['state'] = DISABLED
            planet_entry['state'] = DISABLED
            planet_search_entry['state'] = DISABLED
            return_to_photometry_button['state'] = DISABLED
            return_to_reduction_button['state'] = DISABLED
            fitting_button['state'] = DISABLED
            my_profile_button['state'] = DISABLED
            show_preview_button['state'] = DISABLED

        elif not os.path.isfile(light_curve_file.get()):

            light_curve_file_entry['state'] = NORMAL
            binning_entry['state'] = DISABLED
            scatter_entry['state'] = DISABLED
            iterations_entry['state'] = DISABLED
            burn_entry['state'] = DISABLED
            metallicity_entry['state'] = DISABLED
            temperature_entry['state'] = DISABLED
            logg_entry['state'] = DISABLED
            phot_filter_entry['state'] = DISABLED
            period_entry['state'] = DISABLED
            mid_time_entry['state'] = DISABLED
            rp_over_rs_entry['state'] = DISABLED
            sma_over_rs_entry['state'] = DISABLED
            inclination_entry['state'] = DISABLED
            eccentricity_entry['state'] = DISABLED
            periastron_entry['state'] = DISABLED
            target_ra_dec_entry['state'] = DISABLED
            observer_entry['state'] = DISABLED
            observatory_entry['state'] = DISABLED
            telescope_entry['state'] = DISABLED
            camera_entry['state'] = DISABLED
            planet_entry['state'] = DISABLED
            planet_search_entry['state'] = DISABLED
            return_to_photometry_button['state'] = DISABLED
            return_to_reduction_button['state'] = DISABLED
            fitting_button['state'] = DISABLED
            my_profile_button['state'] = NORMAL
            show_preview_button['state'] = DISABLED

        else:

            light_curve_file_entry['state'] = 'readonly'
            binning_entry['state'] = NORMAL
            scatter_entry['state'] = NORMAL
            iterations_entry['state'] = NORMAL
            burn_entry['state'] = NORMAL
            metallicity_entry['state'] = NORMAL
            temperature_entry['state'] = NORMAL
            logg_entry['state'] = NORMAL
            phot_filter_entry['state'] = 'readonly'
            period_entry['state'] = NORMAL
            mid_time_entry['state'] = NORMAL
            rp_over_rs_entry['state'] = NORMAL
            sma_over_rs_entry['state'] = NORMAL
            inclination_entry['state'] = NORMAL
            eccentricity_entry['state'] = NORMAL
            periastron_entry['state'] = NORMAL
            target_ra_dec_entry['state'] = NORMAL
            observer_entry['state'] = NORMAL
            observatory_entry['state'] = NORMAL
            telescope_entry['state'] = NORMAL
            camera_entry['state'] = NORMAL
            planet_entry['state'] = 'readonly'
            planet_search_entry['state'] = NORMAL
            my_profile_button['state'] = NORMAL
            show_preview_button['state'] = NORMAL

            if isinstance(catalogue.searchPlanet(planet_search.get()), list):
                planet_entry['values'] = tuple([ppp.name for ppp in catalogue.searchPlanet(planet_search.get())])
            elif catalogue.searchPlanet(planet_search.get()):
                planet_entry['values'] = tuple([catalogue.searchPlanet(planet_search.get()).name])
            else:
                planet_entry['values'] = tuple([])

            if update_planet.get():

                try:

                    parameters = plc.find_oec_parameters(planet.get(), catalogue=catalogue)
                    coordinates = plc.find_oec_coordinates(planet.get(), catalogue=catalogue, output='str')
                    planet_search.set(planet.get())
                    logg.set(parameters[1])
                    temperature.set(parameters[2])
                    metallicity.set(parameters[3])
                    rp_over_rs.set(parameters[4])
                    period.set(parameters[6])
                    sma_over_rs.set(parameters[7])
                    eccentricity.set(parameters[8])
                    inclination.set(parameters[9])
                    periastron.set(parameters[10])
                    mid_time.set(parameters[11])
                    target_ra_dec.set(coordinates)

                except:
                    pass

            if target_ra_dec.get() == 'hh:mm:ss +dd:mm:ss':
                try:
                    coordinates = plc.find_oec_coordinates(planet.get(), catalogue=catalogue, output='str')
                    target_ra_dec.set(coordinates)
                except:
                    pass

            if phot_filter.get() == 'default':
                for key in read_local_log_profile('filter_key').split(','):
                    check_filter = test_fits_keyword(observation_files.get(), key)
                    if check_filter[0]:
                        if check_filter[2] in filter_map:
                            phot_filter.set(check_filter[2])
                            break
                if phot_filter.get() == 'default':
                    phot_filter.set(read_local_log_profile('filter'))

            if telescope.get() == 'default':
                for key in read_local_log_profile('telescope_key').split(','):
                    check_telescope = test_fits_keyword(observation_files.get(), key)
                    if check_telescope[0]:
                        telescope.set(check_telescope[2])
                        break
                if telescope.get() == 'default':
                    telescope.set(read_local_log_profile('telescope'))

            if camera.get() == 'default':
                for key in read_local_log_profile('camera_key').split(','):
                    check_camera = test_fits_keyword(observation_files.get(), key)
                    if check_camera[0]:
                        camera.set(check_camera[2])
                        break
                if camera.get() == 'default':
                    camera.set(read_local_log_profile('camera'))

            if observer.get() == 'default':
                for key in read_local_log_profile('observer_key').split(','):
                    check_observer = test_fits_keyword(observation_files.get(), key)
                    if check_observer[0]:
                        observer.set(check_observer[2])
                        break
                if observer.get() == 'default':
                    observer.set(read_local_log_profile('observer'))

            if observatory.get() == 'default':
                for key in read_local_log_profile('observatory_key').split(','):
                    check_observatory = test_fits_keyword(observation_files.get(), key)
                    if check_observatory[0]:
                        observatory.set(check_observatory[2])
                        break
                if observatory.get() == 'default':
                    observatory.set(read_local_log_profile('observatory'))

            enable_buttons = True

            if not os.path.isfile(light_curve_file.get()):
                enable_buttons = False

            check_ra_dec = test_coordinates(target_ra_dec_entry.get())
            target_ra_dec_test.configure(text=check_ra_dec[1])

            if not check_ra_dec[0]:
                enable_buttons = False

            for input_entry in [binning_entry, scatter_entry, iterations_entry, burn_entry, phot_filter_entry,
                                metallicity_entry, temperature_entry, logg_entry, period_entry, mid_time_entry,
                                rp_over_rs_entry, sma_over_rs_entry, inclination_entry, eccentricity_entry,
                                periastron_entry, target_ra_dec_entry]:

                if len(str(input_entry.get())) == 0:
                    enable_buttons = False

            if enable_buttons:
                fitting_button['state'] = NORMAL
                show_preview_button['state'] = NORMAL

            else:
                fitting_button['state'] = DISABLED
                show_preview_button['state'] = DISABLED

            return_to_photometry_button['state'] = NORMAL
            return_to_reduction_button['state'] = NORMAL

        try:

            if update_preview.get():

                light_curve = np.loadtxt(light_curve_file.get(), unpack=True)

                if binning.get() > 1:
                    start = len(light_curve[0]) - (len(light_curve[0]) // binning.get()) * binning.get()
                    light_curve_0 = np.mean(np.reshape(light_curve[0][start:],
                                                       (light_curve[0].size // binning.get(), binning.get())), 1)
                    light_curve_1 = np.mean(np.reshape(light_curve[1][start:],
                                                       (light_curve[1].size // binning.get(), binning.get())), 1)
                else:
                    light_curve_0 = light_curve[0]
                    light_curve_1 = light_curve[1]

                light_curve_0 = light_curve_0[np.where(~np.isnan(light_curve_1))]
                light_curve_1 = light_curve_1[np.where(~np.isnan(light_curve_1))]

                moving_average = []
                for i in range(-5, 6):
                    moving_average.append(np.roll(light_curve_1, i))

                median = np.median(moving_average, 0)
                med = np.median([np.abs(ff - median) for ff in moving_average], 0)

                flag = np.where((np.abs(light_curve_1 - median) < scatter.get() * med))[0]

                limb_darkening_coefficients = plc.clablimb('claret', logg.get(), max(4000, temperature.get()),
                                                           metallicity.get(), filter_map[phot_filter.get()])

                data_delta_t = light_curve_0[flag] - light_curve_0[flag][0]

                def mcmc_f(inputs, detrend_zero, detrend_one, detrend_two, model_rp_over_rs, model_mid_time):

                    if inputs:
                        detrend = detrend_zero * (1 + detrend_one * data_delta_t +
                                                  detrend_two * data_delta_t * data_delta_t)
                        transit_model = plc.transit('claret', limb_darkening_coefficients, model_rp_over_rs,
                                                    period.get(), sma_over_rs.get(), eccentricity.get(),
                                                    inclination.get(), periastron.get(),
                                                    mid_time.get() + model_mid_time,
                                                    time_array=light_curve_0[flag])

                        return detrend * transit_model

                def independent_f(detrend_zero, detrend_one, detrend_two, model_rp_over_rs, model_mid_time):

                    detrend = detrend_zero * (1 + detrend_one * data_delta_t +
                                              detrend_two * data_delta_t * data_delta_t)
                    transit_model = plc.transit('claret', limb_darkening_coefficients, model_rp_over_rs, period.get(),
                                                sma_over_rs.get(), eccentricity.get(), inclination.get(),
                                                periastron.get(), mid_time.get() + model_mid_time,
                                                time_array=light_curve_0[flag])

                    return detrend, transit_model

                popt, pcov = curve_fit(mcmc_f, [1], light_curve_1[flag],
                                       p0=[np.mean(light_curve_1[flag]), 1, 1, rp_over_rs.get(), 0])

                fit_detrend, fit_transit_model = independent_f(*popt)

                predicted_transit_model = plc.transit('claret', limb_darkening_coefficients, rp_over_rs.get(),
                                                      period.get(), sma_over_rs.get(), eccentricity.get(),
                                                      inclination.get(), periastron.get(), mid_time.get(),
                                                      time_array=light_curve_0[flag])

                new_mid_time = (mid_time.get()
                                + round((np.mean(light_curve_0) - mid_time.get()) / period.get()) * period.get()
                                + popt[-1])

                phase = np.array((light_curve_0 - new_mid_time) / period.get())

                ax1.cla()
                ax2.cla()

                ax1.plot(phase, light_curve_1, 'ro', ms=3, mec='r')
                ax1.plot(phase[flag], light_curve_1[flag], 'ko', ms=3)
                ax1.plot(phase[flag], fit_detrend * fit_transit_model, 'r-')
                ax1.set_yticks(ax1.get_yticks()[1:])
                ax1.tick_params(labelbottom=False)
                ax1.set_ylabel(r'$\mathrm{relative} \ \mathrm{flux}$', fontsize=20)

                ax2.plot(phase[flag], light_curve_1[flag] / fit_detrend, 'ko', ms=3)
                ax2.plot(phase[flag], fit_transit_model, 'r-')
                ax2.plot(phase[flag], predicted_transit_model, 'c-')
                ax2.set_title('{0:.1e}'.format(np.std(light_curve_1[flag] - fit_transit_model)))
                ax2.set_ylabel(r'$\mathrm{normalised} \ \mathrm{flux}$', fontsize=20)
                ax2.set_xlabel(r'$\mathrm{phase}$', fontsize=20)

        except:
            pass

        canvas.draw()

    update_window(None)

    # define actions for the different buttons, including calls to the function that updates the window

    def choose_planet(entry):

        if not entry:
            return 0

        update_planet.set(True)
        update_window(None)
        update_planet.set(False)

    def return_to_photometry():

        write_local_log('fitting', light_curve_file.get(), 'light_curve_file')
        write_local_log('fitting', planet_search.get(), 'planet_search')
        write_local_log('fitting', planet.get(), 'planet')
        write_local_log('fitting', binning.get(), 'binning')
        write_local_log('fitting', scatter.get(), 'scatter')
        write_local_log('fitting', iterations.get(), 'iterations')
        write_local_log('fitting', burn.get(), 'burn')
        write_local_log('fitting', metallicity.get(), 'metallicity')
        write_local_log('fitting', temperature.get(), 'temperature')
        write_local_log('fitting', logg.get(), 'logg')
        write_local_log('fitting', period.get(), 'period')
        write_local_log('fitting', mid_time.get(), 'mid_time')
        write_local_log('fitting', rp_over_rs.get(), 'rp_over_rs')
        write_local_log('fitting', sma_over_rs.get(), 'sma_over_rs')
        write_local_log('fitting', inclination.get(), 'inclination')
        write_local_log('fitting', eccentricity.get(), 'eccentricity')
        write_local_log('fitting', periastron.get(), 'periastron')
        write_local_log('fitting', target_ra_dec.get(), 'target_ra_dec')
        write_local_log('fitting', observer.get(), 'observer')
        write_local_log('fitting', observatory.get(), 'observatory')
        write_local_log('fitting', telescope.get(), 'telescope')
        write_local_log('fitting', camera.get(), 'camera')
        write_local_log('fitting', phot_filter.get(), 'phot_filter')

        run.run_from_reduction = False
        run.run_from_photometry = True
        run.run_from_fitting = False

        show_preview_window.close()
        my_profile_window.close()
        root.destroy()

    def return_to_reduction():

        write_local_log('fitting', light_curve_file.get(), 'light_curve_file')
        write_local_log('fitting', planet_search.get(), 'planet_search')
        write_local_log('fitting', planet.get(), 'planet')
        write_local_log('fitting', binning.get(), 'binning')
        write_local_log('fitting', scatter.get(), 'scatter')
        write_local_log('fitting', iterations.get(), 'iterations')
        write_local_log('fitting', burn.get(), 'burn')
        write_local_log('fitting', metallicity.get(), 'metallicity')
        write_local_log('fitting', temperature.get(), 'temperature')
        write_local_log('fitting', logg.get(), 'logg')
        write_local_log('fitting', period.get(), 'period')
        write_local_log('fitting', mid_time.get(), 'mid_time')
        write_local_log('fitting', rp_over_rs.get(), 'rp_over_rs')
        write_local_log('fitting', sma_over_rs.get(), 'sma_over_rs')
        write_local_log('fitting', inclination.get(), 'inclination')
        write_local_log('fitting', eccentricity.get(), 'eccentricity')
        write_local_log('fitting', periastron.get(), 'periastron')
        write_local_log('fitting', target_ra_dec.get(), 'target_ra_dec')
        write_local_log('fitting', observer.get(), 'observer')
        write_local_log('fitting', observatory.get(), 'observatory')
        write_local_log('fitting', telescope.get(), 'telescope')
        write_local_log('fitting', camera.get(), 'camera')
        write_local_log('fitting', phot_filter.get(), 'phot_filter')

        run.run_from_reduction = True
        run.run_from_photometry = False
        run.run_from_fitting = False
        show_preview_window.close()
        my_profile_window.close()
        root.destroy()

    def fitting():

        running.set(True)
        update_window(None)

        write_local_log('fitting', light_curve_file.get(), 'light_curve_file')
        write_local_log('fitting', planet_search.get(), 'planet_search')
        write_local_log('fitting', planet.get(), 'planet')
        write_local_log('fitting', binning.get(), 'binning')
        write_local_log('fitting', scatter.get(), 'scatter')
        write_local_log('fitting', iterations.get(), 'iterations')
        write_local_log('fitting', burn.get(), 'burn')
        write_local_log('fitting', metallicity.get(), 'metallicity')
        write_local_log('fitting', temperature.get(), 'temperature')
        write_local_log('fitting', logg.get(), 'logg')
        write_local_log('fitting', period.get(), 'period')
        write_local_log('fitting', mid_time.get(), 'mid_time')
        write_local_log('fitting', rp_over_rs.get(), 'rp_over_rs')
        write_local_log('fitting', sma_over_rs.get(), 'sma_over_rs')
        write_local_log('fitting', inclination.get(), 'inclination')
        write_local_log('fitting', eccentricity.get(), 'eccentricity')
        write_local_log('fitting', periastron.get(), 'periastron')
        write_local_log('fitting', target_ra_dec.get(), 'target_ra_dec')
        write_local_log('fitting', observer.get(), 'observer')
        write_local_log('fitting', observatory.get(), 'observatory')
        write_local_log('fitting', telescope.get(), 'telescope')
        write_local_log('fitting', camera.get(), 'camera')
        write_local_log('fitting', phot_filter.get(), 'phot_filter')

        ftr_fitting()

        running.set(False)
        update_window(None)

    # connect widgets to functions

    light_curve_file_entry.bind('<<ComboboxSelected>>', update_window)
    planet_entry.bind('<<ComboboxSelected>>', choose_planet)
    planet_search_entry.bind(sequence='<KeyRelease>', func=update_window)
    binning_entry.bind(sequence='<KeyRelease>', func=update_window)
    scatter_entry.bind(sequence='<KeyRelease>', func=update_window)
    iterations_entry.bind(sequence='<KeyRelease>', func=update_window)
    burn_entry.bind(sequence='<KeyRelease>', func=update_window)
    phot_filter_entry.bind('<<ComboboxSelected>>', update_window)
    metallicity_entry.bind(sequence='<KeyRelease>', func=update_window)
    temperature_entry.bind(sequence='<KeyRelease>', func=update_window)
    logg_entry.bind(sequence='<KeyRelease>', func=update_window)
    period_entry.bind(sequence='<KeyRelease>', func=update_window)
    mid_time_entry.bind(sequence='<KeyRelease>', func=update_window)
    rp_over_rs_entry.bind(sequence='<KeyRelease>', func=update_window)
    sma_over_rs_entry.bind(sequence='<KeyRelease>', func=update_window)
    inclination_entry.bind(sequence='<KeyRelease>', func=update_window)
    eccentricity_entry.bind(sequence='<KeyRelease>', func=update_window)
    periastron_entry.bind(sequence='<KeyRelease>', func=update_window)
    target_ra_dec_entry.bind(sequence='<KeyRelease>', func=update_window)
    show_preview_button['command'] = show_preview_window.show
    return_to_photometry_button['command'] = return_to_photometry
    return_to_reduction_button['command'] = return_to_reduction
    fitting_button['command'] = fitting
    my_profile_button['command'] = my_profile_window.show

    # setup window

    Btn = Button(root, text="HOPS UPDATES &\nUSER MANUAL", command=openweb)

    photo = PhotoImage(file=holomon_logo)
    logo_label = Label(root, image=photo)
    window_label = Label(root, text='Fitting')
    created_by_label = Label(root, text=read_log('windows', 'created_by').replace(',', '\n'))

    setup_window(root, [
        [],
        [[logo_label, 0, 1, 6], [window_label, 1, 4, 1, 'title']],
        [],
        [[light_curve_file_label, 1], [light_curve_file_entry, 2, 3, 1]],
        [[binning_label, 1], [binning_entry, 2], [scatter_label, 3], [scatter_entry, 4]],
        [],
        [[phot_filter_label, 1], [phot_filter_entry, 2], [telescope_label, 3], [telescope_entry, 4]],
        [[created_by_label, 0, 1, 3], [camera_label, 1], [camera_entry, 2], [observatory_label, 3],
         [observatory_entry, 4]],
        [[observer_label, 3], [observer_entry, 4]],
        [],
        [[my_profile_button, 0], [planet_label, 1], [planet_search_entry, 2], [target_ra_dec_label, 3], [target_ra_dec_entry, 4]],
        [[Btn, 0, 1, 2], [planet_entry, 2], [target_ra_dec_test, 4]],
        [],
        [[period_label, 1], [period_entry, 2], [metallicity_label, 3], [metallicity_entry, 4]],
        [[mid_time_label, 1], [mid_time_entry, 2], [temperature_label, 3], [temperature_entry, 4]],
        [[rp_over_rs_label, 1], [rp_over_rs_entry, 2], [logg_label, 3], [logg_entry, 4]],
        [[sma_over_rs_label, 1], [sma_over_rs_entry, 2]],
        [[inclination_label, 1], [inclination_entry, 2], [iterations_label, 3], [iterations_entry, 4]],
        [[eccentricity_label, 1], [eccentricity_entry, 2], [burn_label, 3], [burn_entry, 4]],
        [[periastron_label, 1], [periastron_entry, 2]],
        [[show_preview_button, 2], [fitting_button, 3, 2]],
        [[return_to_photometry_button, 3, 2]],
        [[return_to_reduction_button, 3, 2]],
        []
    ])

    # finalise and show  window

    finalise_window(root, position=5)
    show_preview_window.mainloop()
    my_profile_window.mainloop()
    root.mainloop()


class Run:

    def __init__(self):
        self.run_from_reduction = True
        self.run_from_photometry = False
        self.run_from_fitting = False
        self.exit = False

    def run(self):

        if not self.exit and self.run_from_reduction:
            reduction_alignment_window(self)
            return None

        if not self.exit and self.run_from_photometry:
            photometry_window(self)
            return None

        if not self.exit and self.run_from_fitting:
            fitting_window(self)
            return None


def run_app():
    print('Loading... Please wait for the main window to appear.')

    run = Run()

    while not run.exit:
        run.run()
