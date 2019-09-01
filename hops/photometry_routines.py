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


def finalise_window(window, center=True):

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
    window.after_idle(window.attributes, '-topmost', 0)
    window.deiconify()


def photometry():

    # get variables

    reduction_directory = read_local_log('pipeline', 'reduction_directory')
    light_curve_aperture_file = read_local_log('pipeline', 'light_curve_aperture_file')
    photometry_directory = read_local_log('pipeline', 'photometry_directory')
    photometry_file = read_local_log('pipeline', 'photometry_file')
    light_curve_gauss_file = read_local_log('pipeline', 'light_curve_gauss_file')
    results_figure = read_local_log('pipeline', 'results_figure')
    fov_figure = read_local_log('pipeline', 'fov_figure')
    mean_key = read_local_log('pipeline_keywords', 'mean_key')
    std_key = read_local_log('pipeline_keywords', 'std_key')
    align_x0_key = read_local_log('pipeline_keywords', 'align_x0_key')
    align_y0_key = read_local_log('pipeline_keywords', 'align_y0_key')
    align_u0_key = read_local_log('pipeline_keywords', 'align_u0_key')
    # exposure_time_key = read_log('pipeline_keywords', 'exposure_time_key')
    observation_date_key = read_local_log('pipeline_keywords', 'observation_date_key')
    observation_time_key = read_local_log('pipeline_keywords', 'observation_time_key')
    star_std = read_local_log('alignment', 'star_std')
    search_window_std = read_local_log('alignment', 'search_window_std')
    target_ra_dec = read_local_log('photometry', 'target_ra_dec')
    sky_inner_aperture = read_local_log('photometry', 'sky_inner_aperture')
    sky_outer_aperture = read_local_log('photometry', 'sky_outer_aperture')
    max_comparisons = read_local_log('photometry', 'max_comparisons')
    bin_fits = int(read_local_log('reduction', 'bin_fits'))
    burn_limit = int(read_local_log('alignment', 'burn_limit')) * bin_fits * bin_fits
    targets_r_position = [read_local_log('photometry', 'target_r_position')]
    targets_u_position = [read_local_log('photometry', 'target_u_position')]
    targets_aperture = [read_local_log('photometry', 'target_aperture')]
    for comparison in range(max_comparisons):
        targets_r_position.append(read_local_log('photometry', 'comparison_{0}_r_position'.format(comparison)))
        targets_u_position.append(read_local_log('photometry', 'comparison_{0}_u_position'.format(comparison)))
        targets_aperture.append(read_local_log('photometry', 'comparison_{0}_aperture'.format(comparison)))

    science = find_fits_files(os.path.join(reduction_directory, '*'))

    root = Tk()

    exit_var = BooleanVar(value=False)

    def break_and_exit():
        exit_var.set(True)

    initialise_window(root, exit_command=break_and_exit)

    label1 = Label(root, text='PHOTOMETRY')
    label2 = Label(root, text='FILE:')
    label3 = Label(root, text=' ')
    label4 = Label(root, text='COMPLETE:')
    label5 = Label(root, text=' ',)
    label6 = Label(root, text='TIME LEFT:')
    label7 = Label(root, text=' ')

    setup_window(root, [
        [[label1, 1, 2]],
        [[label2, 1], [label3, 2]],
        [[label4, 1], [label5, 2]],
        [[label6, 1], [label7, 2]],
    ])

    gauss_targets_files = []
    gauss_targets_jd = []
    gauss_targets_x_position = []
    gauss_targets_y_position = []
    gauss_targets_x_std = []
    gauss_targets_y_std = []
    gauss_targets_flux = []
    gauss_targets_sky = []
    apperture_targets_files = []
    apperture_targets_jd = []
    apperture_targets_x_position = []
    apperture_targets_y_position = []
    apperture_targets_flux = []
    apperture_targets_sky = []

    # TODO exclude live points

    # for each science_file
    percent = 0
    lt0 = time.time()
    for counter, science_file in enumerate(science):

        fits = pf.open(science_file)

        if fits[1].header[align_x0_key]:

            # calculate heliocentric julian date

            if observation_date_key == observation_time_key:
                local_time = ' '.join(fits[1].header[observation_date_key].split('T'))
            else:
                local_time = ' '.join([fits[1].header[observation_date_key], fits[1].header[observation_time_key]])

            julian_date = plc.ut_to_jd(local_time)

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
            gauss_targets_sky_test = []

            for target in range(max_comparisons + 1):

                if targets_aperture[target] > 0:

                    norm, floor, x_mean, y_mean, x_std, y_std = \
                        fit_2d_gauss_point(fits[1].data,
                                           predicted_x_mean=(ref_x_position + targets_r_position[target] *
                                                             np.cos(ref_u_position + targets_u_position[target])),
                                           predicted_y_mean=(ref_y_position + targets_r_position[target] *
                                                             np.sin(ref_u_position + targets_u_position[target])),
                                           search_window=search_window_std * 3 * star_std, stde=star_std, snr_lim=False)

                    gauss_targets_x_position_test.append(x_mean)
                    gauss_targets_y_position_test.append(y_mean)
                    gauss_targets_x_std_test.append(x_std)
                    gauss_targets_y_std_test.append(y_std)
                    gauss_targets_flux_test.append(2 * np.pi * norm * x_std * y_std)
                    gauss_targets_sky_test.append(floor)

            if np.nan in gauss_targets_x_position_test:

                gauss_targets_files_test = []
                gauss_targets_jd_test = []
                gauss_targets_x_position_test = []
                gauss_targets_y_position_test = []
                gauss_targets_x_std_test = []
                gauss_targets_y_std_test = []
                gauss_targets_flux_test = []
                gauss_targets_sky_test = []

            gauss_targets_files += gauss_targets_files_test
            gauss_targets_jd += gauss_targets_jd_test
            gauss_targets_x_position += gauss_targets_x_position_test
            gauss_targets_y_position += gauss_targets_y_position_test
            gauss_targets_x_std += gauss_targets_x_std_test
            gauss_targets_y_std += gauss_targets_y_std_test
            gauss_targets_flux += gauss_targets_flux_test
            gauss_targets_sky += gauss_targets_sky_test

            apperture_targets_files_test = [science_file]
            apperture_targets_jd_test = [julian_date]
            apperture_targets_x_position_test = []
            apperture_targets_y_position_test = []
            apperture_targets_flux_test = []
            apperture_targets_sky_test = []

            skip = False
            for target in range(max_comparisons + 1):

                if targets_aperture[target] > 0:

                    x_mean = (ref_x_position + targets_r_position[target] *
                              np.cos(ref_u_position + targets_u_position[target]))
                    y_mean = (ref_y_position + targets_r_position[target] *
                              np.sin(ref_u_position + targets_u_position[target]))

                    centroids = find_centroids(fits[1].data,
                                               x_low=int(x_mean - search_window_std * star_std),
                                               x_upper=int(x_mean + search_window_std * star_std + 1),
                                               y_low=int(y_mean - search_window_std * star_std),
                                               y_upper=int(y_mean + search_window_std * star_std + 1),
                                               x_centre=int(x_mean), y_centre=int(y_mean),
                                               mean=fits[1].header[mean_key], std=fits[1].header[std_key],
                                               star_std=star_std)

                    if len(centroids) == 0:
                        skip = True
                        break

                    x_mean = centroids[0][1]
                    y_mean = centroids[0][2]

                    flux_area = fits[1].data[int(y_mean) - targets_aperture[target]:
                                             int(y_mean) + targets_aperture[target] + 1,
                                             int(x_mean) - targets_aperture[target]:
                                             int(x_mean) + targets_aperture[target] + 1]
                    flux_pixels = (2 * targets_aperture[target] + 1) ** 2
                    flux = np.sum(flux_area)

                    sky_area_1 = int(sky_inner_aperture * targets_aperture[target])
                    sky_area_2 = int(sky_outer_aperture * targets_aperture[target])
                    fits[1].data[int(y_mean) - sky_area_1:int(y_mean) + sky_area_1 + 1,
                                 int(x_mean) - sky_area_1:int(x_mean) + sky_area_1 + 1] = 0
                    sky_area = fits[1].data[int(y_mean) - sky_area_2:int(y_mean) + sky_area_2 + 1,
                                            int(x_mean) - sky_area_2:int(x_mean) + sky_area_2 + 1]
                    sky_area = sky_area[np.where((sky_area > 0) &
                                                 (sky_area < fits[1].header[mean_key] + 3 * fits[1].header[
                                                     std_key]))]
                    sky = np.sum(sky_area)
                    sky_pixels = sky_area.size

                    apperture_targets_x_position_test.append(x_mean)
                    apperture_targets_y_position_test.append(y_mean)
                    apperture_targets_flux_test.append(flux - flux_pixels * sky / sky_pixels)
                    apperture_targets_sky_test.append(sky / sky_pixels)

            if not skip:
                apperture_targets_files += apperture_targets_files_test
                apperture_targets_jd += apperture_targets_jd_test
                apperture_targets_x_position += apperture_targets_x_position_test
                apperture_targets_y_position += apperture_targets_y_position_test
                apperture_targets_flux += apperture_targets_flux_test
                apperture_targets_sky += apperture_targets_sky_test
            else:
                print(science_file, ': Stars not found!')

        # counter

        new_percent = round(100 * (counter + 1) / float(len(science)), 1)
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
            finalise_window(root)

        root.update()

        if exit_var.get():
            break

        if counter + 1 == len(science):
            write_local_log('pipeline', True, 'photometry_complete')

    root.destroy()

    if not exit_var.get():

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
        targets_sky = np.swapaxes(np.reshape(gauss_targets_sky, (measurements_number, targets_number)), 0, 1)

        targets_results = [targets_jd] + (list(targets_x_position) + list(targets_y_position) +
                                          list(targets_x_std) + list(targets_y_std) + list(targets_flux) +
                                          list(targets_sky))

        np.savetxt(photometry_file.replace('.txt', '_g.txt'), np.swapaxes(targets_results, 0, 1))

        np.savetxt(light_curve_gauss_file,
                   np.swapaxes([targets_jd, targets_flux[0] / np.sum(targets_flux[1:], 0)], 0, 1))

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
        targets_sky = np.swapaxes(np.reshape(apperture_targets_sky, (measurements_number, targets_number)), 0, 1)

        targets_results = [targets_jd] + (list(targets_x_position) + list(targets_y_position) + list(targets_flux) +
                                          list(targets_sky))

        np.savetxt(photometry_file.replace('.txt', '_a.txt'), np.swapaxes(targets_results, 0, 1))

        np.savetxt(light_curve_aperture_file,
                   np.swapaxes([targets_jd, targets_flux[0] / np.sum(targets_flux[1:], 0)], 0, 1))

        aflux = targets_flux[0] / np.sum(targets_flux[1:], 0)
        diff = np.abs(aflux[1:] - aflux[:-1])
        apperture_scatter = np.std(diff)

        if not os.path.isdir(photometry_directory):
            os.mkdir(photometry_directory)
        else:
            fi = 2
            while os.path.isdir('{0}_{1}'.format(photometry_directory, str(fi))):
                fi += 1
            photometry_directory = '{0}_{1}'.format(photometry_directory, str(fi))
            os.mkdir(photometry_directory)

        root = Tk()

        if comparisons_number > 1:
            f = Figure()
            f.set_figwidth(7)
            f.set_figheight(0.8 * root.winfo_screenheight() / f.get_dpi())
            ax = f.add_subplot(comparisons_number + 1, 1, 1)
        else:
            f = Figure()
            ax = f.add_subplot(1, 1, 1)

        exit_var_2 = BooleanVar(value=False)

        def break_and_exit():
            exit_var_2.set(True)

        initialise_window(root, exit_command=break_and_exit)

        f.patch.set_facecolor('white')
        canvas = FigureCanvasTkAgg(f, root)
        canvas.get_tk_widget().pack()
        NavigationToolbar2TkAgg(canvas, root)

        ax.plot(targets_jd - targets_jd[0], targets_flux[0] / np.sum(targets_flux[1:], 0)
                / np.median(targets_flux[0] / np.sum(targets_flux[1:], 0)), 'ko', ms=3, label='Aperture')
        ax.plot(gauss_targets_jd - gauss_targets_jd[0], targets_gauss_flux[0] / np.sum(targets_gauss_flux[1:], 0)
                / np.median(targets_gauss_flux[0] / np.sum(targets_gauss_flux[1:], 0)), 'ro', ms=3, mec='r', label='PSF')
        ax.tick_params(labelbottom=False)
        ax.set_title(r'$\mathrm{Target}$')
        ax.legend()

        if comparisons_number > 1:
            for comp in range(comparisons_number):
                test_aperture_flux = list(targets_flux[1:])
                test_gauss_flux = list(targets_gauss_flux[1:])
                del test_aperture_flux[comp]
                del test_gauss_flux[comp]
                ax = f.add_subplot(comparisons_number + 1, 1, comp + 2)
                ax.plot(targets_jd - targets_jd[0],
                        (targets_flux[1:][comp] / np.sum(test_aperture_flux, 0)
                         / np.median(targets_flux[1:][comp] / np.sum(test_aperture_flux, 0))), 'ko', ms=3)
                ax.plot(gauss_targets_jd - gauss_targets_jd[0],
                        (targets_gauss_flux[1:][comp] / np.sum(test_gauss_flux, 0) /
                         np.median(targets_gauss_flux[1:][comp] / np.sum(test_gauss_flux, 0))),
                        'ro', ms=3, mec='r')
                ax.tick_params(labelbottom=False)
                ax.set_title(r'${0}{1}{2}$'.format('\mathrm{', "Comparison \, {0}".format(comp + 1), '}'))

        ax.tick_params(labelbottom=True)
        ax.set_xlabel(r'$\mathrm{\Delta t} \ \mathrm{[days]}$', fontsize=20)
        f.text(0.03, 0.5, r'$\mathrm{relative} \ \mathrm{flux}$', fontsize=20,
               ha='center', va='center', rotation='vertical')

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

        exposure_time_key = read_local_log('pipeline_keywords', 'exposure_time_key')
        fits = pf.open(science[np.random.randint(len(science))])
        exp_time = round(fits[1].header[exposure_time_key], 1)

        catalogue = plc.oec_catalogue()
        catalogue_planets = []
        ra_target, dec_target = ra_dec_string_to_deg(read_local_log('photometry', 'target_ra_dec'))
        for catalogue_planet in catalogue.planets:
            if not np.isnan(catalogue_planet.system.dec):
                catalogue_planets.append([np.sqrt((catalogue_planet.system.dec.deg - dec_target) ** 2
                                                  + (catalogue_planet.system.ra.deg - ra_target) ** 2),
                                          catalogue_planet.name])
        catalogue_planets.sort()

        phot_filter = 'None detected'
        for key in read_local_log_profile('filter_key').split(','):
            check_filter = test_fits_keyword(fits, key)
            if check_filter[0]:
                if check_filter[2] in filter_map:
                    phot_filter = check_filter[2]
                    break
            if read_local_log_profile('filter') != '':
                phot_filter = read_local_log_profile('filter')

        if gauss_scatter < apperture_scatter:
            files_to_upload = ['PHOTOMETRY_GAUSS.txt', 'PHOTOMETRY_APERTURE.txt']
        else:
            files_to_upload = ['PHOTOMETRY_APERTURE.txt', 'PHOTOMETRY_GAUSS.txt']

        w = open('{0}/ExoClock_info.txt'.format(photometry_directory), 'w')
        w.write('\n'.join([
            'The ExoClock Project is an effort to keep the ephemerides of exoplanets as precise as \n'
            'possible for planning future observations. If you have observed an exoplnaet you can\n'
            'contribute your observation at: \n\nhttps://ariel-gbfu.azurewebsites.net\n\n'
            'File to upload: {0} \n(this is a suggestion based on the scatter \nof your light curves, '
            'you can also try \nuploading {1})'.format(*files_to_upload),
            '',
            'Planet: {0} \n(this is the closest known exoplanet found \nin the catalogue, if this is not the '
            'target \nyou observed, please ignore)'.format(str(catalogue_planets[0][1]).replace(' ', '')),
            '',
            'Time  format: JD_UTC \n(UTC-based Julian date)',
            '',
            'Flux format: Flux \n(flux of target over summed flux of comparisons)',
            '',
            'Filter: {0}'.format(phot_filter),
            '',
            'Exposure time in seconds: {0}'.format(exp_time),
        ]))
        w.close()

        finalise_window(root)

        while not exit_var_2.get():
            root.update()

        root.destroy()
        #
        #
        # #     all fov photometry
        #
        # root2 = Tk()
        #
        # exit_var = BooleanVar(value=False)
        #
        # def break_and_exit():
        #     exit_var.set(True)
        #
        # initialise_window(root2, exit_command=break_and_exit)
        #
        # label1 = Label(root2, text='EXTRACTING INDIVIDUAL LIGHT-CURVES')
        # label2 = Label(root2, text='POSITION:')
        # label3 = Label(root2, text=' ')
        # label4 = Label(root2, text='COMPLETE:')
        # label5 = Label(root2, text=' ', )
        # label6 = Label(root2, text='TIME LEFT:')
        # label7 = Label(root2, text=' ')
        #
        # setup_window(root2, [
        #     [[label1, 1, 2]],
        #     [[label2, 1], [label3, 2]],
        #     [[label4, 1], [label5, 2]],
        #     [[label6, 1], [label7, 2]],
        # ])
        #
        # all_targets = plc.open_dict(os.path.join(reduction_directory, 'all_frames.pickle'))
        #
        # percent = 0
        # lt0 = time.time()
        # for counter, target in enumerate(all_targets):
        #
        #     print(target)
        #
        #     gauss_flux = []
        #     gauss_sky = []
        #     aperture_flux = []
        #     aperture_sky = []
        #
        #     for frame in all_targets[target]:
        #
        #         fits = pf.open(os.path.join(reduction_directory, frame))
        #
        #         subframe = all_targets[target][frame]
        #
        #         if len(subframe.shape) == 0:
        #
        #             gauss_flux.append(0.0)
        #             gauss_sky.append(0.0)
        #             aperture_flux.append(0.0)
        #             aperture_sky.append(0.0)
        #
        #         else:
        #
        #             norm, floor, x_mean, y_mean, x_std, y_std = \
        #                 fit_2d_gauss_point(subframe, predicted_x_mean=int(len(subframe[0]) / 2),
        #                                    predicted_y_mean=int(len(subframe) / 2),
        #                                    search_window=search_window_std * star_std, stde=star_std, snr_lim=False)
        #
        #             gauss_flux.append(2 * np.pi * norm * x_std * y_std)
        #             gauss_sky.append(floor)
        #
        #             x_mean = int(len(subframe[0]) / 2)
        #             y_mean = int(len(subframe) / 2)
        #
        #             flux_area = subframe[int(y_mean) - 4 * star_std:
        #                                  int(y_mean) + 4 * star_std + 1,
        #                         int(x_mean) - 4 * star_std:
        #                         int(x_mean) + 4 * star_std + 1]
        #             flux_pixels = (2 * (4 * star_std) + 1) ** 2
        #             flux = np.sum(flux_area)
        #
        #             sky_area_1 = int(sky_inner_aperture * (4 * star_std))
        #             sky_area_2 = int(sky_outer_aperture * (4 * star_std))
        #             subframe[int(y_mean) - sky_area_1:int(y_mean) + sky_area_1 + 1,
        #             int(x_mean) - sky_area_1:int(x_mean) + sky_area_1 + 1] = 0
        #             sky_area = subframe[int(y_mean) - sky_area_2:int(y_mean) + sky_area_2 + 1,
        #                        int(x_mean) - sky_area_2:int(x_mean) + sky_area_2 + 1]
        #             sky_area = sky_area[np.where((sky_area > 0) &
        #                                          (sky_area < fits[1].header[mean_key] + 3 * fits[1].header[
        #                                              std_key]))]
        #             sky = np.sum(sky_area)
        #             sky_pixels = sky_area.size
        #
        #             aperture_flux.append(flux - flux_pixels * sky / sky_pixels)
        #             aperture_sky.append(sky / sky_pixels)
        #
        #     all_targets[target]['gauss_flux'] = gauss_flux
        #     all_targets[target]['gauss_sky'] = gauss_sky
        #     all_targets[target]['aperture_flux'] = aperture_flux
        #     all_targets[target]['aperture_sky'] = aperture_sky
        #
        #     # counter
        #
        #     new_percent = round(100 * (counter + 1) / float(len(all_targets)), 1)
        #     if new_percent != percent:
        #         lt1 = time.time()
        #         rm_time = (100 - new_percent) * (lt1 - lt0) / new_percent
        #         hours = rm_time / 3600.0
        #         minutes = (hours - int(hours)) * 60
        #         seconds = (minutes - int(minutes)) * 60
        #         label3.configure(text='     {0}     '.format(target))
        #         label5.configure(text='     {0}%    '.format(new_percent))
        #         label7.configure(text='     %dh %02dm %02ds     ' % (int(hours), int(minutes), int(seconds)))
        #         percent = new_percent
        #
        #     if counter == 0:
        #         finalise_window(root2)
        #
        #     root2.update()
        #
        #     if exit_var.get():
        #         break
        #
        #     if counter + 1 == len(science):
        #         write_log('pipeline', True, 'extraction_complete')
        #
        # root2.destroy()
        #
        # plc.save_dict(all_targets, os.path.join(photometry_directory, 'all_stars.pickle'))
