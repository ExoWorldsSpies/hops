from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .hops_basics import *


class HOPSTransitAndPolyFitting(plc.TransitAndPolyFitting):

    def plot_hops_corner(self, fitting_directory):

        def correlation(x, y):
            n = len(x)
            mx = np.mean(x)
            sx = np.std(x)
            my = np.mean(y)
            sy = np.std(y)
            return np.round(np.sum((x - mx) * (y - my)) / ((n - 1) * sx * sy), 2)

        def td_distribution(datax, datay, axx):

            datax = np.array(datax)
            median = np.median(datax)
            med = np.sqrt(np.median((datax - median) ** 2))
            xstep = med / 5.0
            xmin = min(datax)
            xmax = max(datax)
            x_size = int(round((xmax - xmin) / xstep)) + 1
            datax = np.int_((datax - xmin) / xstep)
            datay = np.array(datay)
            median = np.median(datay)
            med = np.sqrt(np.median((datay - median) ** 2))
            ystep = med / 5.0
            ymin = min(datay)
            ymax = max(datay)
            y_size = int(round((ymax - ymin) / ystep)) + 1
            datay = np.int_((datay - ymin) / ystep)

            yx_size = x_size * y_size
            yx = datay * x_size + datax

            yx = np.bincount(yx)
            yx = np.insert(yx, len(yx), np.zeros(yx_size - len(yx)))

            xx, yy = np.meshgrid(xmin + np.arange(x_size) * xstep, ymin + np.arange(y_size) * ystep)

            final = np.reshape(yx, (y_size, x_size))
            axx.imshow(np.where(final > 0, np.log(np.where(final > 0, final, 1)), 0),
                       extent=(np.min(xx), np.max(xx), np.min(yy), np.max(yy)),
                       cmap=cm.Greys, origin='lower', aspect='auto')

        if not self.mcmc_run_complete:
            raise RuntimeError('MCMC not completed')

        names = []
        results = []
        print_results = []
        errors1 = []
        print_errors1 = []
        errors2 = []
        print_errors2 = []
        errors = []
        traces = []
        traces_bins = []
        traces_counts = []

        for i in self.names:
            if self.results['parameters'][i]['initial']:
                names.append(self.results['parameters'][i]['print_name'])
                results.append(self.results['parameters'][i]['value'])
                print_results.append(self.results['parameters'][i]['print_value'])
                errors1.append(self.results['parameters'][i]['m_error'])
                print_errors1.append(self.results['parameters'][i]['print_m_error'])
                errors2.append(self.results['parameters'][i]['p_error'])
                print_errors2.append(self.results['parameters'][i]['print_p_error'])
                errors.append(0.5 * (self.results['parameters'][i]['m_error'] +
                                     self.results['parameters'][i]['p_error']))
                traces.append(self.results['parameters'][i]['trace'])
                traces_bins.append(self.results['parameters'][i]['trace_bins'])
                traces_counts.append(self.results['parameters'][i]['trace_counts'])

        all_var = len(traces)
        fig = Figure(figsize=(2.5 * all_var, 2.5 * all_var), tight_layout=False)
        canvas = FigureCanvasBase(fig)
        cmap = cm.get_cmap('brg')

        for var in range(len(names)):

            try:
                ax = fig.add_subplot(all_var, all_var, all_var * var + var + 1, facecolor='w')
            except AttributeError:
                ax = fig.add_subplot(all_var, all_var, all_var * var + var + 1, axisbg='w')

            ax.step(traces_bins[var], traces_counts[var], color='k', where='mid')

            ax.axvline(results[var], c='k')
            ax.axvline(results[var] - errors1[var], c='k', ls='--', lw=0.5)
            ax.axvline(results[var] + errors2[var], c='k', ls='--', lw=0.5)

            ax.set_xticks([results[var]])
            ax.set_yticks([0])
            ax.tick_params(left=False, right=False, top=False, bottom=False, labelbottom=False, labelleft=False)

            ax.set_xlabel('{0}\n{1}\n{2}\n{3}'.format(r'${0}$'.format(names[var]), r'${0}$'.format(print_results[var]),
                                                      r'$-{0}$'.format(print_errors1[var]),
                                                      r'$+{0}$'.format(print_errors2[var])), fontsize=20)

            ax.set_xlim(results[var] - 6 * errors[var], results[var] + 6 * errors[var])
            ax.set_ylim(0, ax.get_ylim()[1])

            for j in range(var + 1, all_var):

                try:
                    ax2 = fig.add_subplot(all_var, all_var, all_var * var + 1 + j, facecolor='w')
                except AttributeError:
                    ax2 = fig.add_subplot(all_var, all_var, all_var * var + 1 + j, axisbg='w')

                td_distribution(traces[j], traces[var], ax2)

                ax2.set_yticks([0])
                ax2.set_xticks([results[j]])
                ax2.tick_params(bottom=False, left=False, right=False, top=False, labelbottom=False,
                                labelleft=False, labelright=False, labeltop=False)

                ax2.set_xlim(results[j] - 6 * errors[j], results[j] + 6 * errors[j])
                ax2.set_ylim(results[var] - 6 * errors[var], results[var] + 6 * errors[var])
                text_x = ax2.get_xlim()[1] - 0.05 * (ax2.get_xlim()[1] - ax2.get_xlim()[0])
                text_y = ax2.get_ylim()[1] - 0.05 * (ax2.get_ylim()[1] - ax2.get_ylim()[0])
                ax2.text(text_x, text_y, '{0}{1}{2}'.format(r'$', str(correlation(traces[j], traces[var])), '$'),
                         color=cmap(abs(correlation(traces[j], traces[var])) / 2.),
                         fontsize=20, ha='right', va='top')

        fig.subplots_adjust(hspace=0, wspace=0)
        fig.savefig(os.path.join(fitting_directory, 'corner.pdf'), bbox_inches='tight', transparent=False)

    def plot_hops_output(self, target, data_dates, observer, observatory, fitting_directory):

        if target is None:
            target = ' '

        if data_dates is None:
            data_dates = map(str, ['set_{0}'.format(str(ff)) for ff in range(1, self.total_sets + 1)])

        for set_number in range(self.total_sets):

            funit = 1.0
            fcol = 7
            frow = 5
            fbottom = 0.11
            fright = 0.05
            fsmain = 10
            fsbig = 15
            fig = Figure(figsize=(funit * fcol / (1 - fright), funit * frow / (1 - fbottom)))
            canvas = FigureCanvasBase(fig)
            try:
                gs = gridspec.GridSpec(frow, fcol, fig, 0, fbottom, 1.0 - fright, 1.0, 0.0, 0.0)
            except TypeError:
                gs = gridspec.GridSpec(frow, fcol, 0, fbottom, 1.0 - fright, 1.0, 0.0, 0.0)

            fig.text(0.5, 0.94,  '{0}{1}{2}'.format('$\mathbf{', target, '}$'), fontsize=24, va='center', ha='center')
            fig.text(0.97, 0.97,  data_dates[set_number], fontsize=fsmain, va='top', ha='right')

            logo_ax = fig.add_subplot(gs[0, 0])
            logo_ax.imshow(holomon_logo_jpg)
            logo_ax.spines['top'].set_visible(False)
            logo_ax.spines['bottom'].set_visible(False)
            logo_ax.spines['left'].set_visible(False)
            logo_ax.spines['right'].set_visible(False)
            logo_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

            self.results = {ff: self.results[ff] for ff in self.results}

            period = self.results['parameters']['P']['value']
            mt = self.results['parameters']['mt']['value']
            mt += round((np.mean(self.data[set_number][0]) - mt) / period) * period

            prediction = (self.mid_time +
                          round((np.mean(self.data[set_number][0]) - self.mid_time) / self.period) * self.period)

            duration = plc.transit_duration(self.rp_over_rs, self.period, self.sma_over_rs,
                                            self.inclination, self.eccentricity, self.periastron)

            ingress = prediction - duration / 2
            egress = prediction + duration / 2

            set_indices = np.where(self.data_set_number == set_number)

            ax1 = fig.add_subplot(gs[1:4, 1:])

            ax1.plot(self.results['detrended_output_series']['phase'][set_indices],
                     self.results['detrended_input_series']['value'][set_indices], 'ko', ms=2)
            ax1.plot(self.results['detrended_output_series']['phase'][set_indices],
                     self.results['detrended_output_series']['model'][set_indices], 'r-')

            fig.text(0.04, fbottom + 2.5 * (1 - fbottom) / frow, 'relative flux', fontsize=fsbig, va='center',
                     ha='center', rotation='vertical')

            data_ymin = (min(self.results['detrended_input_series']['value'][set_indices])
                         - 3 * np.std(self.results['detrended_output_series']['residuals'][set_indices]))

            data_ymax = (max(self.results['detrended_input_series']['value'][set_indices])
                         + 2 * np.std(self.results['detrended_output_series']['residuals'][set_indices]))

            ax1.set_yticks(ax1.get_yticks()[np.where(ax1.get_yticks() > data_ymin)])

            ymin, ymax = data_ymax - 1.05 * (data_ymax - data_ymin), data_ymax

            ax1.set_ylim(ymin, ymax)

            x_max = max(np.abs(self.results['detrended_output_series']['phase'][set_indices]) +
                        0.05 * (max(self.results['detrended_output_series']['phase'][set_indices]) -
                                min(self.results['detrended_output_series']['phase'][set_indices])))

            ax1.set_xlim(-x_max, x_max)
            ax1.tick_params(labelbottom=False, labelsize=fsmain)

            rpstr = '{0}{1}{2}{3}{4}{5}{6}{7}'.format(
                r'$R_\mathrm{p}/R_* = ', self.results['parameters']['rp']['print_value'], '_{-',
                self.results['parameters']['rp']['print_m_error'], '}', '^{+',
                self.results['parameters']['rp']['print_p_error'], '}$')
            mtstr = '{0}{1}{2}{3}{4}{5}{6}{7}'.format(
                r'$T_\mathrm{BJD_{TDB}} = ', self.results['parameters']['mt']['print_value'], '_{-',
                self.results['parameters']['mt']['print_m_error'], '}', '^{+',
                self.results['parameters']['mt']['print_p_error'], '}$')

            ax1.text(0, ymin + 0.1 * (ymax - ymin),
                     '{0}{1}{2}'.format(rpstr, '\n', mtstr), ha='center', va='center', fontsize=fsmain)

            ax1.axvline((ingress - mt) / period, 0.3, 1.0, ls='--', c='k', lw=0.75)
            ax1.text((ingress - mt) / period, ax1.get_ylim()[0] + 0.3 * (ymax - ymin),
                     'predicted\ningress\nstart', ha='right', va='top', fontsize=fsmain)
            ax1.axvline((egress - mt) / period, 0.3, 1.0, ls='--', c='k', lw=0.75)
            ax1.text((egress - mt) / period, ax1.get_ylim()[0] + 0.3 * (ymax - ymin),
                     'predicted\negress\nend', ha='left', va='top', fontsize=fsmain)

            fig.text((1 - fright) / fcol, 1 - (1 - fbottom) / frow, '\n\n{0}\n{1}'.format(
                observer, observatory), fontsize=fsmain, ha='left', va='bottom')

            ax2 = fig.add_subplot(gs[4, 1:])
            ax2.plot(self.results['detrended_output_series']['phase'][set_indices],
                     self.results['detrended_output_series']['residuals'][set_indices], 'ko', ms=2)
            ax2.plot(self.results['detrended_output_series']['phase'][set_indices],
                     np.zeros_like(self.results['detrended_output_series']['phase'][set_indices]), 'r-')

            ax2.set_ylim(- 5 * np.std(self.results['detrended_output_series']['residuals'][set_indices]),
                         5 * np.std(self.results['detrended_output_series']['residuals'][set_indices]))

            ax2.set_xlabel('phase', fontsize=fsbig)
            fig.text(0.04, fbottom + 0.5 * (1 - fbottom) / frow, 'residuals', fontsize=fsbig, va='center', ha='center',
                     rotation='vertical')

            ax2.set_xlim(-x_max, x_max)
            ax2.tick_params(labelsize=fsmain)

            ax2.text(ax2.get_xlim()[0] + 0.02 * (ax2.get_xlim()[-1] - ax2.get_xlim()[0]),
                     ax2.get_ylim()[0] + 0.07 * (ax2.get_ylim()[-1] - ax2.get_ylim()[0]),
                     r'$\mathrm{rms}_\mathrm{res} = %.1e$' %
                     np.std(self.results['detrended_output_series']['residuals'][set_indices]),
                     fontsize=fsmain)

            fig.savefig(os.path.join(fitting_directory, 'detrended_model_300dpi.jpg'), dpi=300, transparent=True)
            fig.savefig(os.path.join(fitting_directory, 'detrended_model_900dpi.jpg'), dpi=900, transparent=True)
            fig.savefig(os.path.join(fitting_directory, 'detrended_model_1200dpi.jpg'), dpi=1200, transparent=True)

            return fig


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


def fitting():

    fitting_directory = read_local_log('pipeline', 'fitting_directory')
    reduction_directory = read_local_log('pipeline', 'reduction_directory')
    observation_date_key = read_local_log('pipeline_keywords', 'observation_date_key')
    observation_time_key = read_local_log('pipeline_keywords', 'observation_time_key')
    exposure_time_key = read_local_log('pipeline_keywords', 'exposure_time_key')
    light_curve_file = read_local_log('fitting', 'light_curve_file')
    binning = read_local_log('fitting', 'binning')
    scatter = read_local_log('fitting', 'scatter')
    iterations = read_local_log('fitting', 'iterations')
    burn = read_local_log('fitting', 'burn')
    planet = read_local_log('fitting', 'planet')
    planet_search = read_local_log('fitting', 'planet_search')
    metallicity = read_local_log('fitting', 'metallicity')
    temperature = read_local_log('fitting', 'temperature')
    logg = read_local_log('fitting', 'logg')
    phot_filter = read_local_log('fitting', 'phot_filter')
    period = read_local_log('fitting', 'period')
    period_fit = read_local_log('fitting', 'period_fit')
    mid_time = read_local_log('fitting', 'mid_time')
    mid_time_fit = read_local_log('fitting', 'mid_time_fit')
    rp_over_rs = read_local_log('fitting', 'rp_over_rs')
    rp_over_rs_fit = read_local_log('fitting', 'rp_over_rs_fit')
    sma_over_rs = read_local_log('fitting', 'sma_over_rs')
    sma_over_rs_fit = read_local_log('fitting', 'sma_over_rs_fit')
    inclination = read_local_log('fitting', 'inclination')
    inclination_fit = read_local_log('fitting', 'inclination_fit')
    eccentricity = read_local_log('fitting', 'eccentricity')
    eccentricity_fit = read_local_log('fitting', 'eccentricity_fit')
    periastron = read_local_log('fitting', 'periastron')
    periastron_fit = read_local_log('fitting', 'periastron_fit')
    observer = read_local_log('fitting', 'observer')
    observatory = read_local_log('fitting', 'observatory')
    telescope = read_local_log('fitting', 'telescope')
    camera = read_local_log('fitting', 'camera')
    target_ra_dec = read_local_log('fitting', 'target_ra_dec')

    limb_darkening_coefficients = plc.clablimb('claret', logg, max(4000, temperature), metallicity,
                                               filter_map[phot_filter])

    light_curve = np.loadtxt(light_curve_file, unpack=True)

    if binning > 1:
        start = len(light_curve[0]) - (len(light_curve[0]) // binning) * binning
        light_curve_0 = np.mean(np.reshape(light_curve[0][start:],
                                           (light_curve[0].size // binning, binning)), 1)
        light_curve_1 = np.mean(np.reshape(light_curve[1][start:],
                                           (light_curve[1].size // binning, binning)), 1)
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

    flag = np.where((np.abs(light_curve_1 - median) < scatter * med))[0]

    light_curve_0 = light_curve_0[flag]
    light_curve_1 = light_curve_1[flag]

    ra_target, dec_target = ra_dec_string_to_deg(target_ra_dec)
    light_curve_0 = np.array([plc.jd_utc_to_bjd_tdb(ra_target, dec_target, ff) for ff in light_curve_0])

    if not os.path.isdir(fitting_directory):
        os.mkdir(fitting_directory)
    else:
        fi = 2
        while os.path.isdir('{0}_{1}'.format(fitting_directory, str(fi))):
            fi += 1
        fitting_directory = '{0}_{1}'.format(fitting_directory, str(fi))
        os.mkdir(fitting_directory)

    if period_fit:
        period_fit = [period + period_fit[0], period + period_fit[1]]
    else:
        period_fit = False
    if mid_time_fit:
        mid_time_fit = [mid_time + mid_time_fit[0], mid_time + mid_time_fit[1]]
    else:
        mid_time_fit = False
    if rp_over_rs_fit:
        rp_over_rs_fit = [rp_over_rs * rp_over_rs_fit[0], rp_over_rs * rp_over_rs_fit[1]]
    else:
        rp_over_rs_fit = False
    if sma_over_rs_fit:
        sma_over_rs_fit = [sma_over_rs * sma_over_rs_fit[0], sma_over_rs * sma_over_rs_fit[1]]
    else:
        sma_over_rs_fit = False
    if inclination_fit:
        inclination_fit = [inclination + inclination_fit[0], inclination + inclination_fit[1]]
    else:
        inclination_fit = False
    if eccentricity_fit:
        eccentricity_fit = [eccentricity + eccentricity_fit[0], eccentricity + eccentricity_fit[1]]
    else:
        eccentricity_fit = False
    if periastron_fit:
        periastron_fit = [periastron + periastron_fit[0], periastron + periastron_fit[1]]
    else:
        periastron_fit = False

    science = find_fits_files(os.path.join(reduction_directory, '*'))
    fits = pf.open(science[0])
    if observation_date_key == observation_time_key:
        date = fits[1].header[observation_date_key].split('T')[0]
        local_time_1 = fits[1].header[observation_date_key].split('T')[1]
    else:
        date = fits[1].header[observation_date_key]
        local_time_1 = fits[1].header[observation_time_key]

    local_time_1 = ':'.join(local_time_1.split(':')[:2])
    obs_duration = round((light_curve_0[-1] - light_curve_0[0]) * 24, 1)

    fits = pf.open(science[np.random.randint(len(science))])
    exp_time = round(fits[1].header[exposure_time_key], 1)

    mcmc_fit = HOPSTransitAndPolyFitting([[light_curve_0, light_curve_1, np.ones_like(light_curve_1) *
                                         np.std(0.5 * (light_curve_1[:-1] - light_curve_1[1:]))]],
                                         method='claret',
                                         limb_darkening_coefficients=limb_darkening_coefficients,
                                         rp_over_rs=rp_over_rs,
                                         period=period,
                                         sma_over_rs=sma_over_rs,
                                         eccentricity=eccentricity,
                                         inclination=inclination,
                                         periastron=periastron,
                                         mid_time=mid_time,
                                         fit_rp_over_rs=rp_over_rs_fit,
                                         iterations=iterations,
                                         walkers=50,
                                         burn=burn,
                                         fit_first_order=True,
                                         fit_second_order=True,
                                         fit_period=period_fit,
                                         fit_sma_over_rs=sma_over_rs_fit,
                                         fit_eccentricity=eccentricity_fit,
                                         fit_inclination=inclination_fit,
                                         fit_periastron=periastron_fit,
                                         fit_mid_time=mid_time_fit,
                                         precision=3,
                                         exp_time=exp_time,
                                         time_factor=int(exp_time / 10),
                                         counter=False,
                                         counter_window='FITTING'
                                         )

    mcmc_fit.run_mcmc()
    mcmc_fit.save_results(os.path.join(fitting_directory, 'results.txt'))
    mcmc_fit.save_models(os.path.join(fitting_directory, 'model.txt'))
    mcmc_fit.save_detrended_models(os.path.join(fitting_directory, 'detrended_model.txt'))
    mcmc_fit.plot_hops_corner(fitting_directory)
    figure = mcmc_fit.plot_hops_output(
        planet_search,
        ['{0} {1} (UT)\nDur: {2}h / Exp: {3}s\nFiler: {4}'.format(date, local_time_1, obs_duration, exp_time,
                                                                  phot_filter)],
        observer, '{0} / {1} / {2}'.format(observatory, telescope, camera), fitting_directory)
    shutil.copy('log.yaml', '{0}{1}log.yaml'.format(fitting_directory, os.sep))

    roott = Tk()
    exit_var_2 = BooleanVar(value=False)

    def break_and_exit():
        exit_var_2.set(True)

    initialise_window(roott, exit_command=break_and_exit)

    canvas = FigureCanvasTkAgg(figure, roott)
    canvas.get_tk_widget().pack()
    NavigationToolbar2TkAgg(canvas, roott)

    finalise_window(roott, topmost=True)

    while not exit_var_2.get():
        roott.update()

    roott.destroy()
