
import warnings
warnings.filterwarnings(
    'ignore', message='Matplotlib is building the font cache using fc-list. This may take a moment.')
warnings.filterwarnings(
    'ignore', message='The installed version of numexpr 2.4.4 is not supported in pandas and will be not be used')

import matplotlib
matplotlib.use('TkAgg')

import os
import numpy as np
import shutil
import hops.pylightcurve3 as plc
import matplotlib.cm as cm

from astropy.io import fits as pf
from scipy.optimize import curve_fit
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec

from matplotlib.backend_bases import FigureCanvasBase

from hops.hops_tools.logs import log
from hops.hops_tools.tests import *
from hops.hops_tools.windows import ProgressWindow


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
        fig.savefig(os.path.join(fitting_directory, 'corner.pdf'), transparent=False)

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
            logo_ax.imshow(log.holomon_logo_jpg)
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

            fig.text(0.04, fbottom + 2.5 * (1 - fbottom) / frow, 'relative flux (de-trended)', fontsize=fsbig, va='center',
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
            mtstr = '${0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11}{12}{13}{14}{15}$'.format(
                'T_\mathrm{BJD_{TDB}} = ',
                self.results['parameters']['mt']['print_value'],
                '_{-', self.results['parameters']['mt']['print_m_error'], '}',
                '^{+', self.results['parameters']['mt']['print_p_error'], '}',
                ' \quad \mathrm{O-C_{minutes}} = ',
                round((self.results['parameters']['mt']['value'] - prediction) * 24 * 60, 1),
                '_{-', round(self.results['parameters']['mt']['m_error']* 24 * 60, 1), '}',
                '^{+', round(self.results['parameters']['mt']['p_error']* 24 * 60, 1), '}'
                )

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

            fig.savefig(os.path.join(fitting_directory, 'detrended_model_300dpi.jpg'), dpi=300, transparent=False)
            fig.savefig(os.path.join(fitting_directory, 'detrended_model_900dpi.jpg'), dpi=900, transparent=False)
            fig.savefig(os.path.join(fitting_directory, 'detrended_model_1200dpi.jpg'), dpi=1200, transparent=False)

            return fig


def fitting():

    fitting_directory_base = log.read_local_log('pipeline', 'fitting_directory')
    reduction_directory = log.read_local_log('pipeline', 'reduction_directory')
    exposure_time_key = log.read_local_log('pipeline_keywords', 'exposure_time_key')
    light_curve_file = log.read_local_log('fitting', 'light_curve_file')
    scatter = log.read_local_log('fitting', 'scatter')
    iterations = log.read_local_log('fitting', 'iterations')
    burn = log.read_local_log('fitting', 'burn')
    planet = log.read_local_log('fitting', 'planet')
    metallicity = log.read_local_log('fitting', 'metallicity')
    temperature = log.read_local_log('fitting', 'temperature')
    logg = log.read_local_log('fitting', 'logg')
    phot_filter = log.read_local_log('fitting', 'phot_filter')
    period = log.read_local_log('fitting', 'period')
    mid_time = log.read_local_log('fitting', 'mid_time')
    mid_time_fit = log.read_local_log('fitting', 'mid_time_fit')
    rp_over_rs = log.read_local_log('fitting', 'rp_over_rs')
    rp_over_rs_fit = log.read_local_log('fitting', 'rp_over_rs_fit')
    sma_over_rs = log.read_local_log('fitting', 'sma_over_rs')
    sma_over_rs_fit = log.read_local_log('fitting', 'sma_over_rs_fit')
    inclination = log.read_local_log('fitting', 'inclination')
    inclination_fit = log.read_local_log('fitting', 'inclination_fit')
    eccentricity = log.read_local_log('fitting', 'eccentricity')
    periastron = log.read_local_log('fitting', 'periastron')
    observer = log.read_local_log('fitting', 'observer')
    observatory = log.read_local_log('fitting', 'observatory')
    telescope = log.read_local_log('fitting', 'telescope')
    camera = log.read_local_log('fitting', 'camera')
    target_ra_dec = log.read_local_log('fitting', 'target_ra_dec')

    def fit():

        if mid_time_fit:
            mid_time_fit_p = [mid_time + mid_time_fit[0], mid_time + mid_time_fit[1]]
        else:
            mid_time_fit_p = False
        if rp_over_rs_fit:
            rp_over_rs_fit_p = [rp_over_rs * rp_over_rs_fit[0], rp_over_rs * rp_over_rs_fit[1]]
        else:
            rp_over_rs_fit_p = False
        if sma_over_rs_fit:
            sma_over_rs_fit_p = [sma_over_rs * sma_over_rs_fit[0], sma_over_rs * sma_over_rs_fit[1]]
        else:
            sma_over_rs_fit_p = False
        if inclination_fit:
            inclination_fit_p = [inclination + inclination_fit[0], inclination + inclination_fit[1]]
        else:
            inclination_fit_p = False

        science = find_fits_files(os.path.join(reduction_directory, '*'))
        exp_time = pf.open(science[np.random.randint(len(science))])[1].header[exposure_time_key]

        light_curve = np.loadtxt(light_curve_file, unpack=True)

        date = plc.JD(light_curve[0][0]).utc.isoformat()[:15].replace('T', ' ')
        obs_duration = round(24 * (light_curve[0][-1] - light_curve[0][0]), 1)

        # filter out outliers

        light_curve_0 = light_curve[0]
        light_curve_1 = light_curve[1]

        light_curve_0 = light_curve_0[np.where(~np.isnan(light_curve_1))]
        light_curve_1 = light_curve_1[np.where(~np.isnan(light_curve_1))]

        moving_average = []
        for i in range(-10, 11):
            moving_average.append(np.roll(light_curve_1, i))

        median = np.median(moving_average, 0)
        med = np.median([np.abs(ff - median) for ff in moving_average], 0)

        flag = np.where((np.abs(light_curve_1 - median) < scatter * med))[0]

        light_curve_0 = light_curve_0[flag]
        light_curve_1 = light_curve_1[flag]

        # fix timing

        ra_dec_string = target_ra_dec.replace(':', ' ').split(' ')
        target = plc.Target(plc.Hours(*ra_dec_string[:3]), plc.Degrees(*ra_dec_string[3:]))
        light_curve_0 = np.array([plc.JD(ff + 0.5 * exp_time / 60.0 / 60.0 / 24.0).bjd_tdb(target) for ff in light_curve_0])

        # predictions

        limb_darkening_coefficients = plc.clablimb('claret', logg, max(4000, temperature), metallicity,
                                                   filter_map[phot_filter])

        predicted_mid_time = (mid_time + round((np.mean(light_curve_0) - mid_time) / period) * period)

        # define models

        def mcmc_f(time_array, detrend_zero, detrend_one, detrend_two, model_rp_over_rs, model_mid_time):

            data_delta_t = time_array - light_curve_0[0]

            detrend = detrend_zero * (1 + detrend_one * data_delta_t +
                                      detrend_two * data_delta_t * data_delta_t)
            transit_model = plc.transit_integrated('claret', limb_darkening_coefficients, model_rp_over_rs,
                                                   period, sma_over_rs, eccentricity,
                                                   inclination, periastron,
                                                   predicted_mid_time + model_mid_time,
                                                   time_array, float(exp_time), max(1, int(float(exp_time) / 10)))

            return detrend * transit_model

        def independent_f(time_array, detrend_zero, detrend_one, detrend_two, model_rp_over_rs, model_mid_time):

            data_delta_t = time_array - light_curve_0[0]

            detrend = detrend_zero * (1 + detrend_one * data_delta_t +
                                      detrend_two * data_delta_t * data_delta_t)
            transit_model = plc.transit_integrated('claret', limb_darkening_coefficients, model_rp_over_rs, period,
                                                   sma_over_rs, eccentricity, inclination,
                                                   periastron, predicted_mid_time + model_mid_time,
                                                   time_array, float(exp_time), max(1, int(float(exp_time) / 10)))

            return detrend, transit_model

        # set noise level

        if len(light_curve) == 3:
            sigma = light_curve[2][flag]
        else:
            sigma = np.array([np.roll(light_curve_1, ff) for ff in range(-10, 10)])
            sigma = np.std(sigma, 0)

        popt, pcov = curve_fit(mcmc_f, light_curve_0, light_curve_1,
                               p0=[np.mean(light_curve_1), 1, -1, rp_over_rs, 0],
                               sigma=sigma, maxfev=10000)

        fit_detrend, fit_transit_model = independent_f(light_curve_0, *popt)

        test = []
        for i in range(-int(len(light_curve_0) / 2), int(len(light_curve_0) / 2)):
            test.append([np.sum((light_curve_1 / fit_detrend - np.roll(fit_transit_model, i)) ** 2), i])
        test.sort()

        popt, pcov = curve_fit(mcmc_f, light_curve_0, light_curve_1, p0=[
            popt[0], popt[1], popt[2], popt[3],
            popt[4] + (test[0][1]) * exp_time / 60.0 / 60.0 / 24.0],
                               sigma=sigma, maxfev=10000)

        residuals = light_curve_1 - mcmc_f(light_curve_0, *popt)

        if len(light_curve) == 3:
            sigma *= np.std(residuals) / np.median(sigma)
        else:
            sigma = np.array([np.roll(residuals, ff) for ff in range(-10, 10)])
            sigma = np.std(sigma, 0)

        def function_to_call(counter):

            if counter.update_now:

                progress_bar_1['value'] = counter.percent
                percent_label_1.configure(
                    text='{0} % ({1}h {2}m {3}s left)'.format(counter.percent,
                                                              int(counter.time_left.split(':')[0]),
                                                              int(counter.time_left.split(':')[1]),
                                                              int(counter.time_left.split(':')[2])
                                                              ))
                show_progress.update()

            if show_progress.exit:
                return False
            else:
                return True

        mcmc_fit = HOPSTransitAndPolyFitting([[light_curve_0, light_curve_1, sigma]],
                                             method='claret',
                                             limb_darkening_coefficients=limb_darkening_coefficients,
                                             rp_over_rs=rp_over_rs,
                                             period=period,
                                             sma_over_rs=sma_over_rs,
                                             eccentricity=eccentricity,
                                             inclination=inclination,
                                             periastron=periastron,
                                             mid_time=mid_time,
                                             fit_rp_over_rs=rp_over_rs_fit_p,
                                             iterations=iterations,
                                             walkers=50,
                                             burn=burn,
                                             fit_first_order=True,
                                             fit_second_order=True,
                                             fit_period=False,
                                             fit_sma_over_rs=sma_over_rs_fit_p,
                                             fit_eccentricity=False,
                                             fit_inclination=inclination_fit_p,
                                             fit_periastron=False,
                                             fit_mid_time=mid_time_fit_p,
                                             precision=3,
                                             exp_time=round(exp_time, 1),
                                             time_factor=int(round(exp_time, 1) / 10),
                                             function_to_call=function_to_call
                                             )
        try:
            mcmc_fit.run_mcmc()
        except ValueError:
            show_progress.exit = True

        if not show_progress.exit:

            fitting_directory = fitting_directory_base

            if not os.path.isdir(fitting_directory):
                os.mkdir(fitting_directory)
            else:
                fi = 2
                while os.path.isdir('{0}_{1}'.format(fitting_directory, str(fi))):
                    fi += 1
                fitting_directory = '{0}_{1}'.format(fitting_directory, str(fi))
                os.mkdir(fitting_directory)

            mcmc_fit.save_results(os.path.join(fitting_directory, 'results.txt'))
            mcmc_fit.save_models(os.path.join(fitting_directory, 'model.txt'))
            mcmc_fit.save_detrended_models(os.path.join(fitting_directory, 'detrended_model.txt'))
            mcmc_fit.plot_hops_corner(fitting_directory)
            figure = mcmc_fit.plot_hops_output(
                planet,
                ['{0} (UT)\nDur: {1}h / Exp: {2}s\nFilter: {3}'.format(date, obs_duration, exp_time, phot_filter)],
                observer, '{0} / {1} / {2}'.format(observatory, telescope, camera), fitting_directory)
            shutil.copy('log.yaml', '{0}{1}log.yaml'.format(fitting_directory, os.sep))

            shutil.copy(log.fitting_output_description, fitting_directory)

            return figure

    def plot_fit(figure):

        root = ProgressWindow('HOPS - Fitting')

        canvas = root.FigureCanvasTkAgg(figure)
        canvas.get_tk_widget().pack()
        root.NavigationToolbar2Tk(canvas)

        root.show()

        while not root.exit:
            root.update()

        root.close()

    def run():

        figure = fit()
        if not show_progress.exit:
            plot_fit(figure)
            if not show_progress.exit:
                log.write_local_log('pipeline', True, 'fitting_complete')

        show_progress.close()

    show_progress = ProgressWindow('HOPS - Fitting', 0, 0, 5)

    label_1 = show_progress.Label(text="Running MCMC fitting")

    progress_bar_1 = show_progress.Progressbar()
    percent_label_1 = show_progress.Label(text='0.0 % (0h 0m 0s left)')

    show_progress.setup_window([
        [],
        [],
        [[label_1, 0]],
        [[progress_bar_1, 0]],
        [[percent_label_1, 0]],
        [],
        []
    ], main_font='Courier')

    show_progress.after(200, run)
    show_progress.loop()

