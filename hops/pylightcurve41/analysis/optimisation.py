
__all__ = ['Fitting', 'EmceeFitting', 'values_to_print']

import emcee
import sys
import numpy as np
import warnings
from scipy.optimize import minimize
from scipy.stats import shapiro

from ..errors import *
from ..processes.counter import Counter
from ..processes.files import save_dict
from ..plots.plots_fitting import plot_mcmc_corner, plot_mcmc_traces, plot_mcmc_fitting
from ..analysis.distributions import one_d_distribution
from ..analysis.gaussian import gaussian
from .curve_fit import curve_fit
from .stats import *


class Fitting:

    def __init__(self, input_data_x, input_data_y, input_data_y_unc,
                 model, initials, limits1, limits2,
                 data_x_name='x', data_y_name='y', data_x_print_name='x', data_y_print_name='y',
                 walkers=None, iterations=None, burn_in=None, strech_prior=0.2,
                 parameters_names=None, parameters_print_names=None,
                 counter=None, window_counter=None,
                 optimise_initial_parameters=True,
                 optimise_initial_parameters_trials=4,
                 scale_uncertainties=False,
                 filter_outliers=False,
                 optimiser='emcee',
                 ):

        if optimiser not in ['emcee', 'curve_fit']:
            raise PyLCInputError('Optimiser {0} in not valid. Please choose between '
                                 'emcee, scipy_minimize.'.format(optimiser))

        self.input_data_x = np.array(input_data_x)
        self.input_data_y = np.array(input_data_y)
        self.input_data_y_unc = np.array(input_data_y_unc)
        self.input_data_y_unc_backup = np.array(input_data_y_unc)

        self.data_x_name = data_x_name
        self.data_y_name = data_y_name
        self.data_x_print_name = data_x_print_name
        self.data_y_print_name = data_y_print_name

        self.model = model
        self.initials = np.array(initials)
        self.limits1 = np.array(limits1)
        self.limits2 = np.array(limits2)
        self.fitted_parameters_indices = np.where(~np.isnan(self.limits1 * self.limits2))[0]
        self.dimensions = len(self.fitted_parameters_indices)
        self.internal_limits1 = self.limits1[self.fitted_parameters_indices]
        self.internal_limits2 = self.limits2[self.fitted_parameters_indices]
        self.internal_initials = (self.initials[self.fitted_parameters_indices] - self.internal_limits1) / (self.internal_limits2 - self.internal_limits1)

        self.names = ['p{0}'.format(ff) for ff in range(len(initials))]
        if parameters_names:
            self.names = parameters_names

        self.print_names = ['p{0}'.format(ff) for ff in range(len(initials))]
        if parameters_print_names:
            self.print_names = parameters_print_names

        self.optimiser = optimiser
        if self.optimiser == 'emcee':

            if not walkers:
                self.walkers = 3 * self.dimensions
                print('setting walkers = ', self.walkers)
            else:
                self.walkers = int(walkers)

            if not iterations:
                self.iterations = 5000
                print('setting iterations = ', self.iterations)
            else:
                self.iterations = int(iterations)

            if not burn_in:
                self.burn_in = int(self.iterations * 0.2)
                print('setting burn_in = ', self.burn_in)
            else:
                self.burn_in = int(burn_in)
                if self.burn_in >= self.iterations:
                    raise PyLCInputError('burn_in must be lower than iterations.')

            self.strech_prior = strech_prior
            self.walkers_initial_positions = None

            self.counter = counter
            if self.counter is True:
                self.counter = 'MCMC'
            self.window_counter = window_counter

            self.sampler = emcee.EnsembleSampler(self.walkers, self.dimensions, self.probability)
            self.progress = 0

        else:
            self.counter = None
            self.window_counter = None
            self.walkers = None
            self.iterations = None
            self.burn_in = None
            self.strech_prior = strech_prior
            self.sampler = None
            self.progress = None

        self.scale_uncertainties = scale_uncertainties
        self.filter_outliers = filter_outliers
        self.optimise_initial_parameters = optimise_initial_parameters
        self.optimise_initial_parameters_trials = optimise_initial_parameters_trials

        self.results = {
            'settings': {},
            'prefit': {},
            'input_series': {},
            'parameters': {},
            'parameters_final': [],
            'output_series': {},
            'statistics': {}}

        self.fitted_parameters = []

        self.mcmc_run_complete = False
        self.prefit_complete = False

        self.results['settings'][self.data_x_name] = np.ones_like(self.input_data_x) * self.input_data_x
        self.results['settings'][self.data_y_name] = np.ones_like(self.input_data_y) * self.input_data_y
        self.results['settings'][self.data_y_name + '_unc'] = np.ones_like(self.input_data_y_unc) * self.input_data_y_unc
        self.results['settings']['initials'] = self.initials
        self.results['settings']['limits1'] = self.limits1
        self.results['settings']['limits2'] = self.limits2
        self.results['settings']['walkers'] = self.walkers
        self.results['settings']['iterations'] = self.iterations
        self.results['settings']['burn_in'] = self.burn_in
        self.results['settings']['data_x_name'] = self.data_x_name
        self.results['settings']['data_x_print_name'] = self.data_x_print_name
        self.results['settings']['data_y_name'] = self.data_y_name
        self.results['settings']['data_y_print_name'] = self.data_y_print_name
        self.results['settings']['parameters_names'] = self.names
        self.results['settings']['parameters_print_names'] = self.print_names
        self.results['settings']['strech_prior'] = self.strech_prior
        self.results['settings']['optimise_initial_parameters'] = self.optimise_initial_parameters
        self.results['settings']['scale_uncertainties'] = self.scale_uncertainties
        self.results['settings']['filter_outliers'] = self.filter_outliers
        self.results['settings']['optimiser'] = self.optimiser

    def _pass(self):
        pass

    def internal_model(self, theta):
        parameters = self.initials
        parameters[self.fitted_parameters_indices] = theta * (self.internal_limits2 - self.internal_limits1) + self.internal_limits1
        return self.model(self.input_data_x, *parameters)

    def internal_model_curve_fit(self, x, *theta):
        theta = np.array(theta)
        parameters = self.initials
        parameters[self.fitted_parameters_indices] = theta * (self.internal_limits2 - self.internal_limits1) + self.internal_limits1
        return self.model(self.input_data_x, *parameters)

    def probability(self, theta):
        if np.prod((0 < theta) * (theta < 1)):
            chi = (self.input_data_y - self.internal_model(theta)) / self.input_data_y_unc
            return -0.5 * (np.sum(chi * chi) +
                           np.sum(np.log(2.0 * np.pi * (self.input_data_y_unc * self.input_data_y_unc))))
        else:
            return -np.inf

    def prefit(self):

        if self.scale_uncertainties or self.filter_outliers or self.optimise_initial_parameters:

            nll = lambda *args: -self.probability(*args)
            soln_test = nll(self.internal_initials)

            test_initials = self.internal_initials
            optimisation_ok = False

            for ii in range(self.optimise_initial_parameters_trials):
                print('Optimising initial parameters, attempt {0}: maximizing likelihood...'.format(ii+1))
                soln_i = minimize(nll, test_initials, method='Nelder-Mead')
                soln_test_i = nll(soln_i.x)

                if soln_i.success and soln_test_i < soln_test:

                    print('Optimisation completed.')
                    self.internal_initials = soln_i.x

                    if self.filter_outliers:
                        norm_res = (self.input_data_y - self.internal_model(self.internal_initials)) /self.input_data_y_unc
                        outliers = len(np.where(np.abs(norm_res) >= 3 * np.std(norm_res))[0])

                        while outliers > 0:

                            print('Filtering outliers...'.format(ii+1))

                            flags = np.where(np.abs(norm_res) >= 3 * np.std(norm_res))[0]
                            self.input_data_y_unc[flags] = 1000000000

                            soln_i = minimize(nll, self.internal_initials, method='Nelder-Mead')
                            self.internal_initials = soln_i.x
                            norm_res = (self.input_data_y - self.internal_model(self.internal_initials)) / self.input_data_y_unc
                            # print(np.std(norm_res), np.median(np.abs(norm_res - np.median(norm_res))))

                            outliers = len(np.where(np.abs(norm_res) >= 3 * np.std(norm_res))[0])
                            # print(outliers)

                    optimisation_ok = True
                    break

                elif soln_test_i < soln_test:
                    test_initials = soln_i.x
                    soln_test = nll(test_initials)
                else:
                    test_initials = self.internal_initials + np.random.normal(0, self.strech_prior/2.0, len(self.internal_initials))
                    test_initials = np.maximum(0, test_initials)
                    test_initials = np.minimum(1, test_initials)

            if not optimisation_ok:
                raise PyLCProcessError('Optimisation failed. You can try re-running, or increasing the strech_prior parameter, or the prior limits.')

        self.results['prefit']['outliers_map'] = self.input_data_y_unc == 1000000000
        outliers = np.where(self.results['prefit']['outliers_map'])
        points_to_use = np.where(~self.results['prefit']['outliers_map'])
        self.results['prefit']['outliers'] = len(outliers[0])
        self.input_data_x = self.input_data_x[points_to_use]
        self.input_data_y = self.input_data_y[points_to_use]
        self.input_data_y_unc = self.input_data_y_unc[points_to_use]

        if self.scale_uncertainties:
            scale_factor = np.sqrt(np.nanmean(((self.input_data_y - self.internal_model(self.internal_initials))**2) / (self.input_data_y_unc**2)))
            # import matplotlib.pyplot as plt
            # plt.figure()
            # plt.plot(self.input_data_x, self.input_data_y - self.internal_model(self.internal_initials), 'ko')
            # plt.plot(self.input_data_x, self.input_data_y, 'ko')
            # plt.plot(self.input_data_x, self.internal_model(self.internal_initials), 'r-')
            # plt.show()
        else:
            scale_factor = 1

        self.input_data_y_unc *= scale_factor
        self.results['prefit']['scale_factor'] = scale_factor

        print('Data-points excluded:', self.results['prefit']['outliers'])
        print('Scaling uncertainties by:', scale_factor)
        print('Initial parameters:')

        self.results['prefit']['initials'] = []
        for var in range(len(self.names)):
            if np.isnan(self.limits1[var]):
                self.results['prefit']['initials'].append(self.initials[var])
            else:
                print(
                    self.names[var], ': ',
                    self.internal_initials[np.where(self.fitted_parameters_indices == var)[0][0]] * (self.limits2[var] - self.limits1[var]) + self.limits1[var]
                )
                self.results['prefit']['initials'].append(self.internal_initials[np.where(self.fitted_parameters_indices == var)[0][0]] * (self.limits2[var] - self.limits1[var]) + self.limits1[var])

        self.results['input_series'][self.data_x_name] = self.input_data_x
        self.results['input_series'][self.data_y_name] = self.input_data_y
        self.results['input_series'][self.data_y_name + '_unc'] = self.input_data_y_unc

        self.prefit_complete = True

    def run(self):

        if not self.prefit_complete:
            self.prefit()

        # run sampler

        if self.optimiser == 'curve_fit':

            popt, pcov = curve_fit(self.internal_model_curve_fit, self.input_data_x,
                                   self.input_data_y, sigma=self.input_data_y_unc, p0=self.internal_initials)

            for var in range(len(self.names)):

                if np.isnan(self.limits1[var]):

                    variable = {'name': self.names[var], 'print_name': self.print_names[var],
                                'initial': None, 'min_allowed': None, 'max_allowed': None,
                                'trace': None, 'trace_bins': None, 'trace_counts': None,
                                'value': self.initials[var], 'm_error': None, 'p_error': None,
                                'print_value': self.initials[var], 'print_m_error': '-', 'print_p_error': '-'}

                else:

                    idx = np.where(self.fitted_parameters_indices == var)[0][0]

                    value = popt[idx]
                    min_value = value - np.sqrt(pcov[idx][idx])
                    max_value = value + np.sqrt(pcov[idx][idx])
                    value = value * (self.limits2[var] - self.limits1[var]) + self.limits1[var]
                    min_value = min_value * (self.limits2[var] - self.limits1[var]) + self.limits1[var]
                    max_value = max_value * (self.limits2[var] - self.limits1[var]) + self.limits1[var]

                    m_error = value - min_value
                    p_error = max_value - value

                    print_value, print_m_error, print_p_error = values_to_print(value, m_error, p_error)

                    variable = {'name': self.names[var], 'print_name': self.print_names[var],
                                'initial': self.initials[var],
                                'min_allowed': self.limits1[var], 'max_allowed': self.limits2[var],
                                'trace': None, 'trace_bins': None, 'trace_counts': None,
                                'value': value, 'm_error': m_error, 'p_error': p_error,
                                'print_value': print_value, 'print_m_error': print_m_error, 'print_p_error': print_p_error}

                    self.fitted_parameters.append(self.names[var])

                self.results['parameters'][self.names[var]] = variable
                self.results['parameters_final'].append(variable['value'])

            self.results['statistics']['corr_matrix'] = pcov
            self.results['statistics']['corr_variables'] = ','.join(self.fitted_parameters)

            self.postfit()

        elif self.optimiser == 'emcee':

            sys.setrecursionlimit(self.iterations)

            if self.counter:
                self.counter = Counter(self.counter, self.iterations, show_every=10, increment=10)
            else:
                self.counter = Counter('MCMC', self.iterations, show_every=self.iterations+1000, increment=10)

            if self.window_counter:
                from tkinter import Tk, Label
                self.root, self.mainloop_on, self.exit, self.name, self.jobs, self.widgets = Tk(), False, False, 'MCMC progress', [], []
                self.root.wm_title('MCMC progress'),  self.root.protocol('WM_DELETE_WINDOW', self.close)
                label1, label2, self.label = Label(self.root, textvar='   ...   '), Label(self.root, textvar='   ...   '), Label(self.root, textvar='')
                label1.grid(row=0, column=0), label2.grid(row=0, column=2),self.label.grid(row=0, column=1),self.root.update_idletasks()
                x, y = (self.root.winfo_screenwidth() - self.root.winfo_reqwidth()) / 2, (self.root.winfo_screenheight() - self.root.winfo_reqheight()) / 2
                self.root.geometry('+%d+%d' % (int(x), int(y))), self.root.update_idletasks(), self.root.wm_attributes("-topmost", 1), self.root.after_idle(self.root.attributes, '-topmost', 0), self.root.deiconify(), self.root.update_idletasks(),self.root.after(100, self.emcee_run)
                self.root.mainloop()

            else:
                self.emcee_run()

    def close(self):
        self.exit = True

    def emcee_run(self):

        if self.progress == 0:
            while self.progress == 0:
                try:
                    self.walkers_initial_positions = np.random.uniform(
                        (self.internal_initials - 0.5 * self.strech_prior)[:, None] * np.ones(self.walkers),
                        (self.internal_initials + 0.5 * self.strech_prior)[:, None] * np.ones(self.walkers))
                    self.walkers_initial_positions = np.swapaxes(self.walkers_initial_positions,0, 1)
                    self.walkers_initial_positions = np.minimum(self.walkers_initial_positions, 1)
                    self.walkers_initial_positions = np.maximum(self.walkers_initial_positions, 0)

                    self.sampler.run_mcmc(self.walkers_initial_positions, 10)
                    self.progress += 10
                    self.counter.update()

                    if self.window_counter:
                        if self.exit:
                            self.root.quit(), self.root.destroy()
                        else:
                            self.label['text'] = '\n{0}\n'.format(self.counter.text.replace(',', '\n'))
                            self.root.update_idletasks(), self.root.after(5, self.emcee_run)
                    else:
                        self.emcee_run()

                except ValueError:
                    pass

        elif self.progress < self.iterations:

            self.sampler.run_mcmc(None, 10)
            self.progress += 10
            self.counter.update()

            if self.window_counter:
                if self.exit:
                    self.root.quit(), self.root.destroy()
                else:
                    self.label['text'] = '\n{0}\n'.format(self.counter.text.replace(',', '\n'))
                    self.root.update_idletasks(), self.root.after(5, self.emcee_run)
            else:
                self.emcee_run()

        else:
            if self.window_counter:
                self.root.quit(), self.root.destroy()

            mcmc_results = self.sampler.get_chain()

            trace_to_analyse = 0
            vars_check = 0
            for var in range(len(self.names)):

                if not np.isnan(self.limits1[var]):
                    trace = mcmc_results[self.burn_in:, :, np.where(self.fitted_parameters_indices == var)[0][0]]
                    trace = trace.flatten()

                    median = np.median(trace)
                    mad = np.sqrt(np.median((trace - median) ** 2))

                    trace_to_analyse += (trace > (median - 10 * mad)) * (trace < (median + 10 * mad))

                    vars_check += 1

            trace_to_analyse = np.where(trace_to_analyse == vars_check)

            for var in range(len(self.names)):

                if np.isnan(self.limits1[var]):

                    variable = {'name': self.names[var], 'print_name': self.print_names[var],
                                'initial': None, 'min_allowed': None, 'max_allowed': None,
                                'trace': None, 'trace_bins': None, 'trace_counts': None,
                                'value': self.initials[var], 'm_error': None, 'p_error': None,
                                'print_value': str(self.initials[var]), 'print_m_error': '-', 'print_p_error': '-'}

                else:
                    trace = mcmc_results[self.burn_in:, :, np.where(self.fitted_parameters_indices == var)[0][0]]
                    trace = trace.flatten()

                    trace = trace[trace_to_analyse] * (self.limits2[var] - self.limits1[var]) + self.limits1[var]

                    bins, counts = one_d_distribution(trace)

                    min_value, value, max_value = np.quantile(trace, [0.16, 0.5, 0.84])
                    m_error = value - min_value
                    p_error = max_value - value

                    print_value, print_m_error, print_p_error = values_to_print(value, m_error, p_error)

                    variable = {'name': self.names[var], 'print_name': self.print_names[var],
                                'initial': self.initials[var],
                                'min_allowed': self.limits1[var], 'max_allowed': self.limits2[var],
                                'trace': trace, 'trace_bins': bins, 'trace_counts': counts,
                                'value': value, 'm_error': m_error, 'p_error': p_error,
                                'print_value': print_value, 'print_m_error': print_m_error, 'print_p_error': print_p_error}

                    self.fitted_parameters.append(self.names[var])

                self.results['parameters'][self.names[var]] = variable
                self.results['parameters_final'].append(variable['value'])

            to_correlate = []
            for parameter in self.fitted_parameters:
                to_correlate.append(self.results['parameters'][parameter]['trace'])
            correlation_matrix = np.corrcoef(to_correlate)
            self.results['statistics']['corr_matrix'] = correlation_matrix
            self.results['statistics']['corr_variables'] = ','.join(self.fitted_parameters)

            self.postfit()

    def postfit(self):

        self.results['output_series']['model'] = self.model(self.input_data_x, *self.results['parameters_final'])
        self.results['output_series']['residuals'] = self.input_data_y - self.results['output_series']['model']

        norm_residuals = self.results['output_series']['residuals'] / self.input_data_y_unc
        try:
            norm_residuals = np.swapaxes([self.input_data_x, norm_residuals], 0, 1)
            norm_residuals = sorted(norm_residuals, key=lambda x: x[0])
            norm_residuals = np.swapaxes(norm_residuals, 0, 1)[1]
        except ValueError:
            pass

        res_autocorr = np.correlate(norm_residuals, norm_residuals, mode='full')
        res_autocorr = res_autocorr[res_autocorr.size // 2:]
        res_autocorr /= res_autocorr[0]

        limit3_autocorr = gaussian(np.log10(len(norm_residuals)), 1.08401, 0.03524, -0.26884, 1.49379)

        res_shapiro = shapiro(norm_residuals)

        limit3_shapiro = gaussian(np.log10(len(norm_residuals)), 0.65521, 0.00213, -0.21983, 0.96882)

        self.results['statistics']['res_autocorr'] = res_autocorr
        self.results['statistics']['res_max_autocorr'] = np.max(np.abs(res_autocorr[1:]))
        self.results['statistics']['res_max_autocorr_flag'] = np.max(np.abs(res_autocorr[1:])) > limit3_autocorr
        self.results['statistics']['res_shapiro'] = res_shapiro[0]
        self.results['statistics']['res_shapiro_flag'] = (1 - res_shapiro[0]) > limit3_shapiro
        self.results['statistics']['res_mean'] = np.mean(self.results['output_series']['residuals'])
        self.results['statistics']['res_std'] = np.std(self.results['output_series']['residuals'])
        self.results['statistics']['res_rms'] = np.sqrt(np.mean(self.results['output_series']['residuals']**2))
        self.results['statistics']['res_chi_sqr'] = np.sum(norm_residuals ** 2)
        self.results['statistics']['res_red_chi_sqr'] = (
                self.results['statistics']['res_chi_sqr'] / (len(self.input_data_y_unc) - len(self.fitted_parameters)))

        self.mcmc_run_complete = True

    def save_all(self, export_file):

        if not self.mcmc_run_complete:
            raise PyLCProcessError('MCMC not completed')

        save_dict(self.results, export_file)

    def save_results(self, export_file):

        if not self.mcmc_run_complete:
            raise PyLCProcessError('MCMC not completed')

        cols = [
            ['# variable'],
            ['fix/fit'],
            ['value'],
            ['uncertainty'],
            ['initial'],
            ['min. allowed'],
            ['max. allowed']
        ]

        for i in self.names:

            cols[0].append(self.results['parameters'][i]['name'])
            if self.results['parameters'][i]['initial'] is None:
                cols[1].append('fix')
                cols[2].append(self.results['parameters'][i]['print_value'])
                cols[3].append(' ')
                cols[4].append(' ')
                cols[5].append(' ')
                cols[6].append(' ')
            else:
                cols[1].append('fit')
                cols[2].append(self.results['parameters'][i]['print_value'])
                cols[3].append(
                    '-{0} +{1}'.format(
                        self.results['parameters'][i]['print_m_error'],
                        self.results['parameters'][i]['print_p_error'])
                )
                cols[4].append(str(self.results['parameters'][i]['initial']))
                cols[5].append(str(self.results['parameters'][i]['min_allowed']))
                cols[6].append(str(self.results['parameters'][i]['max_allowed']))

        for col in cols:
            col_length = np.max([len(ff) for ff in col])
            for ff in range(len(col)):
                col[ff] = col[ff] + ' ' * (col_length - len(col[ff]))

        lines = []

        for row in range(len(cols[0])):
            lines.append('  '.join([col[row] for col in cols]))

        lines.append('')
        lines.append('#Pre-fit:')
        lines.append('#Number of outliers removed: {0}'.format(self.results['prefit']['outliers']))
        lines.append('#Uncertainties scale factor: {0}'.format(self.results['prefit']['scale_factor']))
        lines.append('')
        lines.append('#Residuals:')
        lines.append('#Mean: {0}'.format(self.results['statistics']['res_mean']))
        lines.append('#STD: {0}'.format(self.results['statistics']['res_std']))
        lines.append('#RMS: {0}'.format(self.results['statistics']['res_rms']))
        lines.append('#Chi squared: {0}'.format(self.results['statistics']['res_chi_sqr']))
        lines.append('#Reduced chi squared: {0}'.format(self.results['statistics']['res_red_chi_sqr']))
        lines.append('#Max auto-correlation: {0}'.format(self.results['statistics']['res_max_autocorr']))
        lines.append('#Max auto-correlation flag: {0}'.format(self.results['statistics']['res_max_autocorr_flag']))
        lines.append('#Shapiro test: {0}'.format(self.results['statistics']['res_shapiro']))
        lines.append('#Shapiro test flag: {0}'.format(self.results['statistics']['res_shapiro_flag']))

        w = open(export_file, 'w')
        w.write('\n'.join(lines))
        w.close()

    def plot_fitting(self, export_file):

        plot_mcmc_fitting(self, export_file)

    def plot_corner(self, export_file):

        plot_mcmc_corner(self, export_file)

    def plot_traces(self, export_file):

        plot_mcmc_traces(self, export_file)


# decimal points and rounding


def values_to_print(value, error_minus, error_plus):

    value = float(value)
    error_minus = float(error_minus)
    error_plus = float(error_plus)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        if error_minus >= 1.0 or error_minus == 0.0:
            digit1 = 2
        else:
            str_error_minus = '{0:.{test}f}'.format(error_minus, test=10 + abs(int(np.log10(error_minus))))
            digit1 = np.where([ff not in ['0', '.'] for ff in str_error_minus])[0][0]

        if error_plus >= 1.0 or error_plus == 0.0:
            digit2 = 2
        else:
            str_error_plus = '{0:.{test}f}'.format(error_plus, test=10 + abs(int(np.log10(error_plus))))
            digit2 = np.where([ff not in ['0', '.'] for ff in str_error_plus])[0][0]

    width = max(2, digit1, digit2)

    print_value = '{0:.{width}f}'.format(round(value, width), width=width)
    print_m_error = '{0:.{width}f}'.format(round(error_minus, width), width=width)
    print_p_error = '{0:.{width}f}'.format(round(error_plus, width), width=width)

    return print_value, print_m_error, print_p_error


# for compatibility reasons
class EmceeFitting(Fitting):

    def __init__(self, *args, **kwargs):

        try:
            kwargs['walkers'] = args[7]
            kwargs['iterations'] = args[8]
            kwargs['burn_in'] = args[9]
            if args[8] > 10000:
                print('In version 4.1 of PyLightcurve, the total number of calculations is given by iterations x walkers.')
                kwargs['iterations'] /= kwargs['walkers']
                kwargs['burn_in'] /= kwargs['walkers']
        except:
            kwargs['walkers'] = None
            kwargs['iterations'] = None
            kwargs['burn_in'] = None

        print('Warning plc.EmceeFitting will not be supported in future versions.')
        print('Use plc.Fitting instead, and set optimiser="emcee".')

        Fitting.__init__(self, *args[:7], **kwargs)

        self.run_mcmc = self.run
