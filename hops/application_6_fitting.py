
import os
import glob
import numpy as np
import shutil
from scipy.optimize import curve_fit as scipy_curve_fit
import hops.pylightcurve3 as plc
import matplotlib.image as mpimg
import warnings
from hops.application_windows import MainWindow, AddOnWindow


def curve_fit(*args, **kwargs):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message='Covariance of the parameters could not be estimated')
        return scipy_curve_fit(*args, **kwargs)


filter_map = {'Clear': 'V', 'Luminance': 'V',
              'U': 'U', 'B': 'B', 'V': 'V', 'R': 'R', 'I': 'I', 'H': 'H', 'J': 'J', 'K': 'K',
              'u': 'u', 'b': 'b', 'v': 'v', 'y': 'y',
              'u\'': 'u,', 'g\'': 'g,', 'r\'': 'r,', 'i\'': 'i,', 'z\'': 'z,',
              'Astrodon ExoPlanet-BB': 'R',
              'UV': 'U', 'Rc': 'R', 'Ic': 'I', 'Re': 'R', 'Ie': 'I', 'Y': 'y,', 'r': 'r,', 'z': 'z,', 'i': 'i,',
              }


class FittingWindow(MainWindow):

    def __init__(self, log):

        MainWindow.__init__(self, log, name='HOPS - Fitting', position=2)

        # data

        self.all_frames = plc.open_dict(self.log.all_frames)
        self.science_files = []
        self.exposure_time = []
        for science_file in self.all_frames:
            if not self.all_frames[science_file][self.log.skip_key]:
                self.science_files.append(science_file)
                self.exposure_time.append(self.all_frames[science_file][self.log.get_param('exposure_time_key')])

        self.exposure_time = np.median(self.exposure_time)
        self.science_files.sort()
        self.fits = plc.open_fits(os.path.join(self.log.reduction_directory, self.science_files[0]))[1]

        # extra windows

        self.export_window = AddOnWindow(self, name='Export Results', position=1)
        self.export_window_database = self.export_window.DropDown(initial='ExoClock', options=['ExoClock', 'ETD'],
                                                           command=self.export_window_update)
        self.export_window_camera_gain = self.export_window.Entry(value='1', instance=float)
        self.export_window_time_shift = 0
        self.export_window_message = self.export_window.Label(text=' ')

        self.export_window.setup_window([
            [],
            [[self.export_window.Label(text='Database:'), 0]],
            [[self.export_window_database, 0]],
            [],
            [[self.export_window.Label(text='Camera gain:'), 0]],
            [[self.export_window_camera_gain, 0]],
            [],
            [[self.export_window_message, 0]],
            [],
            [[self.export_window.Button(text='EXPORT', command=self.export), 0]],
            []
        ])

        self.export_window_update()

        # main window

        photometry_files = (glob.glob(os.path.join('{0}*'.format(self.log.photometry_directory_base), self.log.light_curve_aperture_file)) +
                            glob.glob(os.path.join('{0}*'.format(self.log.photometry_directory_base), self.log.light_curve_gauss_file)))
        photometry_files_sort = []
        for ff in photometry_files:
            folder = os.path.split(ff)[0]
            if '_' not in folder:
                num = 1
            else:
                num = int(folder.split('_')[1])
            
            photometry_files_sort.append([num, ff])
        
        photometry_files_sort.sort()
        
        photometry_files = [ff[1] for ff in photometry_files_sort]
        photometry_files = ['Choose Light-curve file'] + photometry_files

        if self.log.get_param('light_curve_file') not in photometry_files:
            self.log.set_param('light_curve_file', 'Choose Light-curve file')
        
        self.light_curve_file = self.DropDown(initial=self.log.get_param('light_curve_file'), options=photometry_files,
                                              width=40, command=self.update_window)

        self.scatter = self.Entry(value=self.log.get_param('scatter'), instance=float, command=self.update_window)
        self.iterations = self.Entry(value=self.log.get_param('iterations'), instance=int, command=self.update_window)
        self.burn = self.Entry(value=self.log.get_param('burn'), instance=int, command=self.update_window)
        self.manual_planet = self.CheckButton(initial=self.log.get_param('manual_planet'),
                                              text='Enter param. manually', command=self.update_window)

        self.planet = self.Entry(value=self.log.get_param('planet'), command=self.update_window_no_refit_if_auto)
        self.target_ra_dec = self.Entry(value=self.log.get_param('target_ra_dec'), command=self.update_window_no_refit_if_auto)
        self.metallicity = self.Entry(value=self.log.get_param('metallicity'), instance=float, command=self.update_window_no_refit_if_auto)
        self.temperature = self.Entry(value=self.log.get_param('temperature'), instance=float, command=self.update_window_no_refit_if_auto)
        self.logg = self.Entry(value=self.log.get_param('logg'), instance=float, command=self.update_window_no_refit_if_auto)
        self.period = self.Entry(value=self.log.get_param('period'), instance=float, command=self.update_window_no_refit_if_auto)
        self.mid_time = self.Entry(value=self.log.get_param('mid_time'), instance=float, command=self.update_window_no_refit_if_auto)
        self.rp_over_rs = self.Entry(value=self.log.get_param('rp_over_rs'), instance=float, command=self.update_window_no_refit_if_auto)
        self.sma_over_rs = self.Entry(value=self.log.get_param('sma_over_rs'), instance=float, command=self.update_window_no_refit_if_auto)
        self.inclination = self.Entry(value=self.log.get_param('inclination'), instance=float, command=self.update_window_no_refit_if_auto)
        self.eccentricity = self.Entry(value=self.log.get_param('eccentricity'), instance=float, command=self.update_window_no_refit_if_auto)
        self.periastron = self.Entry(value=self.log.get_param('periastron'), instance=float, command=self.update_window_no_refit_if_auto)

        ra_target, dec_target = self.log.get_param('target_ra_dec').split(' ')
        ecc_planet = plc.find_nearest(plc.Target(plc.Hours(ra_target), plc.Degrees(dec_target)))

        self.auto_planet = self.Label(text=ecc_planet.planet.name)
        self.auto_target_ra_dec = self.Label(text='{0} {1}'.format(ecc_planet.star.ra, ecc_planet.star.dec))
        self.auto_metallicity = self.Label(text=ecc_planet.star.met, instance=float)
        self.auto_temperature = self.Label(text=ecc_planet.star.teff, instance=float)
        self.auto_logg = self.Label(text=ecc_planet.star.logg, instance=float)
        self.auto_period = self.Label(text=ecc_planet.planet.period, instance=float)
        self.auto_mid_time = self.Label(text=0, instance=float)
        self.auto_rp_over_rs = self.Label(text=ecc_planet.planet.rp_over_rs, instance=float)
        self.auto_sma_over_rs = self.Label(text=ecc_planet.planet.sma_over_rs, instance=float)
        self.auto_eccentricity = self.Label(text=ecc_planet.planet.eccentricity, instance=float)
        self.auto_inclination = self.Label(text=ecc_planet.planet.inclination, instance=float)
        self.auto_periastron = self.Label(text=ecc_planet.planet.periastron, instance=float)
        target = plc.Target(plc.Hours(ecc_planet.star.ra), plc.Degrees(ecc_planet.star.dec))
        if ecc_planet.planet.time_format in ['BJD_TDB', 'BJD_TT']:
            self.auto_mid_time.set(ecc_planet.planet.mid_time)
        elif ecc_planet.planet.time_format == 'BJD_UTC':
            self.auto_mid_time.set(plc.BJDUTC(ecc_planet.planet.mid_time, target).bjd_tdb(target))
        elif ecc_planet.planet.time_format in ['HJD_TDB', 'HJD_TT']:
            self.auto_mid_time.set(plc.HJDTDB(ecc_planet.planet.mid_time, target).bjd_tdb(target))
        elif ecc_planet.planet.time_format == 'HJD_UTC':
            self.auto_mid_time.set(plc.HJDUTC(ecc_planet.planet.mid_time, target).bjd_tdb(target))
        elif ecc_planet.planet.time_format == 'JD_UTC':
            self.auto_mid_time.set(plc.JD(ecc_planet.planet.mid_time).bjd_tdb(target))

        if self.planet.get() == 'Choose Planet':
            self.planet.set(self.auto_planet.get())
            self.target_ra_dec.set(self.auto_target_ra_dec.get())
            self.metallicity.set(self.auto_metallicity.get())
            self.temperature.set(self.auto_temperature.get())
            self.logg.set(self.auto_logg.get())
            self.period.set(self.auto_period.get())
            self.mid_time.set(self.auto_mid_time.get())
            self.rp_over_rs.set(self.auto_rp_over_rs.get())
            self.sma_over_rs.set(self.auto_sma_over_rs.get())
            self.inclination.set(self.auto_inclination.get())
            self.eccentricity.set(self.auto_eccentricity.get())
            self.periastron.set(self.auto_periastron.get())

        # set progress variables, useful for updating the window

        self.refit = self.BooleanVar(True)

        self.fitting_button = self.Button(text='RUN FITTING', command=self.run_fitting,
                                          bg='green', highlightbackground='green')

        self.export_button = self.Button(text='EXPORT FOR DATABASES (ExoClock / ETD)', command=self.export_window.show)

        # define preview window

        funit = 0.8
        fcol = 7
        frow = 5
        fbottom = 0.11
        fright = 0.02
        self.fsmain = 10
        self.fsbig = 20
        self.preview_figure = self.FigureWindow(figsize=(0.4, 0.7, 5, 5, 0.9))
        try:
            gs = gridspec.GridSpec(frow, fcol, self.preview_figure.figure, 0.03, fbottom, 1.0 - fright, 1.0, 0.0, 0.0)
        except TypeError:
            gs = gridspec.GridSpec(frow, fcol, 0.03, fbottom, 1.0 - fright, 1.0, 0.0, 0.0)

        logo_ax = self.preview_figure.figure.add_subplot(gs[0, 0])
        logo_ax.imshow(mpimg.imread(self.log.files['holomon_logo']))
        logo_ax.spines['top'].set_visible(False)
        logo_ax.spines['bottom'].set_visible(False)
        logo_ax.spines['left'].set_visible(False)
        logo_ax.spines['right'].set_visible(False)
        logo_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

        self.ax1 = self.preview_figure.figure.add_subplot(gs[1:4, 1:])
        self.ax2 = self.preview_figure.figure.add_subplot(gs[4, 1:])
        self.preview_figure.figure.text(0.04, fbottom + 2.5 * (1 - fbottom) / frow, 'relative flux (de-trended)', fontsize=self.fsmain, va='center',
                 ha='center', rotation='vertical')
        self.preview_figure.figure.text(0.04, fbottom + 0.5 * (1 - fbottom) / frow, 'residuals', fontsize=self.fsmain, va='center',
                 ha='center', rotation='vertical')

        self.title1 = self.preview_figure.figure.text(0.5, 0.94, '', fontsize=self.fsbig, va='center', ha='center')
        self.title2 = self.preview_figure.figure.text(0.97, 0.97, '', fontsize=self.fsmain, va='top', ha='right')
        self.title3 = self.preview_figure.figure.text(0.03 + (1 - fright) / fcol, 1 - (1 - fbottom) / frow, '', fontsize=self.fsmain, ha='left', va='bottom')

        self.setup_window([
            [[self.preview_figure, 0, 1, 24]],
            [[self.Label(text='Light-curve file'), 1, 3]],
            [[self.light_curve_file, 1, 3]],
            [[self.export_button, 1, 3]],
            [],
            [[self.Label(text='Catalogue param.'), 2], [self.manual_planet, 3]],

            [[self.Label(text='Planet'), 1], [self.auto_planet, 2], [self.planet, 3]],

            [[self.Label(text='Planet RA DEC\n(hh:mm:ss +/-dd:mm:ss)'), 1], [self.auto_target_ra_dec, 2], [self.target_ra_dec, 3]],

            [[self.Label(text='Period [days]'), 1], [self.auto_period, 2], [self.period, 3]],

            [[self.Label(text='Mid-time [BJD_TDB]'), 1], [self.auto_mid_time, 2], [self.mid_time, 3]],

            [[self.Label(text='Rp/Rs'), 1], [self.auto_rp_over_rs, 2], [self.rp_over_rs, 3]],

            [[self.Label(text='a/Rs'), 1], [self.auto_sma_over_rs, 2], [self.sma_over_rs, 3]],

            [[self.Label(text='Inclination [deg]'), 1], [self.auto_inclination, 2], [self.inclination, 3]],

            [[self.Label(text='Eccentricity'), 1], [self.auto_eccentricity, 2], [self.eccentricity, 3]],

            [[self.Label(text='Periastron [deg]'), 1], [self.auto_periastron, 2], [self.periastron, 3]],

            [[self.Label(text='M* [Fe/H, dex]'), 1], [self.auto_metallicity, 2], [self.metallicity, 3]],

            [[self.Label(text='T* [K]'), 1], [self.auto_temperature, 2], [self.temperature, 3]],

            [[self.Label(text='log(g*) [cm/s^2]'), 1], [self.auto_logg, 2], [self.logg, 3]],

            [[self.Label(text='Scatter limit'), 1], [self.scatter, 2], [self.Label(text='default = 3.0'), 3]],
            [[self.Label(text='MCMC Iterations'), 1], [self.iterations, 2], [self.Label(text='default = 150000'), 3]],
            [[self.Label(text='MCMC Burn-in\n(less than Iterations)'), 1], [self.burn, 2], [self.Label(text='default = 100000'), 3]],

            [[self.Button(text='RETURN TO MAIN MENU', command=self.close), 1, 3]],
            [[self.Button(text='SAVE OPTIONS & RETURN TO MAIN MENU', command=self.save_and_return), 1, 3]],
            [[self.fitting_button, 1], [self.Button(text='RETURN TO PHOTOMETRY', command=self.return_to_photometry), 2, 2]],
            [],

        ])

        self.after(self.update_window)

    def update_window(self, *entry):

        if not entry:
            pass

        if not os.path.isfile(self.light_curve_file.get()):

            self.scatter.disable()
            self.iterations.disable()
            self.burn.disable()
            self.metallicity.disable()
            self.temperature.disable()
            self.logg.disable()
            self.period.disable()
            self.mid_time.disable()
            self.rp_over_rs.disable()
            self.sma_over_rs.disable()
            self.inclination.disable()
            self.eccentricity.disable()
            self.periastron.disable()
            self.target_ra_dec.disable()
            self.planet.disable()
            self.manual_planet.disable()
            self.fitting_button.disable()

        else:

            self.scatter.activate()
            self.iterations.activate()
            self.burn.activate()
            self.manual_planet.activate()
            self.planet.activate()
            self.target_ra_dec.activate()
            self.metallicity.activate()
            self.temperature.activate()
            self.logg.activate()
            self.period.activate()
            self.mid_time.activate()
            self.rp_over_rs.activate()
            self.sma_over_rs.activate()
            self.inclination.activate()
            self.eccentricity.activate()
            self.periastron.activate()
            self.target_ra_dec.activate()

            if self.manual_planet.get():
                self.planet.activate()
                self.target_ra_dec.activate()
                self.metallicity.activate()
                self.temperature.activate()
                self.logg.activate()
                self.period.activate()
                self.mid_time.activate()
                self.rp_over_rs.activate()
                self.sma_over_rs.activate()
                self.eccentricity.activate()
                self.inclination.activate()
                self.periastron.activate()
                planet_to_plot = self.planet.get()
                target_ra_dec_to_plot = self.target_ra_dec.get()
                metallicity_to_plot = self.metallicity.get()
                temperature_to_plot = self.temperature.get()
                logg_to_plot = self.logg.get()
                period_to_plot = self.period.get()
                mid_time_to_plot = self.mid_time.get()
                rp_over_rs_to_plot = self.rp_over_rs.get()
                sma_over_rs_to_plot = self.sma_over_rs.get()
                inclination_to_plot = self.inclination.get()
                eccentricity_to_plot = self.eccentricity.get()
                periastron_to_plot = self.periastron.get()
            else:
                self.planet.disable()
                self.target_ra_dec.disable()
                self.metallicity.disable()
                self.temperature.disable()
                self.logg.disable()
                self.period.disable()
                self.mid_time.disable()
                self.rp_over_rs.disable()
                self.sma_over_rs.disable()
                self.eccentricity.disable()
                self.inclination.disable()
                self.periastron.disable()
                planet_to_plot = self.auto_planet.get()
                target_ra_dec_to_plot = self.auto_target_ra_dec.get()
                metallicity_to_plot = self.auto_metallicity.get()
                temperature_to_plot = self.auto_temperature.get()
                logg_to_plot = self.auto_logg.get()
                period_to_plot = self.auto_period.get()
                mid_time_to_plot = self.auto_mid_time.get()
                rp_over_rs_to_plot = self.auto_rp_over_rs.get()
                sma_over_rs_to_plot = self.auto_sma_over_rs.get()
                inclination_to_plot = self.auto_inclination.get()
                eccentricity_to_plot = self.auto_eccentricity.get()
                periastron_to_plot = self.auto_periastron.get()

            light_curve = np.loadtxt(self.light_curve_file.get(), unpack=True)

            date = plc.JD(light_curve[0][0]).utc.isoformat()[:16].replace('T', ' ')
            self.observation_date = plc.JD(light_curve[0][0]).utc.isoformat()[:16].split('T')[0]
            obs_duration = round(24 * (light_curve[0][-1] - light_curve[0][0]), 1)

            self.title1.set_text('{0}{1}{2}'.format('$\mathbf{', planet_to_plot, '}$'))
            self.title2.set_text('{0} (UT)\nDur: {1}h / Exp: {2}s\nFilter: {3}'.format(date, obs_duration, self.exposure_time,
                                                                                       self.log.get_param('filter')))
            self.title3.set_text(
                '\n\n{0}\n{1}'.format(
                    self.log.get_param('observer'), '{0} / {1} / {2}'.format(
                        self.log.get_param('observatory'), self.log.get_param('telescope'), self.log.get_param('camera'))))

            if self.refit.get():

                self.ax1.cla()
                self.ax2.cla()
                try:

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

                    flag = np.where((np.abs(light_curve_1 - median) < self.scatter.get() * med))[0]
                    outliers = np.where((np.abs(light_curve_1 - median) >= self.scatter.get() * med))[0]

                    light_curve_0_outliers = light_curve_0[outliers]
                    light_curve_1_outliers = light_curve_1[outliers]
                    light_curve_0 = light_curve_0[flag]
                    light_curve_1 = light_curve_1[flag]

                    # fix timing

                    ra_dec_string = target_ra_dec_to_plot.replace(':', ' ').split(' ')
                    target = plc.Target(plc.Hours(*ra_dec_string[:3]), plc.Degrees(*ra_dec_string[3:]))
                    light_curve_0 = np.array([plc.JD(ff + 0.5 * self.exposure_time / 60.0 / 60.0 / 24.0).bjd_tdb(target) for ff in light_curve_0])

                    # predictions

                    limb_darkening_coefficients = plc.clablimb('claret', logg_to_plot, max(4000, temperature_to_plot),
                                                               metallicity_to_plot, filter_map[self.log.get_param('filter')])

                    predicted_mid_time = (mid_time_to_plot +
                                          round((np.mean(light_curve_0) - mid_time_to_plot) / period_to_plot) * period_to_plot)

                    predicted_transit_model = plc.transit_integrated('claret', limb_darkening_coefficients, rp_over_rs_to_plot,
                                                          period_to_plot, sma_over_rs_to_plot, eccentricity_to_plot,
                                                          inclination_to_plot, periastron_to_plot, predicted_mid_time,
                                                          light_curve_0, float(self.exposure_time), max(1, int(float(self.exposure_time) / 10)))

                    # define models

                    def mcmc_f(time_array, detrend_zero, detrend_one, detrend_two, model_rp_over_rs, model_mid_time):

                        data_delta_t = time_array - light_curve_0[0]

                        detrend = detrend_zero * (1 + detrend_one * data_delta_t +
                                                  detrend_two * data_delta_t * data_delta_t)
                        transit_model = plc.transit_integrated('claret', limb_darkening_coefficients, model_rp_over_rs,
                                                    period_to_plot, sma_over_rs_to_plot, eccentricity_to_plot,
                                                    inclination_to_plot, periastron_to_plot,
                                                    predicted_mid_time + model_mid_time,
                                                    time_array, float(self.exposure_time), max(1, int(float(self.exposure_time) / 10)))

                        return detrend * transit_model

                    def independent_f(time_array, detrend_zero, detrend_one, detrend_two, model_rp_over_rs, model_mid_time):

                        data_delta_t = time_array - light_curve_0[0]

                        detrend = detrend_zero * (1 + detrend_one * data_delta_t +
                                                  detrend_two * data_delta_t * data_delta_t)
                        transit_model = plc.transit_integrated('claret', limb_darkening_coefficients, model_rp_over_rs, period_to_plot,
                                                    sma_over_rs_to_plot, eccentricity_to_plot, inclination_to_plot,
                                                    periastron_to_plot, predicted_mid_time + model_mid_time,
                                                    time_array, float(self.exposure_time), max(1, int(float(self.exposure_time) / 10)))

                        return detrend, transit_model

                    # set noise level

                    if len(light_curve) == 3:
                        sigma = light_curve[2][flag]
                    else:
                        sigma = np.array([np.roll(light_curve_1, ff) for ff in range(-10, 10)])
                        sigma = np.std(sigma, 0)

                    popt, pcov = curve_fit(mcmc_f, light_curve_0, light_curve_1,
                                           p0=[np.mean(light_curve_1), 1, -1, rp_over_rs_to_plot, 0],
                                           sigma=sigma, maxfev=10000)

                    fit_detrend, fit_transit_model = independent_f(light_curve_0, *popt)

                    test = []
                    for i in range(-int(len(light_curve_0) / 2), int(len(light_curve_0) / 2)):
                        test.append([np.sum((light_curve_1 / fit_detrend - np.roll(fit_transit_model, i)) ** 2), i])
                    test.sort()

                    popt, pcov = curve_fit(mcmc_f, light_curve_0, light_curve_1,
                                           p0=[popt[0], popt[1], popt[2], popt[3],
                                               popt[4] + (test[0][1]) * self.exposure_time / 60.0 / 60.0 / 24.0],
                                           sigma=sigma, maxfev=10000)

                    residuals = light_curve_1 - mcmc_f(light_curve_0, *popt)

                    sigma = np.array([np.roll(residuals, ff) for ff in range(-10, 10)])
                    sigma = np.std(sigma, 0)

                    # results

                    popt, pcov = curve_fit(mcmc_f, light_curve_0, light_curve_1, p0=popt, sigma=sigma, maxfev=10000)
                    new_mid_time = predicted_mid_time + popt[-1]

                    fit_detrend, fit_transit_model = independent_f(light_curve_0, *popt)
                    phase = np.array((light_curve_0 - new_mid_time) / period_to_plot)
                    detrended_data = light_curve_1 / fit_detrend
                    detrended_uncertainties = sigma / fit_detrend
                    residuals = detrended_data - fit_transit_model
                    std_res = np.std(residuals)
                    res_autocorr = np.correlate(residuals, residuals, mode='full')
                    res_autocorr = round(np.max(np.abs(
                        res_autocorr[res_autocorr.size // 2:] / res_autocorr[res_autocorr.size // 2:][0])[1:]), 2)

                    fit_detrend_outliers, fit_transit_model_outliers = independent_f(light_curve_0_outliers, *popt)
                    phase_outliers = np.array((light_curve_0_outliers - new_mid_time) / period_to_plot)
                    detrended_data_outliers = light_curve_1_outliers / fit_detrend_outliers

                    # plot

                    self.ax1.errorbar(phase, detrended_data, detrended_uncertainties, c='k', fmt='o', ms=2, lw=0.5)
                    self.ax1.plot(phase, detrended_data, 'ko', ms=2, label='De-trended data')

                    data_ymin = min(detrended_data) - 5 * std_res
                    data_ymax = max(detrended_data) + 2 * std_res
                    self.ax1.set_yticks(self.ax1.get_yticks()[np.where(self.ax1.get_yticks() > data_ymin)])
                    ymin, ymax = data_ymax - 1.3 * (data_ymax - data_ymin), data_ymax
                    self.ax1.set_ylim(ymin, ymax)
                    x_max = max(np.abs(phase) + 0.05 * (max(phase) - min(phase)))

                    self.ax1.set_xlim(-x_max, x_max)
                    self.ax1.tick_params(labelbottom=False, labelsize=self.fsmain)

                    self.ax1.plot(phase_outliers, detrended_data_outliers, 'ro', ms=2, label='Filtered outliers')
                    self.ax1.plot(phase, fit_transit_model, 'r-', label='Best-fit model\n(T$_0$ = {0}, R$_p$/R$_s$ = {1})'.format(
                        round(new_mid_time, 5), round(popt[3], 5)))
                    self.ax1.plot(phase, predicted_transit_model, 'c-', label='Expected model\n(T$_0$ = {0}, R$_p$/R$_s$ = {1})'.format(
                        round(predicted_mid_time, 5), round(rp_over_rs_to_plot, 5)))
                    self.ax1.legend(loc=3, fontsize='small')

                    self.ax1.text(x_max - 0.4*x_max, ymin + 0.1 * (ymax - ymin),
                             'O-C:\n{0} min'.format(round((new_mid_time - predicted_mid_time) * 24 * 60, 1)))

                    self.ax2.errorbar(phase, residuals, detrended_uncertainties, c='k', fmt='o', ms=2, lw=0.5)
                    self.ax2.plot(phase, residuals, 'ko', ms=2)
                    self.ax2.plot(phase, np.zeros_like(phase), 'r-')

                    self.ax2.set_ylim(- 5 * std_res, 5 * std_res)

                    self.ax2.set_xlabel('phase', fontsize=self.fsmain)

                    self.ax2.set_xlim(-x_max, x_max)
                    self.ax2.tick_params(labelsize=self.fsmain)

                    self.ax2.text(self.ax2.get_xlim()[0] + 0.02 * (self.ax2.get_xlim()[-1] - self.ax2.get_xlim()[0]),
                             self.ax2.get_ylim()[1] - 0.15 * (self.ax2.get_ylim()[-1] - self.ax2.get_ylim()[0]),
                             r'STD = %.1f $‰$' % round((std_res * 1000), 1),
                             fontsize=self.fsmain)
                    self.ax2.text(self.ax2.get_xlim()[0] + 0.02 * (self.ax2.get_xlim()[-1] - self.ax2.get_xlim()[0]),
                             self.ax2.get_ylim()[0] + 0.07 * (self.ax2.get_ylim()[-1] - self.ax2.get_ylim()[0]),
                             r'AutoCorr = %.1f' %res_autocorr,
                             fontsize=self.fsmain)

                    self.mcmc_fit = HOPSTransitAndPolyFitting([[light_curve_0, light_curve_1, sigma]],
                                                              method='claret',
                                                              limb_darkening_coefficients=limb_darkening_coefficients,
                                                              rp_over_rs=rp_over_rs_to_plot,
                                                              period=period_to_plot,
                                                              sma_over_rs=sma_over_rs_to_plot,
                                                              eccentricity=eccentricity_to_plot,
                                                              inclination=inclination_to_plot,
                                                              periastron=periastron_to_plot,
                                                              mid_time=mid_time_to_plot,
                                                              fit_rp_over_rs=np.array(self.log.get_param('rp_over_rs_fit')) * rp_over_rs_to_plot,
                                                              iterations=self.iterations.get(),
                                                              walkers=50,
                                                              burn=self.burn.get(),
                                                              fit_first_order=True,
                                                              fit_second_order=True,
                                                              fit_period=False,
                                                              fit_sma_over_rs=self.log.get_param('sma_over_rs_fit'),
                                                              fit_eccentricity=False,
                                                              fit_inclination=self.log.get_param('inclination_fit'),
                                                              fit_periastron=False,
                                                              fit_mid_time=list(np.array(self.log.get_param('mid_time_fit')) + mid_time_to_plot),
                                                              precision=3,
                                                              exp_time=round(self.exposure_time, 1),
                                                              time_factor=int(round(self.exposure_time, 1) / 10),
                                                              counter=False,
                                                              )

                    self.fitting_button.activate()

                    self.preview_figure.draw()

                except Exception as e:
                    # print(e)
                    self.fitting_button.disable()

    def update_window_no_refit(self, *entry):
        self.refit.set(False)
        self.update_window()
        self.refit.set(True)

    def update_window_no_refit_if_auto(self, *entry):
        if self.manual_planet.get():
            self.refit.set(True)
        else:
            self.refit.set(False)
        self.update_window()
        self.refit.set(True)

    def save_settings(self):

        if self.manual_planet.get():
            planet_to_plot = self.planet.get()
            target_ra_dec_to_plot = self.target_ra_dec.get()
            metallicity_to_plot = self.metallicity.get()
            temperature_to_plot = self.temperature.get()
            logg_to_plot = self.logg.get()
            period_to_plot = self.period.get()
            mid_time_to_plot = self.mid_time.get()
            rp_over_rs_to_plot = self.rp_over_rs.get()
            sma_over_rs_to_plot = self.sma_over_rs.get()
            inclination_to_plot = self.inclination.get()
            eccentricity_to_plot = self.eccentricity.get()
            periastron_to_plot = self.periastron.get()
        else:
            planet_to_plot = self.auto_planet.get()
            target_ra_dec_to_plot = self.auto_target_ra_dec.get()
            metallicity_to_plot = self.auto_metallicity.get()
            temperature_to_plot = self.auto_temperature.get()
            logg_to_plot = self.auto_logg.get()
            period_to_plot = self.auto_period.get()
            mid_time_to_plot = self.auto_mid_time.get()
            rp_over_rs_to_plot = self.auto_rp_over_rs.get()
            sma_over_rs_to_plot = self.auto_sma_over_rs.get()
            inclination_to_plot = self.auto_inclination.get()
            eccentricity_to_plot = self.auto_eccentricity.get()
            periastron_to_plot = self.auto_periastron.get()

        self.log.set_param('light_curve_file', self.light_curve_file.get())
        self.log.set_param('scatter', self.scatter.get())
        self.log.set_param('iterations', self.iterations.get())
        self.log.set_param('burn', self.burn.get())
        self.log.set_param('manual_planet', self.manual_planet.get())
        self.log.set_param('planet', planet_to_plot)
        self.log.set_param('target_ra_dec', target_ra_dec_to_plot)
        self.log.set_param('metallicity', metallicity_to_plot)
        self.log.set_param('temperature', temperature_to_plot)
        self.log.set_param('logg', logg_to_plot)
        self.log.set_param('period', period_to_plot)
        self.log.set_param('mid_time', mid_time_to_plot)
        self.log.set_param('rp_over_rs', rp_over_rs_to_plot)
        self.log.set_param('sma_over_rs', sma_over_rs_to_plot)
        self.log.set_param('inclination', inclination_to_plot)
        self.log.set_param('eccentricity', eccentricity_to_plot)
        self.log.set_param('periastron', periastron_to_plot)

        # save log

        self.log.save_local_log()

    def save_and_return(self):

        self.save_settings()

        self.close()

    def run_fitting(self):

        self.save_settings()

        self.log.set_param('proceed', 'run_fitting')
        self.close()

    def return_to_photometry(self):

        self.save_settings()

        self.log.set_param('proceed', 'return_to_photometry')
        self.close()

    def export(self, *entry):

        light_curve = np.loadtxt(self.light_curve_file.get(), unpack=True)

        w = open('HOPS_{5}_{0}_{1}_{2}_{3}s_for_{4}.txt'.format(
            self.observation_date,
            self.log.get_param('planet'),
            self.log.get_param('filter'),
            self.exposure_time,
            self.export_window_database.get(),
            self.light_curve_file.get().replace(os.sep, '_').replace('.txt', '')
        ), 'w')

        for i in range(len(light_curve[0])):
            w.write('{0:.10f}\t{1:.10f}\t{2:.10f}\n'.format(light_curve[0][i] + self.export_window_time_shift,
                                             light_curve[1][i],
                                             light_curve[2][i] * (1 / np.sqrt(self.export_window_camera_gain.get())),
                                             ))

        return None

    def export_window_update(self, *entry):

        if self.export_window_database.get() == 'ExoClock':

            self.export_window_camera_gain.set(1.0)
            self.export_window_time_shift = 0

            self.export_window_message.set('No transformation is needed to upload your results to\nExoClock ' \
                                         '(www.exoclock.space).\n\nWhen uploading you will need to use the following ' \
                                         'details:\n' \
                                         'Planet: {0}\n'\
                                         'Time format: JD_UTC\n'\
                                         'Time stamp: Exposure start\n'\
                                         'Flux format: Flux\n'\
                                         'Filter: {1}\n'\
                                         'Exposure time in seconds: {2}'.format(self.log.get_param('planet'),
                                                                                self.log.get_param('filter'),
                                                                                self.exposure_time))

        else:
            self.export_window_time_shift = self.exposure_time / 60 / 60 / 24 / 2

            self.export_window_message.set('To upload your results to ETD (var2.astro.cz/ETD)\nthe following modifications ' \
                                         'are applied:\n\n' \
                                         '1. Half of the exposure time is added as the time-stamp required on ETD ' \
                                         'referes to the exposure mid-time\n' \
                                         '2. The uncertainties are  mupliplied by 1/SQRT(gain) to represent the ' \
                                         'uncertainies in electrons rather than counts.\n\n'\
                                         'When uploading you will need to use the following details:\n' \
                                         'Planet: {0}\n'\
                                         'JD format: geocentric\n'\
                                         'Brightness column: in flux'.format(self.log.get_param('planet')))


class FittingProgressWindow(MainWindow):

    def __init__(self, log):

        MainWindow.__init__(self, log, name='HOPS - Fitting progress', position=5)

        # load variabvles)

        self.all_frames = plc.open_dict(self.log.all_frames)
        self.science_files = []
        self.exposure_time = []
        for science_file in self.all_frames:
            if not self.all_frames[science_file][self.log.skip_key]:
                self.science_files.append(science_file)
                self.exposure_time.append(self.all_frames[science_file][self.log.get_param('exposure_time_key')])

        self.exposure_time = np.median(self.exposure_time)
        self.science_files.sort()
        self.fits = plc.open_fits(os.path.join(self.log.reduction_directory, self.science_files[0]))[1]

        self.light_curve_file = log.get_param('light_curve_file')
        self.scatter = self.log.get_param('scatter')
        self.iterations = self.log.get_param('iterations')
        self.burn = self.log.get_param('burn')
        self.planet = self.log.get_param('planet')
        self.target_ra_dec = self.log.get_param('target_ra_dec')
        self.metallicity = self.log.get_param('metallicity')
        self.temperature = self.log.get_param('temperature')
        self.logg = self.log.get_param('logg')
        self.period = self.log.get_param('period')
        self.mid_time = self.log.get_param('mid_time')
        self.rp_over_rs = self.log.get_param('rp_over_rs')
        self.sma_over_rs = self.log.get_param('sma_over_rs')
        self.inclination = self.log.get_param('inclination')
        self.eccentricity = self.log.get_param('eccentricity')
        self.periastron = self.log.get_param('periastron')

        light_curve = np.loadtxt(self.light_curve_file, unpack=True)

        date = plc.JD(light_curve[0][0]).utc.isoformat()[:16].replace('T', ' ')
        self.observation_date = plc.JD(light_curve[0][0]).utc.isoformat()[:16].split('T')[0]
        obs_duration = round(24 * (light_curve[0][-1] - light_curve[0][0]), 1)

        self.title1 = '{0}{1}{2}'.format('$\mathbf{', self.planet, '}$')
        self.title2 = '{0} (UT)\nDur: {1}h / Exp: {2}s\nFilter: {3}'.format(date, obs_duration, self.exposure_time,
                                                                  self.log.get_param('filter'))
        self.title3 = '\n\n{0}\n{1}'.format(
                self.log.get_param('observer'), '{0} / {1} / {2}'.format(
                    self.log.get_param('observatory'), self.log.get_param('telescope'), self.log.get_param('camera')))

        light_curve_0 = light_curve[0]
        light_curve_1 = light_curve[1]

        light_curve_0 = light_curve_0[np.where(~np.isnan(light_curve_1))]
        light_curve_1 = light_curve_1[np.where(~np.isnan(light_curve_1))]

        moving_average = []
        for i in range(-10, 11):
            moving_average.append(np.roll(light_curve_1, i))

        median = np.median(moving_average, 0)
        med = np.median([np.abs(ff - median) for ff in moving_average], 0)

        flag = np.where((np.abs(light_curve_1 - median) < self.scatter * med))[0]
        outliers = np.where((np.abs(light_curve_1 - median) >= self.scatter * med))[0]

        light_curve_0_outliers = light_curve_0[outliers]
        light_curve_1_outliers = light_curve_1[outliers]
        light_curve_0 = light_curve_0[flag]
        light_curve_1 = light_curve_1[flag]

        # fix timing

        ra_dec_string = self.target_ra_dec.replace(':', ' ').split(' ')
        target = plc.Target(plc.Hours(*ra_dec_string[:3]), plc.Degrees(*ra_dec_string[3:]))
        light_curve_0 = np.array(
            [plc.JD(ff + 0.5 * self.exposure_time / 60.0 / 60.0 / 24.0).bjd_tdb(target) for ff in light_curve_0])

        # predictions

        limb_darkening_coefficients = plc.clablimb('claret', self.logg, max(4000, self.temperature),
                                                   self.metallicity, filter_map[self.log.get_param('filter')])

        self.predicted_mid_time = (self.mid_time +
                              round((np.mean(light_curve_0) - self.mid_time) / self.period) * self.period)

        predicted_transit_model = plc.transit_integrated('claret', limb_darkening_coefficients, self.rp_over_rs,
                                                         self.period, self.sma_over_rs, self.eccentricity,
                                                         self.inclination, self.periastron, self.predicted_mid_time,
                                                         light_curve_0, float(self.exposure_time),
                                                         max(1, int(float(self.exposure_time) / 10)))

        # define models

        def mcmc_f(time_array, detrend_zero, detrend_one, detrend_two, model_rp_over_rs, model_mid_time):

            data_delta_t = time_array - light_curve_0[0]

            detrend = detrend_zero * (1 + detrend_one * data_delta_t +
                                      detrend_two * data_delta_t * data_delta_t)
            transit_model = plc.transit_integrated('claret', limb_darkening_coefficients, model_rp_over_rs,
                                                   self.period, self.sma_over_rs, self.eccentricity,
                                                   self.inclination, self.periastron,
                                                   self.predicted_mid_time + model_mid_time,
                                                   time_array, float(self.exposure_time),
                                                   max(1, int(float(self.exposure_time) / 10)))

            return detrend * transit_model

        def independent_f(time_array, detrend_zero, detrend_one, detrend_two, model_rp_over_rs, model_mid_time):

            data_delta_t = time_array - light_curve_0[0]

            detrend = detrend_zero * (1 + detrend_one * data_delta_t +
                                      detrend_two * data_delta_t * data_delta_t)
            transit_model = plc.transit_integrated('claret', limb_darkening_coefficients, model_rp_over_rs,
                                                   self.period,
                                                   self.sma_over_rs, self.eccentricity, self.inclination,
                                                   self.periastron, self.predicted_mid_time + model_mid_time,
                                                   time_array, float(self.exposure_time),
                                                   max(1, int(float(self.exposure_time) / 10)))

            return detrend, transit_model

        # set noise level

        if len(light_curve) == 3:
            sigma = light_curve[2][flag]
        else:
            sigma = np.array([np.roll(light_curve_1, ff) for ff in range(-10, 10)])
            sigma = np.std(sigma, 0)

        popt, pcov = curve_fit(mcmc_f, light_curve_0, light_curve_1,
                               p0=[np.mean(light_curve_1), 1, -1, self.rp_over_rs, 0],
                               sigma=sigma, maxfev=10000)

        fit_detrend, fit_transit_model = independent_f(light_curve_0, *popt)

        test = []
        for i in range(-int(len(light_curve_0) / 2), int(len(light_curve_0) / 2)):
            test.append([np.sum((light_curve_1 / fit_detrend - np.roll(fit_transit_model, i)) ** 2), i])
        test.sort()

        popt, pcov = curve_fit(mcmc_f, light_curve_0, light_curve_1,
                               p0=[popt[0], popt[1], popt[2], popt[3],
                                   popt[4] + (test[0][1]) * self.exposure_time / 60.0 / 60.0 / 24.0],
                               sigma=sigma, maxfev=10000)

        residuals = light_curve_1 - mcmc_f(light_curve_0, *popt)

        sigma = np.array([np.roll(residuals, ff) for ff in range(-10, 10)])
        sigma = np.std(sigma, 0)

        self.mcmc_fit = HOPSTransitAndPolyFitting([[light_curve_0, light_curve_1, sigma]],
                                                  method='claret',
                                                  limb_darkening_coefficients=limb_darkening_coefficients,
                                                  rp_over_rs=self.rp_over_rs,
                                                  period=self.period,
                                                  sma_over_rs=self.sma_over_rs,
                                                  eccentricity=self.eccentricity,
                                                  inclination=self.inclination,
                                                  periastron=self.periastron,
                                                  mid_time=self.mid_time,
                                                  fit_rp_over_rs=np.array(self.log.get_param('rp_over_rs_fit')) * self.rp_over_rs,
                                                  iterations=int(self.iterations / (self.iterations/1000)),
                                                  walkers=50,
                                                  burn=self.burn,
                                                  fit_first_order=True,
                                                  fit_second_order=True,
                                                  fit_period=False,
                                                  fit_sma_over_rs=self.log.get_param('sma_over_rs_fit'),
                                                  fit_eccentricity=False,
                                                  fit_inclination=self.log.get_param('inclination_fit'),
                                                  fit_periastron=False,
                                                  fit_mid_time=list(np.array(self.log.get_param('mid_time_fit')) + self.mid_time),
                                                  precision=3,
                                                  exp_time=round(self.exposure_time, 1),
                                                  time_factor=int(round(self.exposure_time, 1) / 10),
                                                  counter=False,
                                                  )

        self.progress_fitting = self.Progressbar(task="MCMC Fitting")

        self.stop_and_return_button = self.Button(
            text='STOP FITTING', command=self.trigger_exit)
        self.return_to_fitting_button = self.Button(text='RETURN TO FITTING MENU', command=self.close)
        self.proceed_button = self.Button(text='EXIT', command=self.proceed)

        # structure

        structure = []

        structure.append([])

        structure.append([[self.progress_fitting, 0, 4]])
        structure.append([])
        structure.append([[self.stop_and_return_button, 0],
                          [self.return_to_fitting_button, 2], [self.proceed_button, 3]])
        structure.append([])

        self.set_close_button_function(self.trigger_exit)

        self.setup_window(structure, main_font='Courier')

        self.stop_and_return_button.activate()
        self.return_to_fitting_button.disable()
        self.proceed_button.disable()

    def run_fitting(self):

        self.log.set_param('proceed', 'return_to_fitting')

        self.round = 0
        self.max_round = int(self.iterations / 1000)

        self.progress_fitting.initiate(self.max_round)

        self.after(self.measure)

    def measure(self):

        if self.exit:
            self.after(self.save)

        else:
            if self.round == 0:
                self.mcmc_fit.run_mcmc()
            else:
                self.mcmc_fit.rerun_mcmc()

            self.round += 1
            self.progress_fitting.update()

            if self.round >= self.max_round:
                self.progress_fitting.show_message('Saving results...')
                self.after(self.save)
            else:
                self.after(self.measure)

    def save(self):

        if self.exit:
            self.after(self.options)

        else:

            self.mcmc_fit.get_results()

            fitting_directory = '.'.join(self.light_curve_file.split('.')[:-1]) + '_' + self.log.fitting_directory_base

            if not os.path.isdir(fitting_directory):
                os.mkdir(fitting_directory)
            else:
                fi = 2
                while os.path.isdir('{0}_{1}'.format(fitting_directory, str(fi))):
                    fi += 1
                fitting_directory = '{0}_{1}'.format(fitting_directory, str(fi))
                os.mkdir(fitting_directory)

            # saving

            shutil.copy(self.log.files['fitting_output_description'], fitting_directory)

            self.mcmc_fit.results['parameters']['mt']['print_name'] += '(BJD_{TDB})'

            self.mcmc_fit.save_results(os.path.join(fitting_directory, 'results.txt'))
            self.mcmc_fit.save_models(os.path.join(fitting_directory, 'full_model.txt'))
            self.mcmc_fit.save_detrended_models(os.path.join(fitting_directory, 'detrended_model.txt'))

            self.mcmc_fit.plot_hops_corner(fitting_directory)
            self.mcmc_fit.plot_hops_output(
                self.title1, self.title2, self.title3, fitting_directory, self.log.logo_jpg)

            self.log.set_param('fitting_complete', True)
            self.log.set_param('fitting_version', self.log.version)

            self.log.save_local_log()

            self.log.export_local_log(fitting_directory)

            self.progress_fitting.show_message(' ')

            self.showinfo('Results saved successfully',
                          'Your results have been saved successfully in {0}'.format(fitting_directory))
            self.show()

            self.after(self.options)

    def options(self):

        if self.exit:
            self.def_close()
        else:
            self.set_close_button_function(self.close)

            self.stop_and_return_button.disable()

            self.return_to_fitting_button.activate()
            self.proceed_button.activate()

    def proceed(self):

        self.log.set_param('proceed', 'exit')
        self.close()







import matplotlib.cm as cm
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
from matplotlib.backend_bases import FigureCanvasBase
try:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
except ImportError:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
    NavigationToolbar2Tk = NavigationToolbar2TkAgg


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

    def plot_hops_output(self, title1, title2, title3, fitting_directory, logo_file):

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

            fig.text(0.5, 0.94, title1, fontsize=24, va='center', ha='center')
            fig.text(0.97, 0.97, title2, fontsize=fsmain, va='top', ha='right')
            fig.text((1 - fright) / fcol, 1 - (1 - fbottom) / frow, title3,
                                                          fontsize=fsmain, ha='left', va='bottom')

            logo_ax = fig.add_subplot(gs[0, 0])
            logo_ax.imshow(logo_file)
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

            fig.savefig(os.path.join(fitting_directory, 'detrended_model.jpg'), dpi=1200, transparent=False)


            fig = Figure(figsize=(funit * fcol / (1 - fright), funit * frow / (1 - fbottom)))
            canvas = FigureCanvasBase(fig)
            try:
                gs = gridspec.GridSpec(frow, fcol, fig, 0, fbottom, 1.0 - fright, 1.0, 0.0, 0.0)
            except TypeError:
                gs = gridspec.GridSpec(frow, fcol, 0, fbottom, 1.0 - fright, 1.0, 0.0, 0.0)

            fig.text(0.5, 0.94, title1, fontsize=24, va='center', ha='center')
            fig.text(0.97, 0.97, title2, fontsize=fsmain, va='top', ha='right')
            fig.text((1 - fright) / fcol, 1 - (1 - fbottom) / frow, title3,
                                                          fontsize=fsmain, ha='left', va='bottom')

            logo_ax = fig.add_subplot(gs[0, 0])
            logo_ax.imshow(logo_file)
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
                     self.results['input_series']['value'][set_indices], 'ko', ms=2)
            ax1.plot(self.results['detrended_output_series']['phase'][set_indices],
                     self.results['output_series']['model'][set_indices], 'r-', label='Transit + Quadratic de-trending')

            fig.text(0.04, fbottom + 2.5 * (1 - fbottom) / frow, 'relative flux (raw)', fontsize=fsbig, va='center',
                     ha='center', rotation='vertical')

            data_ymin = (min(self.results['input_series']['value'][set_indices])
                         - 3 * np.std(self.results['output_series']['residuals'][set_indices]))

            data_ymax = (max(self.results['input_series']['value'][set_indices])
                         + 2 * np.std(self.results['output_series']['residuals'][set_indices]))

            ax1.set_yticks(ax1.get_yticks()[np.where(ax1.get_yticks() > data_ymin)])

            ymin, ymax = data_ymax - 1.05 * (data_ymax - data_ymin), data_ymax

            ax1.set_ylim(ymin, ymax)

            x_max = max(np.abs(self.results['output_series']['phase'][set_indices]) +
                        0.05 * (max(self.results['output_series']['phase'][set_indices]) -
                                min(self.results['output_series']['phase'][set_indices])))

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

            ax2 = fig.add_subplot(gs[4, 1:])
            ax2.plot(self.results['output_series']['phase'][set_indices],
                     self.results['output_series']['residuals'][set_indices], 'ko', ms=2)
            ax2.plot(self.results['output_series']['phase'][set_indices],
                     np.zeros_like(self.results['output_series']['phase'][set_indices]), 'r-')

            ax2.set_ylim(- 5 * np.std(self.results['output_series']['residuals'][set_indices]),
                         5 * np.std(self.results['output_series']['residuals'][set_indices]))

            ax2.set_xlabel('phase', fontsize=fsbig)
            fig.text(0.04, fbottom + 0.5 * (1 - fbottom) / frow, 'residuals', fontsize=fsbig, va='center', ha='center',
                     rotation='vertical')

            ax2.set_xlim(-x_max, x_max)
            ax2.tick_params(labelsize=fsmain)

            ax2.text(ax2.get_xlim()[0] + 0.02 * (ax2.get_xlim()[-1] - ax2.get_xlim()[0]),
                     ax2.get_ylim()[0] + 0.07 * (ax2.get_ylim()[-1] - ax2.get_ylim()[0]),
                     r'$\mathrm{rms}_\mathrm{res} = %.1e$' %
                     np.std(self.results['output_series']['residuals'][set_indices]),
                     fontsize=fsmain)

            ax1.legend(loc=1, fontsize='x-small')

            fig.savefig(os.path.join(fitting_directory, 'full_model.jpg'), dpi=1200, transparent=False)

