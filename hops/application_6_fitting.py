
import os
import glob
import numpy as np
import shutil
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
import hops.pylightcurve41 as plc

from matplotlib.figure import Figure
from matplotlib.backend_bases import MouseEvent as mpl_MouseEvent

from hops.application_windows import MainWindow, AddOnWindow

filter_map = {'Clear': 'clear',
              'Luminance': 'luminance',
              'U': 'JOHNSON_U',
              'B': 'JOHNSON_B',
              'V': 'JOHNSON_V',
              'R': 'COUSINS_R',
              'I': 'COUSINS_I',
              'H': '2mass_h',
              'J': '2mass_j',
              'K': '2mass_ks',
              'Astrodon ExoPlanet-BB': 'exoplanets_bb',
              'u\'': 'sdss_u',
              'g\'': 'sdss_g',
              'r\'': 'sdss_r',
              'z\'': 'sdss_z',
              }

__location__ = os.path.abspath(os.path.dirname(__file__))


class FittingWindow(MainWindow):

    def __init__(self, log):

        # define windows

        MainWindow.__init__(self, log, name='HOPS - Fitting', position=2)
        self.export_window = AddOnWindow(self, name='Export Results', position=1)

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
                                              width=40, command=self.update_lc)

        self.iterations = self.Entry(value=self.log.get_param('iterations'), instance=int)
        self.burn = self.Entry(value=self.log.get_param('burn'), instance=int)
        self.a_i_fit = self.CheckButton(initial=self.log.get_param('a_i_fit'), text='Fit for a/Rs, i.',
                                        command=self.update_planet)
        self.detrending = self.DropDown(initial=self.log.get_param('detrending'),
                                        options=['Airmass', 'Quadratic', 'Linear'],
                                        instance=str, command=self.update_planet)
        self.manual_planet = self.CheckButton(initial=self.log.get_param('manual_planet'),
                                              text='Enter param. manually', command=self.update_planet)

        self.planet = self.Entry(value=self.log.get_param('planet'), command=self.update_planet)
        self.target_ra_dec = self.Entry(value=self.log.get_param('target_ra_dec'), command=self.update_planet)
        self.metallicity = self.Entry(value=self.log.get_param('metallicity'), instance=float, command=self.update_planet)
        self.temperature = self.Entry(value=self.log.get_param('temperature'), instance=float, command=self.update_planet)
        self.logg = self.Entry(value=self.log.get_param('logg'), instance=float, command=self.update_planet)
        self.period = self.Entry(value=self.log.get_param('period'), instance=float, command=self.update_planet)
        self.mid_time = self.Entry(value=self.log.get_param('mid_time'), instance=float, command=self.update_planet)
        self.rp_over_rs = self.Entry(value=self.log.get_param('rp_over_rs'), instance=float, command=self.update_planet)
        self.sma_over_rs = self.Entry(value=self.log.get_param('sma_over_rs'), instance=float, command=self.update_planet)
        self.inclination = self.Entry(value=self.log.get_param('inclination'), instance=float, command=self.update_planet)
        self.eccentricity = self.Entry(value=self.log.get_param('eccentricity'), instance=float, command=self.update_planet)
        self.periastron = self.Entry(value=self.log.get_param('periastron'), instance=float, command=self.update_planet)

        try:
            ra_target, dec_target = self.log.get_param('target_ra_dec').split(' ')
            planet = plc.locate_planet(plc.Hours(ra_target), plc.Degrees(dec_target))
            ecc_data = planet.all_data
            self.auto_planet = self.Label(text=planet.name)
            self.auto_target_ra_dec = self.Label(text='{0} {1}'.format(ecc_data['star']['ra'], ecc_data['star']['dec']))
            self.auto_metallicity = self.Label(text=ecc_data['planet']['meta'], instance=float)
            self.auto_temperature = self.Label(text=ecc_data['planet']['teff'], instance=float)
            self.auto_logg = self.Label(text=ecc_data['planet']['logg'], instance=float)
            self.auto_period = self.Label(text=ecc_data['planet']['ephem_period'], instance=float)
            self.auto_mid_time = self.Label(text=ecc_data['planet']['ephem_mid_time'], instance=float)
            self.auto_rp_over_rs = self.Label(text=ecc_data['planet']['rp_over_rs'], instance=float)
            self.auto_sma_over_rs = self.Label(text=ecc_data['planet']['sma_over_rs'], instance=float)
            self.auto_eccentricity = self.Label(text=ecc_data['planet']['eccentricity'], instance=float)
            self.auto_inclination = self.Label(text=ecc_data['planet']['inclination'], instance=float)
            self.auto_periastron = self.Label(text=ecc_data['planet']['periastron'], instance=float)
        except:
            planet = None
            self.auto_planet = self.Label(text='None')
            self.auto_target_ra_dec = self.Label(text='None')
            self.auto_metallicity = self.Label(text=0, instance=float)
            self.auto_temperature = self.Label(text=0, instance=float)
            self.auto_logg = self.Label(text=0, instance=float)
            self.auto_period = self.Label(text=0, instance=float)
            self.auto_mid_time = self.Label(text=0, instance=float)
            self.auto_rp_over_rs = self.Label(text=0, instance=float)
            self.auto_sma_over_rs = self.Label(text=0, instance=float)
            self.auto_eccentricity = self.Label(text=0, instance=float)
            self.auto_inclination = self.Label(text=0, instance=float)
            self.auto_periastron = self.Label(text=0, instance=float)

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

        self.test_button = self.Button(text='RUN TEST', command=[self.disable, self.run_test],
                                          bg='green', highlightbackground='green')
        self.fitting_button = self.Button(text='RUN FITTING',
                                          command=[self.disable, self.save_settings, self.run_fitting, self.activate],
                                          bg='green', highlightbackground='green')

        self.export_button = self.Button(text='EXPORT FOR DATABASES (ExoClock / ETD)', command=self.export_window.show)

        # define preview window

        funit = 0.8
        fcol = 7
        frow = 5
        fbottom = 0.11
        fright = 0.02
        self.fssmall = 9
        self.fsmain = 9
        self.fsbig = 20
        self.preview_figure = self.FigureWindow(figsize=(6, 7), show_nav=True)
        self.preview_figure.figure.canvas.callbacks.connect('button_press_event', self.show_point_info)
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
        self.preview_figure.figure.text(0.04, fbottom + 2.5 * (1 - fbottom) / frow, 'Relative flux (norm)', fontsize=self.fsmain, va='center',
                 ha='center', rotation='vertical')
        self.preview_figure.figure.text(0.04, fbottom + 0.5 * (1 - fbottom) / frow, 'Residuals (norm)', fontsize=self.fsmain, va='center',
                 ha='center', rotation='vertical')

        self.title1 = self.preview_figure.figure.text(0.5, 0.94, '', fontsize=self.fsbig, va='center', ha='center')
        self.title2 = self.preview_figure.figure.text(0.97, 0.97, '', fontsize=self.fsmain, va='top', ha='right')
        self.title3 = self.preview_figure.figure.text(0.03 + (1 - fright) / fcol, 1 - (1 - fbottom) / frow, '', fontsize=self.fsmain, ha='left', va='bottom')
        self.dblclick = False

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

            [[self.Label(text='MCMC options:'), 1, 2], [self.Label(text='De-trending:'), 3]],

            [[self.Label(text='Iterations (def=5000)'), 1], [self.iterations, 2], [self.detrending, 3]],
            [[self.Label(text='Burn-in (def=1000)'), 1], [self.burn, 2], [self.a_i_fit, 3]],

            [[self.Button(text='RETURN TO\nMAIN MENU', command=self.close), 1],
             [self.Button(text='SAVE OPTIONS &\nRETURN TO MAIN MENU', command=[self.save_settings, self.close]), 2],
             [self.Button(text='RETURN TO\nPHOTOMETRY', command=self.return_to_photometry), 3]
             ],
            [[self.test_button, 1],[self.fitting_button, 2]],
            [],

        ])

        self.update_lc()
        self.update_planet()

        # extra windows

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

    def update_lc(self, *entry):

        if not entry:
            pass

        self.fitting_button.disable()

        if not os.path.isfile(self.light_curve_file.get()):

            self.export_button.disable()
            self.iterations.disable()
            self.burn.disable()
            self.a_i_fit.disable()
            self.detrending.disable()
            self.fitting_button.disable()
            self.test_button.disable()

            self.manual_planet.disable()
            self.target_ra_dec.disable()
            self.planet.disable()
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

            self.ax1.cla()
            self.ax2.cla()
            self.title2.set_text('')
            self.title3.set_text('')

            self.preview_figure.draw(update_level=2)

        else:
            self.export_button.activate()
            self.iterations.activate()
            self.burn.activate()
            self.a_i_fit.activate()
            self.detrending.activate()
            self.test_button.activate()

            self.manual_planet.activate()
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
                self.inclination.activate()
                self.eccentricity.activate()
                self.periastron.activate()
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

            light_curve = np.loadtxt(self.light_curve_file.get(), unpack=True)

            self.ax1.cla()
            self.ax2.cla()
            self.plot_data = self.ax1.errorbar(light_curve[0], light_curve[1], light_curve[2], c='k', fmt='o', ms=2, lw=0.5, zorder=0)
            self.plot_data = self.plot_data.lines[0].get_xdata()

            date = plc.JD(light_curve[0][0]).utc.isoformat()[:16].replace('T', ' ')
            self.observation_date = plc.JD(light_curve[0][0]).utc.isoformat()[:16].split('T')[0]
            obs_duration = round(24 * (light_curve[0][-1] - light_curve[0][0]), 1)

            self.title2.set_text('{0} (UT)\nDur: {1}h / Exp: {2}s\nFilter: {3}'.format(date, obs_duration, self.exposure_time,
                                                                                       self.log.get_param('filter')))
            self.title3.set_text(
                '\n\n{0}\n{1}'.format(
                    self.log.get_param('observer'), '{0} / {1} / {2}'.format(
                        self.log.get_param('observatory'), self.log.get_param('telescope'), self.log.get_param('camera'))))

            self.preview_figure.draw(update_level=2)

    def update_planet(self, *entry):

        if not entry:
            pass

        self.fitting_button.disable()

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

            self.title1.set_text('{0}{1}{2}'.format('$\mathbf{', self.planet.get(), '}$'))

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

            self.title1.set_text('{0}{1}{2}'.format('$\mathbf{', self.auto_planet.get(), '}$'))

        self.preview_figure.draw()

    def run_test(self):
        try:
            planet, fitting = self.call_plc_fitting(optimiser='curve_fit')
            self.plot_results_to_axes(planet, fitting, self.ax1, self.ax2, self.title1, self.title2, self.title3, 'Raw')
            self.preview_figure.draw()
            self.activate()
            self.manual_planet.activate()
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
                self.inclination.activate()
                self.eccentricity.activate()
                self.periastron.activate()
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

        except:
            if self.log.get_param('location') == '+dd:mm:ss +dd:mm:ss':
                self.showinfo('Test failed', 'Test failed! Please update the location of your telescope '
                                             'in the SELECT DATA & TARGET window.')
            elif self.auto_planet.get() == 'None':
                self.showinfo('Test failed', 'Test failed! Please update the target coordinates in the '
                                             'SELECT DATA & TARGET window or add the planet parameters manually.')
            else:
                self.showinfo('Test failed', 'Test failed! Please change the input parameters.')
            self.activate()
            self.fitting_button.disable()

    def run_fitting(self):

        planet, fitting = self.call_plc_fitting(optimiser='emcee')

        if fitting.mcmc_run_complete:

            fitting_directory = '.'.join(self.light_curve_file.get().split('.')[:-1]) + '_' + self.log.fitting_directory_base

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

            fitting.save_results(os.path.join(fitting_directory, 'results.txt'))
            fitting.plot_corner(os.path.join(fitting_directory, 'corner.pdf'))
            fitting.plot_traces(os.path.join(fitting_directory, 'traces.pdf'))

            np.savetxt(
                os.path.join(fitting_directory, 'model.txt'),
                np.swapaxes(
                    [
                        fitting.results['input_series']['time'],
                        (fitting.results['input_series']['time'] - fitting.results['parameters']['T_mid']['value'])/fitting.results['parameters']['P']['value'],
                        fitting.results['input_series']['flux'],
                        fitting.results['input_series']['flux_unc'],
                        fitting.results['output_series']['model'],
                        fitting.results['output_series']['residuals'],
                        ], 0, 1
                )
            )

            np.savetxt(
                os.path.join(fitting_directory, 'detrended_model.txt'),
                np.swapaxes(
                    [
                        fitting.results['input_series']['time'],
                        (fitting.results['input_series']['time'] - fitting.results['parameters']['T_mid']['value'])/fitting.results['parameters']['P']['value'],
                        fitting.results['detrended_input_series']['flux'],
                        fitting.results['detrended_input_series']['flux_unc'],
                        fitting.results['detrended_output_series']['model'],
                        fitting.results['detrended_output_series']['residuals'],
                        ], 0, 1
                )
            )

            w = open(os.path.join(fitting_directory, 'results.txt'), 'a')

            w.write('\n\n#Detrended Residuals:\n')
            w.write('#Mean: {0}\n'.format(fitting.results['detrended_statistics']['res_mean']))
            w.write('#STD: {0}\n'.format(fitting.results['detrended_statistics']['res_std']))
            w.write('#RMS: {0}\n'.format(fitting.results['detrended_statistics']['res_rms']))

            w.close()

            funit = 0.8
            fcol = 7
            frow = 5
            fbottom = 0.11
            fright = 0.02
            fssmall = 8
            fsmain = 10
            fsbig = 20
            figure = Figure(figsize=(funit * fcol / (1 - fright), funit * frow / (1 - fbottom)))
            try:
                gs = gridspec.GridSpec(frow, fcol, figure, 0.03, fbottom, 1.0 - fright, 1.0, 0.0, 0.0)
            except TypeError:
                gs = gridspec.GridSpec(frow, fcol, 0.03, fbottom, 1.0 - fright, 1.0, 0.0, 0.0)

            logo_ax = figure.add_subplot(gs[0, 0])
            logo_ax.imshow(mpimg.imread(self.log.files['holomon_logo']))
            logo_ax.spines['top'].set_visible(False)
            logo_ax.spines['bottom'].set_visible(False)
            logo_ax.spines['left'].set_visible(False)
            logo_ax.spines['right'].set_visible(False)
            logo_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

            ax1 = figure.add_subplot(gs[1:4, 1:])
            ax2 = figure.add_subplot(gs[4, 1:])
            figure.text(0.04, fbottom + 2.5 * (1 - fbottom) / frow, 'Relative flux (norm)', fontsize=fsmain, va='center',
                        ha='center', rotation='vertical')
            figure.text(0.04, fbottom + 0.5 * (1 - fbottom) / frow, 'Residuals (norm)', fontsize=fsmain, va='center',
                        ha='center', rotation='vertical')

            title1 = figure.text(0.5, 0.94, '', fontsize=self.fsbig, va='center', ha='center')
            title2 = figure.text(0.97, 0.97, '', fontsize=self.fsmain, va='top', ha='right')
            title3 = figure.text(0.03 + (1 - fright) / fcol, 1 - (1 - fbottom) / frow, '', fontsize=self.fsmain, ha='left', va='bottom')

            self.plot_results_to_axes(planet, fitting, ax1, ax2, title1, title2, title3, 'De-trended')

            figure.savefig(os.path.join(fitting_directory, 'detrended_model.jpg'), dpi=1200, transparent=False)

            self.log.set_param('fitting_complete', True)
            self.log.set_param('fitting_version', self.log.version)

            self.log.save_local_log()

            self.log.export_local_log(fitting_directory)

            self.showinfo('Results saved successfully',
                          'Your results have been saved successfully in {0}'.format(fitting_directory))

    def call_plc_fitting(self, optimiser):

        if self.manual_planet.get():
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

            planet_to_plot = self.planet.get()

        else:
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

            planet_to_plot = self.auto_planet.get()

        light_curve = np.loadtxt(self.light_curve_file.get(), unpack=True)

        location_string = self.log.get_param('location').split(' ')
        observatory = plc.Observatory(plc.Degrees(location_string[0]), plc.Degrees(location_string[1]))

        ra_dec_string = target_ra_dec_to_plot.split(' ')
        planet = plc.Planet(
            planet_to_plot, plc.Hours(ra_dec_string[0]), plc.Degrees(ra_dec_string[1]),
            logg_to_plot, temperature_to_plot, metallicity_to_plot,
            rp_over_rs_to_plot, period_to_plot, sma_over_rs_to_plot, eccentricity_to_plot, inclination_to_plot,
            periastron_to_plot, mid_time_to_plot, 'BJD_TDB')

        planet.add_observation(
            time=light_curve[0],
            time_format='JD_UTC',
            exp_time=float(self.exposure_time),
            time_stamp='start',
            flux=light_curve[1],
            flux_unc=light_curve[2],
            flux_format='flux',
            filter_name=filter_map[self.log.get_param('filter')],
            observatory_latitude=observatory.latitude,
            observatory_longitude=observatory.longitude,
        )

        if self.a_i_fit.get():
            fit_sma_over_rs=True
            fit_inclination=True
        else:
            fit_sma_over_rs=False
            fit_inclination=False

        if self.detrending.get() == 'Airmass':
            detrending = 'airmass'
            detrending_order = 1
        elif self.detrending.get() == 'Linear':
            detrending = 'time'
            detrending_order = 1
        else:
            detrending = 'time'
            detrending_order = 2

        return planet, planet.transit_fitting(
            detrending_order=detrending_order, detrending=detrending,
            optimise_initial_parameters=True, scale_uncertainties=True, filter_outliers=True,
            fit_sma_over_rs=fit_sma_over_rs, fit_inclination=fit_inclination,
            counter=None, window_counter=True,
            iterations=self.iterations.get(), burn_in=self.burn.get(),
            optimiser=optimiser,
        )

    def plot_results_to_axes(self, planet, fitting, ax1, ax2, title1, title2, title3, lc_type):

        light_curve = np.loadtxt(self.light_curve_file.get(), unpack=True)
        date = plc.JD(light_curve[0][0]).utc.isoformat()[:16].replace('T', ' ')
        obs_duration = round(24 * (light_curve[0][-1] - light_curve[0][0]), 1)

        title1.set_text('{0}{1}{2}'.format('$\mathbf{', planet.name, '}$'))

        title2.set_text('{0} (UT)\nDur: {1}h / Exp: {2}s\nFilter: {3}'.format(date, obs_duration, self.exposure_time,
                                                                                   self.log.get_param('filter')))
        title3.set_text(
            '\n\n{0}\n{1}'.format(
                self.log.get_param('observer'), '{0} / {1} / {2}'.format(
                    self.log.get_param('observatory'), self.log.get_param('telescope'), self.log.get_param('camera'))))

        new_mid_time = fitting.results['parameters']['T_mid']['value']
        time = fitting.results['input_series']['time']
        phase = (time - new_mid_time)/planet.period

        predicted_transit_model = planet.transit_integrated(
            time, 'BJD_TDB', float(self.exposure_time), 'start', filter_map[self.log.get_param('filter')]
        )
        epoch = round((new_mid_time-planet.mid_time)/planet.period)
        predicted_mid_time = planet.mid_time + epoch * planet.period
        predicted_rp_over_rs = planet.rp_over_rs

        ax1.cla()

        if lc_type == 'Raw':
            flux = fitting.results['input_series']['flux']
            flux_unc = fitting.results['input_series']['flux_unc']
            model = fitting.results['output_series']['model']
            detrended_model = fitting.results['detrended_output_series']['model']
            trend = model / detrended_model
            std_res = fitting.results['statistics']['res_std']
        else:
            flux = fitting.results['detrended_input_series']['flux']
            flux_unc = fitting.results['detrended_input_series']['flux_unc']
            model = fitting.results['detrended_output_series']['model']
            trend = 1
            std_res = fitting.results['detrended_statistics']['res_std']

        self.plot_data = ax1.errorbar(phase, flux, flux_unc, c='k', fmt='o', ms=2, lw=0.5, zorder=0)
        self.plot_data = self.plot_data.lines[0].get_xdata()
        ax1.plot(phase, flux, 'ko', ms=2, label=lc_type + ' data', zorder=0)

        if lc_type == 'Raw':
            ax1.errorbar(
                (fitting.results['settings']['time'][np.where(fitting.results['prefit']['outliers_map'])]- new_mid_time)/planet.period,
                fitting.results['settings']['flux'][np.where(fitting.results['prefit']['outliers_map'])],
                fitting.results['settings']['flux_unc'][np.where(fitting.results['prefit']['outliers_map'])] * fitting.results['prefit']['scale_factor'],
                c='k', fmt='o', ms=2, lw=0.5, zorder=0)
            ax1.plot(
                (fitting.results['settings']['time'][np.where(fitting.results['prefit']['outliers_map'])]- new_mid_time)/planet.period,
                fitting.results['settings']['flux'][np.where(fitting.results['prefit']['outliers_map'])],
                'rx', ms=7, label='Outliers (not fitted)', zorder=1)

        data_ymin = min(flux) - 5 * std_res
        data_ymax = max(flux) + 2 * std_res

        ax1.set_yticks(ax1.get_yticks()[np.where(ax1.get_yticks() > data_ymin)])
        ymin, ymax = data_ymax - 1.3 * (data_ymax - data_ymin), data_ymax
        ax1.set_ylim(ymin, ymax)
        x_max = max(np.abs(phase) + 0.05 * (max(phase) - min(phase)))

        ax1.set_xlim(-x_max, x_max)
        ax1.tick_params(labelbottom=False, labelsize=self.fsmain)

        ax1.plot(phase, model, 'r-', label='Best-fit model (De-trending: {11}, a/Rs,i: {12})\n{0}{3}{0}={0}{4}_{1}-{5}{2}^{1}+{6}{2}{0}, {0}{7}{0}={0}{8}_{1}-{9}{2}^{1}+{10}{2}{0}'.format(
            '$', '{', '}',
            fitting.results['parameters']['T_mid']['print_name'],
            fitting.results['parameters']['T_mid']['print_value'],
            fitting.results['parameters']['T_mid']['print_m_error'],
            fitting.results['parameters']['T_mid']['print_p_error'],
            fitting.results['parameters']['rp']['print_name'],
            fitting.results['parameters']['rp']['print_value'],
            fitting.results['parameters']['rp']['print_m_error'],
            fitting.results['parameters']['rp']['print_p_error'],
            self.detrending.get(),
            ['Fixed', 'Free'][self.a_i_fit.get()],
        ),zorder=3)

        ax1.plot(phase, predicted_transit_model * trend, 'c-', label='Expected model\n{0}{3}{0}={0}{4}{0}, {0}{5}{0}={0}{6}{0}, O-C={0}{7}_{1}-{8}{2}^{1}+{9}{2}{0} min'.format(
            '$','{', '}',
            fitting.results['parameters']['T_mid']['print_name'],
            round(predicted_mid_time, 5),
            fitting.results['parameters']['rp']['print_name'],
            round(predicted_rp_over_rs, 5),
            round((new_mid_time - predicted_mid_time) * 24 * 60, 1),
            round(fitting.results['parameters']['T_mid']['m_error'] * 24 * 60, 1),
            round(fitting.results['parameters']['T_mid']['p_error'] * 24 * 60, 1),
        ), zorder=2)

        ax1.legend(loc=3, fontsize=self.fssmall)

        ax2.cla()

        detrended_flux_unc = fitting.results['detrended_input_series']['flux_unc']
        detrended_residuals = fitting.results['detrended_output_series']['residuals']
        detrended_std_res = fitting.results['detrended_statistics']['res_std']
        res_autocorr = fitting.results['statistics']['res_max_autocorr']
        res_autocorr_flag = fitting.results['statistics']['res_max_autocorr_flag']

        ax2.errorbar(phase, detrended_residuals,  detrended_flux_unc, c='k', fmt='o', ms=2, lw=0.5)
        ax2.plot(phase, detrended_residuals, 'ko', ms=2)
        ax2.plot(phase, np.zeros_like(phase), 'r-')

        ax2.set_ylim(- 6 * detrended_std_res, 6 * detrended_std_res)

        ax2.set_xlabel('Phase', fontsize=self.fsmain)

        ax2.set_xlim(-x_max, x_max)
        ax2.tick_params(labelsize=self.fsmain)

        ax2.text(ax2.get_xlim()[0] + 0.02 * (ax2.get_xlim()[-1] - ax2.get_xlim()[0]),
                 ax2.get_ylim()[1] - 0.15 * (ax2.get_ylim()[-1] - ax2.get_ylim()[0]),
                 r'STD = %.2f $‰$' % (detrended_std_res * 1000),
                 fontsize=self.fssmall)
        ax2.text(ax2.get_xlim()[0] + 0.02 * (ax2.get_xlim()[-1] - ax2.get_xlim()[0]),
                 ax2.get_ylim()[0] + 0.07 * (ax2.get_ylim()[-1] - ax2.get_ylim()[0]),
                 r'AutoCorr = %.2f' %res_autocorr + [' (acceptable)', ' (too high)'][res_autocorr_flag],
                 fontsize=self.fssmall)

    def show_point_info(self, event=None):
        if isinstance(event, mpl_MouseEvent):
            if event.inaxes:
                if event.dblclick and not self.dblclick:
                    if event.y > 0.33 * self.preview_figure.canvas.get_tk_widget().winfo_reqheight():
                        idx = np.argmin((event.xdata - self.plot_data) ** 2)
                        self.showinfo('Point information', self.science_files[idx])
                    self.dblclick = True
                else:
                    self.dblclick = False

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
        self.log.set_param('iterations', self.iterations.get())
        self.log.set_param('burn', self.burn.get())
        self.log.set_param('a_i_fit', self.a_i_fit.get())
        self.log.set_param('detrending', self.detrending.get())
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

    def return_to_photometry(self):

        self.save_settings()
        self.log.set_param('proceed', 'return_to_photometry')
        self.close()

    def export(self):

        if self.manual_planet.get():
            planet_to_plot = self.planet.get()
        else:
            planet_to_plot = self.auto_planet.get()

        light_curve = np.loadtxt(self.light_curve_file.get(), unpack=True)

        file_name = 'HOPS_{5}_{0}_{1}_{2}_{3}s_for_{4}.txt'.format(
            self.observation_date,
            planet_to_plot,
            self.log.get_param('filter'),
            self.exposure_time,
            self.export_window_database.get(),
            self.light_curve_file.get().replace(os.sep, '_').replace('.txt', '')
        )

        w = open(file_name, 'w')

        w.write(self.export_window_header)

        for i in range(len(light_curve[0])):
            w.write('{0:.10f}\t{1:.10f}\t{2:.10f}\n'.format(light_curve[0][i] + self.export_window_time_shift,
                                             light_curve[1][i],
                                             light_curve[2][i] * (1 / np.sqrt(self.export_window_camera_gain.get())),
                                             ))

        self.showinfo('Export', 'File saved at: ' + file_name)

    def export_window_update(self):

        if self.manual_planet.get():
            planet_to_plot = self.planet.get()
        else:
            planet_to_plot = self.auto_planet.get()

        if self.export_window_database.get() == 'ExoClock':

            self.export_window_camera_gain.set(1.0)
            self.export_window_time_shift = 0

            self.export_window_header = '\n'.join([
                '#Planet: {0}'.format(planet_to_plot),
                '#Time format: JD_UTC',
                '#Time stamp: Exposure start',
                '#Flux format: Flux',
                '#Filter: {0}'.format(self.log.get_param('filter')),
                '#Exposure time in seconds: {0}'.format(self.exposure_time),
                '#Time                   Flux            Flux uncertainty\n'
            ])

            self.export_window_message.set('No transformation is needed to upload your results to\nExoClock ' \
                                         '(www.exoclock.space).\n\nWhen uploading you will need to use the following ' \
                                         'details:\n' \
                                         'Planet: {0}\n'\
                                         'Time format: JD_UTC\n'\
                                         'Time stamp: Exposure start\n'\
                                         'Flux format: Flux\n'\
                                         'Filter: {1}\n'\
                                         'Exposure time in seconds: {2}'.format(planet_to_plot,
                                                                                self.log.get_param('filter'),
                                                                                self.exposure_time))

        else:
            self.export_window_time_shift = self.exposure_time / 60 / 60 / 24 / 2

            self.export_window_header = '\n'.join([
                '#Planet: {0}'.format(planet_to_plot),
                '#Time format: JD_UTC (geocentric)',
                '#Time stamp: Mid-exposure',
                '#Flux format: Flux (in flux)',
                '#Filter: {0}'.format(self.log.get_param('filter')),
                '#Exposure time in seconds: {0}'.format(self.exposure_time),
                '#Time                   Flux            Flux uncertainty\n'
            ])

            self.export_window_message.set('To upload your results to ETD (var2.astro.cz/ETD)\nthe following modifications ' \
                                         'are applied:\n\n' \
                                         '1. Half of the exposure time is added as the time-stamp required on ETD ' \
                                         'referes to the exposure mid-time\n' \
                                         '2. The uncertainties are  mupliplied by 1/SQRT(gain) to represent the ' \
                                         'uncertainies in electrons rather than counts.\n\n'\
                                         'When uploading you will need to use the following details:\n' \
                                         'Planet: {0}\n'\
                                         'JD format: geocentric\n'\
                                         'Brightness column: in flux'.format(planet_to_plot))
