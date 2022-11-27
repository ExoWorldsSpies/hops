__all__ = ['Planet', 'get_filter', 'add_filter']

import os
import shutil
import warnings

import numpy as np

from pylightcurve.errors import *
from pylightcurve.__databases__ import plc_data

from pylightcurve.models.exoplanet_lc import planet_orbit, planet_star_projected_distance, planet_phase, \
    transit, transit_integrated, \
    transit_duration, transit_depth, eclipse, eclipse_integrated, eclipse_mid_time, eclipse_duration, eclipse_depth,\
    fp_over_fs, transit_t12, exotethys
from pylightcurve.analysis.optimisation import Fitting
from pylightcurve.processes.files import open_dict, save_dict,  copy_dict
from pylightcurve.plots.plots_fitting import plot_transit_fitting_models
from pylightcurve.spacetime.angles import Degrees, _request_angle, _reformat_or_request_angle
from pylightcurve.spacetime.times import JD
from pylightcurve.spacetime.targets import FixedTarget
from pylightcurve.plots.plots_fitting import in_brackets

from pylightcurve.spacetime.times import now
from pylightcurve.spacetime.observing import Observatory

from inspect import signature


class Filter:

    def __init__(self, name, passband=None):

        if name in plc_data.all_filters():
            photometry_path = plc_data.photometry()
            passband = np.loadtxt(os.path.join(photometry_path, name + '.pass'))
        else:
            try:
                passband = np.loadtxt(passband)
            except:
                raise PyLCInputError('Wrong passband format or file path.')

        self.name = name
        self.passband = passband
        self.ldcs_mem = {}
        self.fp_over_fs_mem = {}

    def ldcs(self, stellar_logg, stellar_temperature, stellar_metallicity, wlrange=None, method='claret', stellar_model='Phoenix_2018'):

        if wlrange is None:
            ldcs_id = '_'.join([str(ff) for ff in [float(stellar_logg), float(stellar_temperature), float(stellar_metallicity), method, stellar_model, 'None']])
        else:
            ldcs_id = '_'.join([str(ff) for ff in [float(stellar_logg), float(stellar_temperature), float(stellar_metallicity), method, stellar_model, float(wlrange[0]), float(wlrange[1])]])

        try:
            return self.ldcs_mem[ldcs_id]
        except KeyError:
            pass

        calc = exotethys(
            stellar_logg, stellar_temperature, stellar_metallicity, self.passband, wlrange, method, stellar_model
        )

        self.ldcs_mem[ldcs_id] = calc
        return calc

    def fp_over_fs(self, rp_over_rs, sma_over_rs, albedo, emissivity, stellar_temperature, wlrange=None):

        if wlrange is None:
            fpfs_id = '_'.join([str(ff) for ff in [rp_over_rs, sma_over_rs, albedo, emissivity, stellar_temperature] + ['None']])
        else:
            fpfs_id = '_'.join([str(ff) for ff in [rp_over_rs, sma_over_rs, albedo, emissivity, stellar_temperature] + wlrange])

        try:
            return self.fp_over_fs_mem[fpfs_id]
        except KeyError:
            pass

        culc = fp_over_fs(rp_over_rs, sma_over_rs, albedo, emissivity, stellar_temperature, self.passband, wlrange)

        self.fp_over_fs_mem[fpfs_id] = culc
        return culc


class Trend:

    def __init__(self, function, trend_priors=None):

        self.function = function
        self.coefficients = len(str(signature(function))[1:-1].split(','))

        self.names = ['c{0}'.format(ff) for ff in range(self.coefficients)]
        self.print_names = ['c{0}'.format(ff) for ff in range(self.coefficients)]

        if trend_priors is None:
            self.limits1 = [0] + [-2 for ff in range(self.coefficients - 1)]
        else:
            self.limits1 = [0] + list([ff[0] for ff in trend_priors])

        if trend_priors is None:
            self.limits2 = [2] + [2 for ff in range(self.coefficients - 1)]
        else:
            self.limits2 = [2] + list([ff[2] for ff in trend_priors])

        if trend_priors is None:
            self.initial = [1] + [0 for ff in range(self.coefficients - 1)]
        else:
            self.initial = [1] + list([ff[1] for ff in trend_priors])

    def adjust_normalisation_factor(self, flux):

        df = max(np.median(flux) - np.min(flux), np.max(flux) - np.median(flux))
        self.initial[0] = np.median(flux)
        self.limits1[0] = np.min(flux) - 2 * df
        self.limits2[0] = np.max(flux) + 2 * df

    def evaluate(self, auxiliary_data, *coeff):

        coefficients = list(coeff)

        return coefficients[0] * self.function(auxiliary_data, *coefficients[1:])


built_in_filters = {ff: Filter(ff) for ff in plc_data.all_filters()}


def get_filter(filter_name):
    if filter_name in built_in_filters:
        return built_in_filters[filter_name]
    else:
        raise PyLCInputError('{0} is not available. Available filters are: {1}. '
                             '\nAlternatively you can define your own filter as follows: plc.add_filter(filter_name, passbandfile). '
                             '\nThe passband file should contain the passband of the filter '
                             '(two columns txt file, column 1: wavelength in A, column 2: total throughput in electrons/photons)'.format(
            filter_name, ','.join(plc_data.all_filters())))


def add_filter(filter_name, passband):

    built_in_filters[filter_name] = Filter(filter_name, passband)


class Planet:

    def __init__(self, name, ra, dec, stellar_logg, stellar_temperature, stellar_metallicity,
                 rp_over_rs, period, sma_over_rs, eccentricity, inclination,
                 periastron, mid_time, mid_time_format, ww=0,
                 albedo=0.15, emissivity=1.0):

        ra = _reformat_or_request_angle(ra)
        dec = _reformat_or_request_angle(dec)
        inclination = _reformat_or_request_angle(inclination)
        inclination = inclination.deg()
        periastron = _reformat_or_request_angle(periastron)
        periastron = periastron.deg()

        self.name = name
        self.target = FixedTarget(ra, dec)

        self.stellar_logg = stellar_logg
        self.stellar_temperature = stellar_temperature
        self.stellar_metallicity = stellar_metallicity

        self.rp_over_rs = rp_over_rs
        self.period = period
        self.sma_over_rs = sma_over_rs
        self.eccentricity = eccentricity
        self.inclination = inclination
        self.periastron = periastron
        self.mid_time = self.target.convert_to_bjd_tdb(mid_time, mid_time_format)
        self.mid_time_format = 'BJD_TDB'
        self.eclipse_mid_time = eclipse_mid_time(self.period, self.sma_over_rs, self.eccentricity, self.inclination,
                                                 self.periastron, self.mid_time)
        self.ww = ww
        self.albedo = albedo
        self.emissivity = emissivity

        self.filters = built_in_filters

        self.observations = {}

    # planet calculations

    def ldcs(self, filter_name, wlrange=None, method='claret', stellar_model='Phoenix_2018'):
        return get_filter(filter_name).ldcs(self.stellar_logg, self.stellar_temperature, self.stellar_metallicity, wlrange, method, stellar_model)

    def fp_over_fs(self, filter_name, wlrange=None):
        return get_filter(filter_name).fp_over_fs(self.rp_over_rs, self.sma_over_rs, self.albedo, self.emissivity, self.stellar_temperature, wlrange)

    def planet_orbit(self, time, time_format):
        time = self.target.convert_to_bjd_tdb(time, time_format)
        return planet_orbit(self.period, self.sma_over_rs, self.eccentricity, self.inclination, self.periastron,
                            self.mid_time, time, ww=self.ww)

    def planet_star_projected_distance(self, time, time_format):
        time = self.target.convert_to_bjd_tdb(time, time_format)
        return planet_star_projected_distance(self.period, self.sma_over_rs, self.eccentricity, self.inclination,
                                              self.periastron, self.mid_time, time)

    def planet_phase(self, time, time_format):
        time = np.array([self.target.convert_to_bjd_tdb(ff, time_format) for ff in time])
        return planet_phase(self.period, self.mid_time, time)

    def transit(self, time, time_format, filter_name, wlrange=None, method='claret', stellar_model='Phoenix_2018', precision=3):

        time = self.target.convert_to_bjd_tdb(time, time_format)

        return transit(self.ldcs(filter_name, wlrange, method, stellar_model),
                       self.rp_over_rs, self.period, self.sma_over_rs, self.eccentricity,
                       self.inclination, self.periastron, self.mid_time, time,
                       method=method, precision=precision)

    def transit_integrated(self, time, time_format, exp_time, time_stamp, filter_name, wlrange=None, method='claret', stellar_model='Phoenix_2018', max_sub_exp_time=10, precision=3):

        if time_stamp == 'start':
            time = np.array(time) + 0.5 * exp_time / (60.0 * 60.0 * 24.0)
        elif time_stamp == 'mid':
            time = np.array(time)
        elif time_stamp == 'end':
            time = np.array(time) - 0.5 * exp_time / (60.0 * 60.0 * 24.0)
        else:
            raise PyLCInputError(
                'Not acceptable time stamp {0}. Please choose between "mid", "start", "end".'.format(time_stamp))

        time = self.target.convert_to_bjd_tdb(time, time_format)

        return transit_integrated(self.ldcs(filter_name, wlrange=wlrange, method=method, stellar_model=stellar_model),
                                  self.rp_over_rs, self.period, self.sma_over_rs, self.eccentricity,
                                  self.inclination, self.periastron, self.mid_time, time, exp_time,
                                  max_sub_exp_time=max_sub_exp_time,
                                  method=method, precision=precision)

    def transit_duration(self):
        return transit_duration(self.rp_over_rs, self.period, self.sma_over_rs, self.eccentricity, self.inclination,
                                self.periastron)

    def transit_t12(self):
        return transit_t12(self.rp_over_rs, self.period, self.sma_over_rs, self.eccentricity, self.inclination,
                           self.periastron)

    def transit_depth(self, filter_name, wlrange=None, method='claret', stellar_model='Phoenix_2018', precision=6):

        return transit_depth(self.ldcs(filter_name, wlrange, method, stellar_model),
                             self.rp_over_rs, self.period, self.sma_over_rs, self.eccentricity,
                             self.inclination, self.periastron,
                             method=method, precision=precision)

    def eclipse(self, time, time_format, filter_name, wlrange=None, precision=3):

        time = self.target.convert_to_bjd_tdb(time, time_format)

        return eclipse(self.fp_over_fs(filter_name, wlrange), self.rp_over_rs,
                       self.period, self.sma_over_rs, self.eccentricity,
                       self.inclination, self.periastron, self.eclipse_mid_time, time, precision=precision)

    def eclipse_integrated(self, time, time_format, exp_time, time_stamp, filter_name, wlrange=None,
                           max_sub_exp_time=10, precision=3):

        if time_stamp == 'start':
            time = np.array(time) + 0.5 * exp_time / (60.0 * 60.0 * 24.0)
        elif time_stamp == 'mid':
            time = np.array(time)
        elif time_stamp == 'end':
            time = np.array(time) - 0.5 * exp_time / (60.0 * 60.0 * 24.0)
        else:
            raise PyLCInputError(
                'Not acceptable time stamp {0}. Please choose between "mid", "start", "end".'.format(time_stamp))

        time = self.target.convert_to_bjd_tdb(time, time_format)

        return eclipse_integrated(self.fp_over_fs(filter_name, wlrange=wlrange), self.rp_over_rs,
                                  self.period, self.sma_over_rs, self.eccentricity,
                                  self.inclination, self.periastron, self.eclipse_mid_time, time, exp_time,
                                  max_sub_exp_time=max_sub_exp_time, precision=precision)

    def eclipse_duration(self):
        return eclipse_duration(self.rp_over_rs, self.period, self.sma_over_rs, self.eccentricity,
                                self.inclination, self.periastron)

    def eclipse_depth(self, filter_name, wlrange=None, precision=6):
        return eclipse_depth(self.fp_over_fs(filter_name, wlrange), self.rp_over_rs, self.period,
                             self.sma_over_rs, self.eccentricity, self.inclination, self.periastron,
                             precision=precision)

    # data fitting

    def add_observation(self, time, time_format, exp_time, time_stamp, flux, flux_unc, flux_format, filter_name,
                        wlrange=None, observatory_latitude=None, observatory_longitude=None, auxiliary_data=None):

        _ = get_filter(filter_name)

        original_data = {
            'time': time,
            'time_stamp': time_stamp,
            'time_format': time_format,
            'flux': flux,
            'flux_format': flux_format,
            'flux_unc': flux_unc,
            'exp_time': exp_time,
            'filter_name': filter_name,
            'wlrange': wlrange,
        }

        if observatory_latitude and observatory_longitude:
            observatory_latitude = _reformat_or_request_angle(observatory_latitude)
            observatory_longitude = _reformat_or_request_angle(observatory_longitude)
            original_data['observatory_latitude'] = observatory_latitude.deg_coord()
            original_data['observatory_longitude'] = observatory_longitude.deg()
            observatory = Observatory(observatory_latitude, observatory_longitude)
        else:
            observatory = None
            original_data['observatory_latitude'] = None
            original_data['observatory_longitude'] = None

        if time_stamp == 'start':
            time = np.array(time) + 0.5 * exp_time / (60.0 * 60.0 * 24.0)
        elif time_stamp == 'mid':
            time = np.array(time)
        elif time_stamp == 'end':
            time = np.array(time) - 0.5 * exp_time / (60.0 * 60.0 * 24.0)
        else:
            raise PyLCInputError(
                'Not acceptable time stamp {0}. Please choose between "mid", "start", "end".'.format(time_stamp))

        date = self.target.convert_to_jd_utc(time[0], time_format)

        time = np.array([self.target.convert_to_bjd_tdb(ff, time_format) for ff in time])

        if flux_format == 'mag':
            flux_unc = np.abs(np.array(flux_unc) * (-0.921034 * np.exp(0.921034 * flux[0] - 0.921034 * np.array(flux))))
            flux = 10 ** ((flux[0] - np.array(flux)) / 2.5)
        elif flux_format == 'flux':
            flux = np.array(flux)
            flux_unc = np.array(flux_unc)
        else:
            raise PyLCInputError('Not acceptable flux format {0}. Please choose between "flux" and '
                                 '"mag".'.format(flux_format))

        obs_id = 0
        while obs_id in self.observations:
            obs_id += 1

        check_transit = (np.mean(time) - self.mid_time) / self.period
        check_transit = abs(check_transit - int(check_transit))

        check_eclipse = (np.mean(time) - self.eclipse_mid_time) / self.period
        check_eclipse = abs(check_eclipse - int(check_eclipse))

        if check_transit < check_eclipse:
            observation_type = 'transit'
            epoch = int(round((np.mean(time) - self.mid_time) / self.period, 0))
        else:
            observation_type = 'eclipse'
            epoch = int(round((np.mean(time) - self.eclipse_mid_time) / self.period, 0))

        if observatory:
            airmass = []
            for i in time:
                airmass.append(observatory.airmass(self.target, JD(self.target.convert_to_jd_utc(i, time_format))))
        else:
            airmass = np.zeros_like(time)

        dtime = time-time[0]

        self.observations[obs_id] = {
            'target': self.name,
            'time': time,
            'flux': flux / np.median(flux),
            'flux_unc': flux_unc / np.median(flux),
            'exp_time': exp_time,
            'filter_name': filter_name,
            'wlrange': wlrange,
            'epoch': epoch,
            'date': date,
            'observation_type': observation_type,
            'original_data': original_data,
            'auxiliary_data': {
                'dtime': np.array(dtime),
                'dtimesqr': np.array(dtime)**2,
                'airmass': np.array(airmass),
            }
        }

        if auxiliary_data is not None:
            for auxiliary_timeseries in auxiliary_data:
                self.observations[obs_id]['auxiliary_data'][auxiliary_timeseries] = auxiliary_data[auxiliary_timeseries]

    def add_observation_from_dict(self, dictionary):
        if 'wlrange' not in dictionary:
            dictionary['wlrange'] = None
        if 'observatory_latitude' not in dictionary:
            dictionary['observatory_latitude'] = None
        if 'observatory_longitude' not in dictionary:
            dictionary['observatory_longitude'] = None
        if 'auxiliary_data' not in dictionary:
            dictionary['auxiliary_data'] = None

        self.add_observation(
            time=dictionary['time'],
            time_format=dictionary['time_format'],
            exp_time=dictionary['exp_time'],
            time_stamp=dictionary['time_stamp'],
            flux=dictionary['flux'],
            flux_unc=dictionary['flux_unc'],
            flux_format=dictionary['flux_format'],
            filter_name=dictionary['filter_name'],
            wlrange=dictionary['wlrange'],
            observatory_latitude=dictionary['observatory_latitude'],
            observatory_longitude=dictionary['observatory_longitude'],
            auxiliary_data=dictionary['auxiliary_data'])

    def transit_fitting(self, output_folder=None,

                        method='claret', stellar_model='Phoenix_2018',

                        max_sub_exp_time=10, precision=3,

                        detrending='time', detrending_order=2,

                        trend_function=None, trend_priors=None,

                        iterations=None, walkers=None, burn_in=None,

                        fit_ldc1=False, fit_ldc2=False, fit_ldc3=False, fit_ldc4=False,

                        fit_rp_over_rs=True, fit_individual_rp_over_rs=False,

                        fit_sma_over_rs=False, fit_inclination=False,

                        fit_mid_time=True, fit_period=False,
                        fit_individual_times=True,

                        fit_ldc_limits=[-1.0, 1.0],
                        fit_rp_over_rs_limits=[0.25, 4.0],
                        fit_sma_over_rs_limits=[0.25, 4.0],
                        fit_inclination_limits=[70.0, 90.0],
                        fit_mid_time_limits=[-0.2, 0.2],
                        fit_period_limits=[0.99, 1.01],

                        counter='Transit fitting',
                        optimise_initial_parameters=True,
                        optimise_initial_parameters_trials=3,
                        scale_uncertainties=False,
                        filter_outliers=False,
                        optimiser='emcee',
                        window_counter=False,
                        strech_prior=0.2,
                        return_traces=True,
                        ):

        parameters_map = [[] for observation in self.observations]

        names = []
        print_names = []
        limits1 = []
        limits2 = []
        initial = []

        if trend_function is None:

            if detrending == 'time':
                if detrending_order == 2:
                    def trend_function(auxilary_data, c1, c2):
                        return 1 + c1 * auxilary_data['dtime'] + c2 * auxilary_data['dtimesqr']
                elif detrending_order == 1:
                    def trend_function(auxilary_data, c1):
                        return 1 + c1 * auxilary_data['dtime']
                else:
                    def trend_function(auxilary_data):
                        return np.ones_like(auxilary_data['dtime'])
            else:
                if detrending_order == 2:
                    def trend_function(auxilary_data, c1, c2):
                        trend_x = auxilary_data['airmass']
                        return 1 + c1 * trend_x + c2 * trend_x * trend_x
                elif detrending_order == 1:
                    def trend_function(auxilary_data, c1):
                        trend_x = auxilary_data['airmass']
                        return 1 + c1 * trend_x
                else:
                    def trend_function(auxilary_data):
                        return np.ones_like(auxilary_data['dtime'])

        trend = Trend(function=trend_function, trend_priors=trend_priors)

        # de-trending parameters

        for observation_num, observation in enumerate(self.observations):

            trend.adjust_normalisation_factor(self.observations[observation]['flux'])

            for coefficient in range(trend.coefficients):

                if len(self.observations) == 1:
                    names.append('{0}'.format(trend.names[coefficient]))
                    print_names.append('{0}'.format(trend.print_names[coefficient]))
                else:
                    names.append('{0}_{1}'.format(trend.names[coefficient], observation_num + 1))
                    print_names.append('{0}_{1}'.format(trend.print_names[coefficient], observation_num + 1))

                initial.append(trend.initial[coefficient])
                limits1.append(trend.limits1[coefficient])
                limits2.append(trend.limits2[coefficient])

                parameters_map[observation_num].append(len(names) - 1)

        # limb-darkening and rp_over_rs parameters

        unique_filters = []
        for observation in self.observations:
            if self.observations[observation]['wlrange'] is not None:
                ff = '{0}>>>{1}>>>{2}'.format(self.observations[observation]['filter_name'],
                                          self.observations[observation]['wlrange'][0],
                                          self.observations[observation]['wlrange'][1])
            else:
                ff = self.observations[observation]['filter_name']
            if ff not in unique_filters:
                unique_filters.append(ff)

        if len(unique_filters) == 1:
            fit_individual_rp_over_rs = False

        for phot_filter in unique_filters:

            if len(phot_filter.split('>>>')) == 1:
                ldc1, ldc2, ldc3, ldc4 = self.ldcs(phot_filter, method=method, stellar_model=stellar_model)
            else:
                ldc1, ldc2, ldc3, ldc4 = self.ldcs(phot_filter.split('>>>')[0],
                                                   wlrange=[float(phot_filter.split('>>>')[1]),
                                                            float(phot_filter.split('>>>')[2])],
                                                   method=method, stellar_model=stellar_model)

            rp_over_rs = self.rp_over_rs

            names.append('LDC1_{0}'.format(phot_filter.replace('>>>', '_')))
            print_names.append('LDC1_\mathrm{0}'.format(in_brackets(phot_filter.replace('>>>', '-'))))
            initial.append(ldc1)
            if fit_ldc1:
                limits1.append(ldc1 + fit_ldc_limits[0])
                limits2.append(ldc1 + fit_ldc_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                if self.observations[observation]['wlrange'] is not None:
                    ff = '{0}>>>{1}>>>{2}'.format(self.observations[observation]['filter_name'],
                                                  self.observations[observation]['wlrange'][0],
                                                  self.observations[observation]['wlrange'][1])
                else:
                    ff = self.observations[observation]['filter_name']
                if ff == phot_filter:
                    parameters_map[observation_num].append(len(names) - 1)

            names.append('LDC2_{0}'.format(phot_filter.replace('>>>', '_')))
            print_names.append('LDC2_\mathrm{0}'.format(in_brackets(phot_filter.replace('>>>', '-'))))
            initial.append(ldc2)
            if fit_ldc2:
                limits1.append(ldc2 + fit_ldc_limits[0])
                limits2.append(ldc2 + fit_ldc_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                if self.observations[observation]['wlrange'] is not None:
                    ff = '{0}>>>{1}>>>{2}'.format(self.observations[observation]['filter_name'],
                                                  self.observations[observation]['wlrange'][0],
                                                  self.observations[observation]['wlrange'][1])
                else:
                    ff = self.observations[observation]['filter_name']
                if ff == phot_filter:
                    parameters_map[observation_num].append(len(names) - 1)

            names.append('LDC3_{0}'.format(phot_filter.replace('>>>', '_')))
            print_names.append('LDC3_\mathrm{0}'.format(in_brackets(phot_filter.replace('>>>', '-'))))
            initial.append(ldc3)
            if fit_ldc3:
                limits1.append(ldc3 + fit_ldc_limits[0])
                limits2.append(ldc3 + fit_ldc_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                if self.observations[observation]['wlrange'] is not None:
                    ff = '{0}>>>{1}>>>{2}'.format(self.observations[observation]['filter_name'],
                                                  self.observations[observation]['wlrange'][0],
                                                  self.observations[observation]['wlrange'][1])
                else:
                    ff = self.observations[observation]['filter_name']
                if ff == phot_filter:
                    parameters_map[observation_num].append(len(names) - 1)

            names.append('LDC4_{0}'.format(phot_filter.replace('>>>', '_')))
            print_names.append('LDC4_\mathrm{0}'.format(in_brackets(phot_filter.replace('>>>', '-'))))
            initial.append(ldc4)
            if fit_ldc4:
                limits1.append(ldc4 + fit_ldc_limits[0])
                limits2.append(ldc4 + fit_ldc_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                if self.observations[observation]['wlrange'] is not None:
                    ff = '{0}>>>{1}>>>{2}'.format(self.observations[observation]['filter_name'],
                                                  self.observations[observation]['wlrange'][0],
                                                  self.observations[observation]['wlrange'][1])
                else:
                    ff = self.observations[observation]['filter_name']
                if ff == phot_filter:
                    parameters_map[observation_num].append(len(names) - 1)

            if fit_individual_rp_over_rs:

                names.append('rp_{0}'.format(phot_filter.replace('>>>', '_')))
                print_names.append('(R_\mathrm{p}/R_*)_\mathrm{' + phot_filter.replace('>>>', '-') + '}')
                initial.append(rp_over_rs)
                limits1.append(rp_over_rs * fit_rp_over_rs_limits[0])
                limits2.append(rp_over_rs * fit_rp_over_rs_limits[1])

                for observation_num, observation in enumerate(self.observations):
                    if self.observations[observation]['wlrange'] is not None:
                        ff = '{0}>>>{1}>>>{2}'.format(self.observations[observation]['filter_name'],
                                                      self.observations[observation]['wlrange'][0],
                                                      self.observations[observation]['wlrange'][1])
                    else:
                        ff = self.observations[observation]['filter_name']
                    if ff == phot_filter:
                        parameters_map[observation_num].append(len(names) - 1)

        if not fit_individual_rp_over_rs:

            names.append('rp')
            print_names.append('R_\mathrm{p}/R_*')
            initial.append(self.rp_over_rs)
            if fit_rp_over_rs:
                limits1.append(self.rp_over_rs * fit_rp_over_rs_limits[0])
                limits2.append(self.rp_over_rs * fit_rp_over_rs_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                parameters_map[observation_num].append(len(names) - 1)

        # orbital parameters

        names.append('P')
        print_names.append('P')
        initial.append(self.period)
        if fit_period:
            limits1.append(self.period * fit_period_limits[0])
            limits2.append(self.period * fit_period_limits[1])
        else:
            limits1.append(np.nan)
            limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        names.append('a')
        print_names.append('a/R_*')
        initial.append(self.sma_over_rs)
        if fit_sma_over_rs:
            limits1.append(self.sma_over_rs * fit_sma_over_rs_limits[0])
            limits2.append(self.sma_over_rs * fit_sma_over_rs_limits[1])
        else:
            limits1.append(np.nan)
            limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        names.append('e')
        print_names.append('e')
        initial.append(self.eccentricity)
        limits1.append(np.nan)
        limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        names.append('i')
        print_names.append('i')
        initial.append(self.inclination)
        if fit_inclination:
            limits1.append(fit_inclination_limits[0])
            limits2.append(fit_inclination_limits[1])
        else:
            limits1.append(np.nan)
            limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        names.append('w')
        print_names.append('\omega')
        initial.append(self.periastron)
        limits1.append(np.nan)
        limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        # time parameters

        test_epochs = []
        test_epochs_weights = []

        for observation in self.observations:
            test_epochs.append((self.observations[observation]['time'] - self.mid_time) / self.period)
            norm_errors = self.observations[observation]['flux_unc'] / self.observations[observation]['flux']
            test_epochs_weights.append(1 / (norm_errors * norm_errors))

        test_epochs = np.concatenate(test_epochs)
        test_epochs_weights = np.concatenate(test_epochs_weights)

        new_epoch = np.round(np.sum(test_epochs * test_epochs_weights) / np.sum(test_epochs_weights), 0)
        new_mid_time = self.mid_time + new_epoch * self.period

        for observation in self.observations:
            self.observations[observation]['epoch'] = int(round((np.mean(self.observations[observation]['time']) -new_mid_time) / self.period, 0))

        unique_epochs = []
        for observation in self.observations:
            if self.observations[observation]['epoch'] not in unique_epochs:
                unique_epochs.append(self.observations[observation]['epoch'])

        if len(unique_epochs) == 1:
            fit_individual_times = False

        if not fit_individual_times:

            names.append('T_mid')
            print_names.append('T_\mathrm{mid}')
            initial.append(new_mid_time)
            if fit_mid_time:
                limits1.append(new_mid_time + fit_mid_time_limits[0])
                limits2.append(new_mid_time + fit_mid_time_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                parameters_map[observation_num].append(len(names) - 1)

        else:

            if fit_period:
                raise PyLCInputError('Period and individual mid times cannot be fitted simultaneously.')

            for epoch in unique_epochs:

                names.append('T_mid_{0}'.format(epoch))
                print_names.append('T_\mathrm{mid_{' + str(epoch) + '}}')
                initial.append(new_mid_time + epoch * self.period)
                limits1.append(new_mid_time + epoch * self.period + fit_mid_time_limits[0])
                limits2.append(new_mid_time + epoch * self.period + fit_mid_time_limits[1])

                for observation_num, observation in enumerate(self.observations):
                    if self.observations[observation]['epoch'] == epoch:
                        parameters_map[observation_num].append(len(names) - 1)

        model_ids = np.array([])
        model_flux = np.array([])
        model_flux_unc = np.array([])
        for observation in self.observations:
            model_ids = np.append(model_ids, np.int_(10000*observation + np.arange(0, len(self.observations[observation]['time']))))
            model_flux = np.append(model_flux, self.observations[observation]['flux'])
            model_flux_unc = np.append(model_flux_unc, self.observations[observation]['flux_unc'])

        parameters_map = np.array(parameters_map)

        def detrend_model(model_ids, *model_variables):

            model = np.array([])

            for observation in self.observations:
                coefficients = np.array(model_variables)[parameters_map[observation]][:trend.coefficients]
                ids = np.int_(model_ids[np.where(np.int_(model_ids/10000.0) == observation)] - observation * 10000.0)
                model = np.append(model, trend.evaluate(self.observations[observation]['auxiliary_data'], *coefficients)[ids])

            return model

        def full_model(model_ids, *model_variables):

            model_variables = np.array(model_variables)

            model = np.array([])

            for observation in self.observations:

                coefficients = np.array(model_variables)[parameters_map[observation]][:trend.coefficients]
                ldc1, ldc2, ldc3, ldc4, r, p, a, e, i, w, mt = np.array(model_variables)[parameters_map[observation]][trend.coefficients:]
                ids = np.int_(model_ids[np.where(np.int_(model_ids/10000.0) == observation)] - observation * 10000.0)

                trend_model = trend.evaluate(
                    self.observations[observation]['auxiliary_data'],
                    *coefficients)

                transit_model = transit_integrated(
                    [ldc1, ldc2, ldc3, ldc4], r, p, a, e, i, w, mt,
                    time_array=self.observations[observation]['time'],
                    exp_time=self.observations[observation]['exp_time'],
                    max_sub_exp_time=max_sub_exp_time,
                    method=method,
                    precision=precision
                )

                model = np.append(model, (trend_model * transit_model)[ids])

            return model

        fitting = Fitting(model_ids, model_flux, model_flux_unc,
                          full_model, initial, limits1, limits2,
                          data_x_name='time', data_y_name='flux',
                          data_x_print_name='t_{BJD_TDB}',  data_y_print_name='Relative Flux',
                          parameters_names=names, parameters_print_names=print_names,
                          walkers=walkers, iterations=iterations, burn_in=burn_in,
                          counter=counter, window_counter=window_counter,
                          optimise_initial_parameters=optimise_initial_parameters,
                          optimise_initial_parameters_trials=optimise_initial_parameters_trials,
                          scale_uncertainties=scale_uncertainties,
                          filter_outliers=filter_outliers,
                          optimiser=optimiser,
                          strech_prior=strech_prior,
                          )

        fitting.run()

        if fitting.mcmc_run_complete:

            fitting.results['settings']['max_sub_exp_time'] = max_sub_exp_time
            fitting.results['settings']['precision'] = precision
            fitting.results['settings']['detrending_order'] = detrending_order
            fitting.results['settings']['iterations'] = iterations
            fitting.results['settings']['walkers'] = walkers
            fitting.results['settings']['burn_in'] = burn_in
            fitting.results['settings']['fit_ldc1'] = fit_ldc1
            fitting.results['settings']['fit_ldc2'] = fit_ldc2
            fitting.results['settings']['fit_ldc3'] = fit_ldc3
            fitting.results['settings']['fit_ldc4'] = fit_ldc4
            fitting.results['settings']['fit_rp_over_rs'] = fit_rp_over_rs
            fitting.results['settings']['fit_individual_rp_over_rs'] = fit_individual_rp_over_rs
            fitting.results['settings']['fit_sma_over_rs'] = fit_sma_over_rs
            fitting.results['settings']['fit_inclination'] = fit_inclination
            fitting.results['settings']['fit_mid_time'] = fit_mid_time
            fitting.results['settings']['fit_period'] = fit_period
            fitting.results['settings']['fit_individual_times'] = fit_individual_times
            fitting.results['settings']['fit_ldc_limits'] = fit_ldc_limits
            fitting.results['settings']['fit_rp_over_rs_limits'] = fit_rp_over_rs_limits
            fitting.results['settings']['fit_sma_over_rs_limits'] = fit_sma_over_rs_limits
            fitting.results['settings']['fit_inclination_limits'] = fit_inclination_limits
            fitting.results['settings']['fit_mid_time_limits'] = fit_mid_time_limits
            fitting.results['settings']['fit_period_limits'] = fit_period_limits
            fitting.results['settings']['filter_map'] = self.filters

            original_model_time = np.array([])
            for observation_num, observation in enumerate(self.observations):
                original_model_time = np.append(original_model_time, self.observations[observation]['time'])

            model_ids = fitting.results['input_series']['time']
            model_time = np.array([])
            for observation in self.observations:
                ids = np.int_(model_ids[np.where(np.int_(model_ids/10000.0) == observation)] - observation * 10000.0)
                model_time = np.append(model_time, self.observations[observation]['time'][ids])

            trend = detrend_model(model_ids, *fitting.results['parameters_final'])

            fitting.results['settings']['time'] = original_model_time
            fitting.results['input_series']['time'] = model_time
            fitting.results['output_series']['trend'] = trend

            fitting.results['detrended_input_series'] = {
                'time': model_time,
                'flux': fitting.results['input_series']['flux'] / trend,
                'flux_unc': fitting.results['input_series']['flux_unc'] / trend
            }

            fitting.results['detrended_output_series'] = {
                'model': fitting.results['output_series']['model'] / trend,
                'residuals': fitting.results['output_series']['residuals'] / trend
            }

            fitting.results['detrended_statistics'] = {
                'res_mean': np.mean(fitting.results['detrended_output_series']['residuals']),
                'res_std': np.std(fitting.results['detrended_output_series']['residuals']),
                'res_rms': np.sqrt(np.mean(fitting.results['detrended_output_series']['residuals'] ** 2)),
            }

            for observation in self.observations:

                ids = np.int_(model_ids[np.where(np.int_(model_ids/10000.0) == observation)] - observation * 10000.0)

                self.observations[observation]['time'] = self.observations[observation]['time'][ids]
                self.observations[observation]['flux'] = self.observations[observation]['flux'][ids]
                self.observations[observation]['flux_unc'] = self.observations[observation]['flux_unc'][ids]
                self.observations[observation]['auxiliary_data'] = {ff:self.observations[observation]['auxiliary_data'][ff][ids] for ff in self.observations[observation]['auxiliary_data']}

                id_series = np.where(np.int_(model_ids/10000.0) == observation)

                self.observations[observation]['model'] = fitting.results['output_series']['model'][id_series]
                self.observations[observation]['residuals'] = fitting.results['output_series']['residuals'][id_series]
                self.observations[observation]['detrended_flux'] = fitting.results['detrended_input_series']['flux'][id_series]
                self.observations[observation]['detrended_flux_unc'] = fitting.results['detrended_input_series']['flux_unc'][id_series]
                self.observations[observation]['detrended_model'] = fitting.results['detrended_output_series']['model'][id_series]
                self.observations[observation]['detrended_residuals'] = fitting.results['detrended_output_series']['residuals'][id_series]

                self.observations[observation]['res_mean'] = np.mean(self.observations[observation]['residuals'])
                self.observations[observation]['res_std'] = np.std(self.observations[observation]['residuals'])
                self.observations[observation]['res_rms'] = np.sqrt(np.mean(self.observations[observation]['residuals'] ** 2))
                self.observations[observation]['res_chi_sqr'] = np.sum(
                    (self.observations[observation]['residuals']/self.observations[observation]['flux_unc']) ** 2)

                self.observations[observation]['detrended_res_mean'] = np.mean(self.observations[observation]['detrended_residuals'])
                self.observations[observation]['detrended_res_std'] = np.std(self.observations[observation]['detrended_residuals'])
                self.observations[observation]['detrended_res_rms'] = np.sqrt(np.mean(self.observations[observation]['detrended_residuals'] ** 2))
                self.observations[observation]['detrended_res_chi_sqr'] = np.sum(
                    (self.observations[observation]['detrended_residuals']/self.observations[observation]['detrended_flux_unc']) ** 2)

            fitting.results['observations'] = self.observations

            results_copy = copy_dict(fitting.results)
            if not return_traces:
                for parameter in results_copy['parameters']:
                    results_copy['parameters'][parameter]['trace'] = None

            if output_folder:

                if not os.path.isdir(output_folder):
                    os.mkdir(output_folder)
                else:
                    shutil.rmtree(output_folder)
                    os.mkdir(output_folder)

                save_dict(results_copy, os.path.join(output_folder, 'results.pickle'))

                fitting.save_results(os.path.join(output_folder, 'results.txt'))

                fitting.plot_corner(os.path.join(output_folder, 'correlations.pdf'))

                fitting.plot_traces(os.path.join(output_folder, 'traces.pdf'))

                plot_transit_fitting_models(fitting.results, os.path.join(output_folder, 'lightcurves.pdf'))

                for observation in self.observations:

                    w = open(os.path.join(output_folder, 'diagnostics_dataset_{0}.txt'.format(observation + 1)), 'w')
                    w.write('\n#Residuals:\n')
                    w.write('#Mean: {0}\n'.format(self.observations[observation]['res_mean']))
                    w.write('#STD: {0}\n'.format(self.observations[observation]['res_std']))
                    w.write('#RMS: {0}\n'.format(self.observations[observation]['res_rms']))
                    w.write('#Chi squared: {0}\n'.format(self.observations[observation]['res_chi_sqr']))

                    w.write('\n\n#Detrended Residuals:\n')
                    w.write('#Mean: {0}\n'.format(self.observations[observation]['detrended_res_mean']))
                    w.write('#STD: {0}\n'.format(self.observations[observation]['detrended_res_std']))
                    w.write('#RMS: {0}\n'.format(self.observations[observation]['detrended_res_rms']))
                    w.write('#Chi squared: {0}\n'.format(self.observations[observation]['detrended_res_chi_sqr']))

                    w.close()

            fitting.results = results_copy

        return fitting

    def eclipse_fitting(self, output_folder=None,

                            max_sub_exp_time=10, precision=3,

                            detrending='time', detrending_order=2,

                            trend_function=None, trend_priors=None,

                            iterations=None, walkers=None, burn_in=None,

                            fit_fp_over_fs=True, fit_individual_fp_over_fs=False,
                            fit_rp_over_rs=False, fit_individual_rp_over_rs=False,

                            fit_sma_over_rs=False, fit_inclination=False,

                            fit_mid_time=True, fit_period=False,
                            fit_individual_times=True,

                            fit_rp_over_rs_limits=[0.25, 4.0],
                            fit_fp_over_fs_limits=[0.001, 1000.0],
                            fit_sma_over_rs_limits=[0.25, 4.0],
                            fit_inclination_limits=[70.0, 90.0],
                            fit_mid_time_limits=[-0.2, 0.2],
                            fit_period_limits=[0.99, 1.01],

                            counter='Eclipse fitting',
                            optimise_initial_parameters=True,
                            optimise_initial_parameters_trials=4,
                            scale_uncertainties=False,
                            filter_outliers=False,
                            optimiser='emcee',
                            window_counter=False,
                            strech_prior=0.2,
                            return_traces=True,
                            ):


        parameters_map = [[] for observation in self.observations]

        names = []
        print_names = []
        limits1 = []
        limits2 = []
        initial = []

        if trend_function is None:

            if detrending=='time':
                if detrending_order == 2:
                    def trend_function(auxilary_data, c1, c2):
                        return 1 + c1 * auxilary_data['dtime'] + c2 * auxilary_data['dtimesqr']
                elif detrending_order == 1:
                    def trend_function(auxilary_data, c1):
                        return 1 + c1 * auxilary_data['dtime']
                else:
                    def trend_function(auxilary_data):
                        return np.ones_like(auxilary_data['dtime'])
            else:
                if detrending_order == 2:
                    def trend_function(auxilary_data, c1, c2):
                        trend_x = auxilary_data['airmass']
                        return 1 + c1 * trend_x + c2 * trend_x * trend_x
                elif detrending_order == 1:
                    def trend_function(auxilary_data, c1):
                        trend_x = auxilary_data['airmass']
                        return 1 + c1 * trend_x
                else:
                    def trend_function(auxilary_data):
                        return np.ones_like(auxilary_data['dtime'])

        trend = Trend(function=trend_function, trend_priors=trend_priors)

        # de-trending parameters

        for observation_num, observation in enumerate(self.observations):

            trend.adjust_normalisation_factor(self.observations[observation]['flux'])

            for coefficient in range(trend.coefficients):

                if len(self.observations) == 1:
                    names.append('{0}'.format(trend.names[coefficient]))
                    print_names.append('{0}'.format(trend.print_names[coefficient]))
                else:
                    names.append('{0}_{1}'.format(trend.names[coefficient], observation_num + 1))
                    print_names.append('{0}_{1}'.format(trend.print_names[coefficient], observation_num + 1))

                initial.append(trend.initial[coefficient])
                limits1.append(trend.limits1[coefficient])
                limits2.append(trend.limits2[coefficient])

                parameters_map[observation_num].append(len(names) - 1)

        # limb-darkening and rp_over_rs parameters

        unique_filters = []
        for observation in self.observations:
            if self.observations[observation]['wlrange'] is not None:
                ff = '{0}>>>{1}>>>{2}'.format(self.observations[observation]['filter_name'],
                                              self.observations[observation]['wlrange'][0],
                                              self.observations[observation]['wlrange'][1])
            else:
                ff = self.observations[observation]['filter_name']
            if ff not in unique_filters:
                unique_filters.append(ff)

        if len(unique_filters) == 1:
            fit_individual_rp_over_rs = False
            fit_individual_fp_over_fs = False

        if fit_individual_fp_over_fs:

            for phot_filter in unique_filters:

                if len(phot_filter.split('>>>')) == 1:
                    fp_over_fs = self.fp_over_fs(phot_filter)
                else:
                    fp_over_fs= self.fp_over_fs(phot_filter.split('>>>')[0],
                                                wlrange=[float(phot_filter.split('>>>')[1]),
                                                         float(phot_filter.split('>>>')[2])])


                names.append('fp_{0}'.format(phot_filter.replace('>>>', '_')))
                print_names.append('(F_\mathrm{p}/F_*)_\mathrm{' + phot_filter.replace('>>>', '-') + '}')
                initial.append(fp_over_fs)
                limits1.append(fp_over_fs * fit_fp_over_fs_limits[0])
                limits2.append(fp_over_fs * fit_fp_over_fs_limits[1])

                for observation_num, observation in enumerate(self.observations):
                    if self.observations[observation]['wlrange'] is not None:
                        ff = '{0}>>>{1}>>>{2}'.format(self.observations[observation]['filter_name'],
                                                      self.observations[observation]['wlrange'][0],
                                                      self.observations[observation]['wlrange'][1])
                    else:
                        ff = self.observations[observation]['filter_name']
                    if ff == phot_filter:
                        parameters_map[observation_num].append(len(names) - 1)

        else:
            fp_over_fs = 0.000001

            names.append('fp')
            print_names.append('F_\mathrm{p}/F_*')
            initial.append(fp_over_fs)
            if fit_fp_over_fs:
                limits1.append(fp_over_fs * fit_fp_over_fs_limits[0])
                limits2.append(fp_over_fs * fit_fp_over_fs_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                parameters_map[observation_num].append(len(names) - 1)

        if fit_individual_rp_over_rs:

            for phot_filter in unique_filters:

                rp_over_rs = self.rp_over_rs

                names.append('rp_{0}'.format(phot_filter.replace('>>>', '_')))
                print_names.append('(R_\mathrm{p}/R_*)_\mathrm{' + phot_filter.replace('>>>', '-') + '}')
                initial.append(rp_over_rs)
                limits1.append(rp_over_rs * fit_rp_over_rs_limits[0])
                limits2.append(rp_over_rs * fit_rp_over_rs_limits[1])

                for observation_num, observation in enumerate(self.observations):
                    if self.observations[observation]['wlrange'] is not None:
                        ff = '{0}>>>{1}>>>{2}'.format(self.observations[observation]['filter_name'],
                                                      self.observations[observation]['wlrange'][0],
                                                      self.observations[observation]['wlrange'][1])
                    else:
                        ff = self.observations[observation]['filter_name']
                    if ff == phot_filter:
                        parameters_map[observation_num].append(len(names) - 1)

        else:

            rp_over_rs = self.rp_over_rs

            names.append('rp')
            print_names.append('R_\mathrm{p}/R_*')
            initial.append(rp_over_rs)
            if fit_rp_over_rs:
                limits1.append(rp_over_rs * fit_rp_over_rs_limits[0])
                limits2.append(rp_over_rs * fit_rp_over_rs_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                parameters_map[observation_num].append(len(names) - 1)

        # orbital parameters

        names.append('P')
        print_names.append('P')
        initial.append(self.period)
        if fit_period:
            limits1.append(self.period * fit_period_limits[0])
            limits2.append(self.period * fit_period_limits[1])
        else:
            limits1.append(np.nan)
            limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        names.append('a')
        print_names.append('a/R_*')
        initial.append(self.sma_over_rs)
        if fit_sma_over_rs:
            limits1.append(self.sma_over_rs * fit_sma_over_rs_limits[0])
            limits2.append(self.sma_over_rs * fit_sma_over_rs_limits[1])
        else:
            limits1.append(np.nan)
            limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        names.append('e')
        print_names.append('e')
        initial.append(self.eccentricity)
        limits1.append(np.nan)
        limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        names.append('i')
        print_names.append('i')
        initial.append(self.inclination)
        if fit_inclination:
            limits1.append(fit_inclination_limits[0])
            limits2.append(fit_inclination_limits[1])
        else:
            limits1.append(np.nan)
            limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        names.append('w')
        print_names.append('\omega')
        initial.append(self.periastron)
        limits1.append(np.nan)
        limits2.append(np.nan)

        for observation_num, observation in enumerate(self.observations):
            parameters_map[observation_num].append(len(names) - 1)

        # time parameters

        test_epochs = []
        test_epochs_weights = []

        for observation in self.observations:
            test_epochs.append((self.observations[observation]['time'] - self.eclipse_mid_time) / self.period)
            norm_errors = self.observations[observation]['flux_unc'] / self.observations[observation]['flux']
            test_epochs_weights.append(1 / (norm_errors * norm_errors))

        test_epochs = np.concatenate(test_epochs)
        test_epochs_weights = np.concatenate(test_epochs_weights)

        new_epoch = np.round(np.sum(test_epochs * test_epochs_weights) / np.sum(test_epochs_weights), 0)
        new_mid_time = self.eclipse_mid_time + new_epoch * self.period

        for observation in self.observations:
            self.observations[observation]['epoch'] = int(round((np.mean(self.observations[observation]['time']) -new_mid_time) / self.period, 0))

        unique_epochs = []
        for observation in self.observations:
            if self.observations[observation]['epoch'] not in unique_epochs:
                unique_epochs.append(self.observations[observation]['epoch'])

        if len(unique_epochs) == 1:
            fit_individual_times = False

        if not fit_individual_times:

            names.append('T_mid')
            print_names.append('T_\mathrm{mid}')
            initial.append(new_mid_time)
            if fit_mid_time:
                limits1.append(new_mid_time + fit_mid_time_limits[0])
                limits2.append(new_mid_time + fit_mid_time_limits[1])
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            for observation_num, observation in enumerate(self.observations):
                parameters_map[observation_num].append(len(names) - 1)

        else:

            if fit_period:
                raise PyLCInputError('Period and individual mid times cannot be fitted simultaneously.')

            for epoch in unique_epochs:

                names.append('T_mid_{0}'.format(epoch))
                print_names.append('T_\mathrm{mid_{' + str(epoch) + '}}')
                initial.append(new_mid_time + epoch * self.period)
                limits1.append(new_mid_time + epoch * self.period + fit_mid_time_limits[0])
                limits2.append(new_mid_time + epoch * self.period + fit_mid_time_limits[1])

                for observation_num, observation in enumerate(self.observations):
                    if self.observations[observation]['epoch'] == epoch:
                        parameters_map[observation_num].append(len(names) - 1)

        model_ids = np.array([])
        model_flux = np.array([])
        model_flux_unc = np.array([])
        for observation in self.observations:
            model_ids = np.append(model_ids, np.int_(10000*observation + np.arange(0, len(self.observations[observation]['time']))))
            model_flux = np.append(model_flux, self.observations[observation]['flux'])
            model_flux_unc = np.append(model_flux_unc, self.observations[observation]['flux_unc'])

        parameters_map = np.array(parameters_map)

        def detrend_model(model_ids, *model_variables):

            model = np.array([])

            for observation in self.observations:
                coefficients = np.array(model_variables)[parameters_map[observation]][:trend.coefficients]
                ids = np.int_(model_ids[np.where(np.int_(model_ids/10000.0) == observation)] - observation * 10000.0)
                model = np.append(model, trend.evaluate(self.observations[observation]['auxiliary_data'], *coefficients)[ids])

            return model

        def full_model(model_ids, *model_variables):

            model_variables = np.array(model_variables)

            model = np.array([])

            for observation in self.observations:

                coefficients = np.array(model_variables)[parameters_map[observation]][:trend.coefficients]
                f, r, p, a, e, i, w, mt = np.array(model_variables)[parameters_map[observation]][trend.coefficients:]
                ids = np.int_(model_ids[np.where(np.int_(model_ids/10000.0) == observation)] - observation * 10000.0)

                trend_model = trend.evaluate(
                    self.observations[observation]['auxiliary_data'],
                    *coefficients)

                eclipse_model = eclipse_integrated(
                    f, r, p, a, e, i, w, mt,
                    time_array=self.observations[observation]['time'],
                    exp_time=self.observations[observation]['exp_time'],
                    max_sub_exp_time=max_sub_exp_time,
                    precision=precision
                )

                model = np.append(model, (trend_model * eclipse_model)[ids])

            return model

        fitting = Fitting(model_ids, model_flux, model_flux_unc,
                          full_model, initial, limits1, limits2,
                          data_x_name='time', data_y_name='flux',
                          data_x_print_name='t_{BJD_TDB}',  data_y_print_name='Relative Flux',
                          parameters_names=names, parameters_print_names=print_names,
                          walkers=walkers, iterations=iterations, burn_in=burn_in,
                          counter=counter, window_counter=window_counter,
                          optimise_initial_parameters=optimise_initial_parameters,
                          optimise_initial_parameters_trials=optimise_initial_parameters_trials,
                          scale_uncertainties=scale_uncertainties,
                          filter_outliers=filter_outliers,
                          optimiser=optimiser,
                          strech_prior=strech_prior,
                          )

        fitting.run()

        results_copy = copy_dict(fitting.results)
        for parameter in results_copy['parameters']:
            results_copy['parameters'][parameter]['trace'] = None

        if fitting.mcmc_run_complete:

            fitting.results['settings']['max_sub_exp_time'] = max_sub_exp_time
            fitting.results['settings']['precision'] = precision
            fitting.results['settings']['detrending_order'] = detrending_order
            fitting.results['settings']['iterations'] = iterations
            fitting.results['settings']['walkers'] = walkers
            fitting.results['settings']['burn_in'] = burn_in
            fitting.results['settings']['fit_rp_over_rs'] = fit_rp_over_rs
            fitting.results['settings']['fit_individual_rp_over_rs'] = fit_individual_rp_over_rs
            fitting.results['settings']['fit_fp_over_fs'] = fit_fp_over_fs
            fitting.results['settings']['fit_individual_fp_over_fs'] = fit_individual_fp_over_fs
            fitting.results['settings']['fit_sma_over_rs'] = fit_sma_over_rs
            fitting.results['settings']['fit_inclination'] = fit_inclination
            fitting.results['settings']['fit_mid_time'] = fit_mid_time
            fitting.results['settings']['fit_period'] = fit_period
            fitting.results['settings']['fit_individual_times'] = fit_individual_times
            fitting.results['settings']['fit_fp_over_fs_limits'] = fit_fp_over_fs_limits
            fitting.results['settings']['fit_rp_over_rs_limits'] = fit_rp_over_rs_limits
            fitting.results['settings']['fit_sma_over_rs_limits'] = fit_sma_over_rs_limits
            fitting.results['settings']['fit_inclination_limits'] = fit_inclination_limits
            fitting.results['settings']['fit_mid_time_limits'] = fit_mid_time_limits
            fitting.results['settings']['fit_period_limits'] = fit_period_limits
            fitting.results['settings']['filter_map'] = self.filters

            original_model_time = np.array([])
            for observation_num, observation in enumerate(self.observations):
                original_model_time = np.append(original_model_time, self.observations[observation]['time'])

            model_ids = fitting.results['input_series']['time']
            model_time = np.array([])
            for observation in self.observations:
                ids = np.int_(model_ids[np.where(np.int_(model_ids/10000.0) == observation)] - observation * 10000.0)
                model_time = np.append(model_time, self.observations[observation]['time'][ids])

            trend = detrend_model(model_ids, *fitting.results['parameters_final'])

            fitting.results['settings']['time'] = original_model_time
            fitting.results['input_series']['time'] = model_time
            fitting.results['output_series']['trend'] = trend

            fitting.results['detrended_input_series'] = {
                'time': model_time,
                'flux': fitting.results['input_series']['flux'] / trend,
                'flux_unc': fitting.results['input_series']['flux_unc'] / trend
            }

            fitting.results['detrended_output_series'] = {
                'model': fitting.results['output_series']['model'] / trend,
                'residuals': fitting.results['output_series']['residuals'] / trend
            }

            fitting.results['detrended_statistics'] = {
                'res_mean': np.mean(fitting.results['detrended_output_series']['residuals']),
                'res_std': np.std(fitting.results['detrended_output_series']['residuals']),
                'res_rms': np.sqrt(np.mean(fitting.results['detrended_output_series']['residuals'] ** 2)),
            }

            for observation in self.observations:

                ids = np.int_(model_ids[np.where(np.int_(model_ids/10000.0) == observation)] - observation * 10000.0)

                self.observations[observation]['time'] = self.observations[observation]['time'][ids]
                self.observations[observation]['flux'] = self.observations[observation]['flux'][ids]
                self.observations[observation]['flux_unc'] = self.observations[observation]['flux_unc'][ids]
                self.observations[observation]['auxiliary_data'] = {ff:self.observations[observation]['auxiliary_data'][ff][ids] for ff in self.observations[observation]['auxiliary_data']}

                id_series = np.where(np.int_(model_ids/10000.0) == observation)

                self.observations[observation]['model'] = fitting.results['output_series']['model'][id_series]
                self.observations[observation]['residuals'] = fitting.results['output_series']['residuals'][id_series]
                self.observations[observation]['detrended_flux'] = fitting.results['detrended_input_series']['flux'][id_series]
                self.observations[observation]['detrended_flux_unc'] = fitting.results['detrended_input_series']['flux_unc'][id_series]
                self.observations[observation]['detrended_model'] = fitting.results['detrended_output_series']['model'][id_series]
                self.observations[observation]['detrended_residuals'] = fitting.results['detrended_output_series']['residuals'][id_series]

                self.observations[observation]['res_mean'] = np.mean(self.observations[observation]['residuals'])
                self.observations[observation]['res_std'] = np.std(self.observations[observation]['residuals'])
                self.observations[observation]['res_rms'] = np.sqrt(np.mean(self.observations[observation]['residuals'] ** 2))
                self.observations[observation]['res_chi_sqr'] = np.sum(
                    (self.observations[observation]['residuals']/self.observations[observation]['flux_unc']) ** 2)

                self.observations[observation]['detrended_res_mean'] = np.mean(self.observations[observation]['detrended_residuals'])
                self.observations[observation]['detrended_res_std'] = np.std(self.observations[observation]['detrended_residuals'])
                self.observations[observation]['detrended_res_rms'] = np.sqrt(np.mean(self.observations[observation]['detrended_residuals'] ** 2))
                self.observations[observation]['detrended_res_chi_sqr'] = np.sum(
                    (self.observations[observation]['detrended_residuals']/self.observations[observation]['detrended_flux_unc']) ** 2)

            fitting.results['observations'] = self.observations

            if output_folder:

                if not os.path.isdir(output_folder):
                    os.mkdir(output_folder)
                else:
                    shutil.rmtree(output_folder)
                    os.mkdir(output_folder)

                save_dict(results_copy, os.path.join(output_folder, 'results.pickle'))

                fitting.save_results(os.path.join(output_folder, 'results.txt'))

                fitting.plot_corner(os.path.join(output_folder, 'correlations.pdf'))

                fitting.plot_traces(os.path.join(output_folder, 'traces.pdf'))

                plot_transit_fitting_models(fitting.results, os.path.join(output_folder, 'lightcurves.pdf'))

                for observation in self.observations:

                    w = open(os.path.join(output_folder, 'diagnostics_dataset_{0}.txt'.format(observation + 1)), 'w')
                    w.write('\n#Residuals:\n')
                    w.write('#Mean: {0}\n'.format(self.observations[observation]['res_mean']))
                    w.write('#STD: {0}\n'.format(self.observations[observation]['res_std']))
                    w.write('#RMS: {0}\n'.format(self.observations[observation]['res_rms']))
                    w.write('#Chi squared: {0}\n'.format(self.observations[observation]['res_chi_sqr']))

                    w.write('\n\n#Detrended Residuals:\n')
                    w.write('#Mean: {0}\n'.format(self.observations[observation]['detrended_res_mean']))
                    w.write('#STD: {0}\n'.format(self.observations[observation]['detrended_res_std']))
                    w.write('#RMS: {0}\n'.format(self.observations[observation]['detrended_res_rms']))
                    w.write('#Chi squared: {0}\n'.format(self.observations[observation]['detrended_res_chi_sqr']))

                    w.close()

        fitting.results = results_copy

        return fitting

    def visibility(self, latitude, longitude, time_zone=0, horizon=0, start_time=now(), window=1, oot=2/24.0,
                   target_min_altitude=Degrees(20), sun_max_altitude=Degrees(-18),
                   max_moon_illumination=0.9, min_moon_distance=Degrees(30)):

        latitude = _reformat_or_request_angle(latitude)
        longitude = _reformat_or_request_angle(longitude)
        target_min_altitude = _reformat_or_request_angle(target_min_altitude)
        sun_max_altitude = _reformat_or_request_angle(sun_max_altitude)
        min_moon_distance = _reformat_or_request_angle(min_moon_distance)

        observatory = Observatory(latitude, longitude, time_zone, horizon)

        return observatory.periodic_events_visibility(
            self.target, start_time, window,
            self.mid_time, self.mid_time_format, self.period, self.transit_duration() + oot,
            target_min_altitude, sun_max_altitude, max_moon_illumination, min_moon_distance
        )
