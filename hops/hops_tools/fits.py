
__all__ = ['find_fits_files', 'get_fits_data', 'get_fits_data_and_header']

import glob
import numpy as np
from astropy.io import fits as pf

import hops.pylightcurve41 as plc


def find_fits_files(name_identifier):

    if len(name_identifier) > 0:

        fits_files_names = glob.glob('*{0}*.f*t*'.format(name_identifier)) + glob.glob('*{0}*.F*T*'.format(name_identifier))
        fits_files_names = list(np.unique(fits_files_names))
        fits_files_names.sort()
        return fits_files_names

    else:
        return []


def get_fits_data(fits_file_name):

    with plc.open_fits(fits_file_name) as fits_file:

        try:
            fits_data = [fits_file['SCI']]
        except KeyError:
            sci_id = 0
            for sci_id in range(len(fits_file)):
                try:
                    if fits_file[sci_id].data.all():
                        break
                except:
                    pass
            fits_data = [fits_file[sci_id]]

    return fits_data


def get_fits_data_and_header(path):

    with pf.open(path, memmap=False) as fits_file:

        fits_file.verify('fix')

        try:
            fits_data = np.array(fits_file['SCI'].data)
            fits_header = {ff: fits_file['SCI'].header[ff] for ff in fits_file['SCI'].header}
        except KeyError:
            for idx in range(len(fits_file)):
                try:
                    fits_data = np.array(fits_file[idx].data)
                    fits_header = {ff: fits_file[idx].header[ff] for ff in fits_file[idx].header if ff not in ['HISTORY', 'COMMENT', '']}
                    if isinstance(fits_file[idx].data, type(None)):
                        pass
                    else:
                        break
                except:
                    pass

        # bit10_test = np.sum((fits_data/64.0-np.int_(fits_data/64.0))**2) == 0
        # bit12_test = np.sum((fits_data/16.0-np.int_(fits_data/16.0))**2) == 0
        # bit14_test = np.sum((fits_data/4.0-np.int_(fits_data/4.0))**2) == 0
        #
        # if bit10_test:
        #     fits_data = fits_data/64.0
        #     fits_header['BITPIX'] = 10
        # elif bit12_test:
        #     fits_data = fits_data/16.0
        #     fits_header['BITPIX'] = 12
        # elif bit14_test:
        #     fits_data = fits_data/4.0
        #     fits_header['BITPIX'] = 14

        fits_data[np.where(np.isnan(fits_data))] = 1
        fits_data[np.where(fits_data == 0)] = 1

    return np.array(fits_data, dtype=float), fits_header


