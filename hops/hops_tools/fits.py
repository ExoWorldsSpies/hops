
__all__ = ['find_fits_files', 'get_fits_data']

import glob
import numpy as np

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



