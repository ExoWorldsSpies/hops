
__all__ = ['find_fits_files', 'test_fits_keyword', 'test_coordinates', 'test_file_number', 'test_float_input',
           'test_float_positive_input', 'test_int_input', 'test_int_positive_input', 'test_int_positive_non_zero_input',
           'test_date', 'filter_map']

import numpy as np
import glob
import hops.pylightcurve3 as plc
from astropy.io import fits as pf


def find_fits_files(fits_file):

    fits_list = glob.glob('*{0}*.f*t*'.format(fits_file)) + glob.glob('*{0}*.F*T*'.format(fits_file))
    fits_list = list(np.unique(fits_list))
    fits_list.sort()
    return fits_list


def test_fits_keyword(fits_file, keyword):

    if len(fits_file) == 0:
        return [False, 'No keyword found']

    else:
        try:
            fits_file = find_fits_files(fits_file)[0]

            fits = pf.open(fits_file)

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

            if fits[0].header[str(keyword)]:
                return [True, 'Keyword found', fits[0].header[str(keyword)]]

            else:
                return [False, 'No keyword found']

        except (KeyError, IndexError):
            return [False, 'No keyword found']


def test_file_number(fits_file):

    if len(fits_file) == 0:
        test = 0
    else:
        test = len(find_fits_files(fits_file))

    if test > 0:
        return [True, '{0} files found'.format(test)]
    else:
        return [False, 'No files found']


def test_coordinates(ra_dec_string, single_line=False):

    try:
        ra_dec_string = ra_dec_string.split(' ')[0].split(':') + ra_dec_string.split(' ')[1].split(':')
        target = plc.Target(plc.Hours(ra_dec_string[0], ra_dec_string[1], ra_dec_string[2]),
                            plc.Degrees(ra_dec_string[3], ra_dec_string[4], ra_dec_string[5]))
        if single_line:
            return[True, 'Coordinates accepted']
        else:
            return[True, 'Coordinates\naccepted']
    except:
        if single_line:
            return [False, 'Wrong coordinates']
        else:
            return [False, 'Wrong\ncoordinates']


def test_float_input(input_str, typing):

    if typing == '1':
        try:
            if float(input_str) >= 0:
                return True
            elif float(input_str) < 0:
                return True
            else:
                return False
        except ValueError:
            return False

    else:
        return True


def test_float_positive_input(input_str, typing):

    if typing == '1':
        try:
            if float(input_str) >= 0:
                return True
            else:
                return False
        except ValueError:
            return False

    else:
        return True


def test_int_input(input_str, typing):

    if typing == '1':
        try:
            if int(input_str):
                return True
            else:
                return False
        except ValueError:
            return False

    else:
        return True


def test_int_positive_input(input_str, typing):

    if typing == '1':
        try:
            if int(input_str) >= 0:
                return True
            else:
                return False
        except ValueError:
            return False

    else:
        return True


def test_int_positive_non_zero_input(input_str, typing):

    if typing == '1':
        try:
            if int(input_str) > 0:
                return True
            else:
                return False
        except ValueError:
            return False

    else:
        return True

def test_date(year_month_string):

    pass_test = True

    if len(year_month_string.split()) != 2:
        pass_test = False

    elif len(year_month_string.split()[0]) != 4:
        pass_test = False

    elif len(year_month_string.split()[1]) != 2:
        pass_test = False
    else:
        try:
            year = int(year_month_string.split()[0])
            month = int(year_month_string.split()[1])
            if int(month) < 1 or int(month) > 12:
                pass_test = False
            if int(year) < 0:
                pass_test = False
        except:
            pass_test = False

    if pass_test:
        return[True, 'Date\naccepted']
    else:
        return [False, 'Wrong\ndate']

filter_map = {'Clear': 'V', 'Luminance': 'V',
              'U': 'U', 'B': 'B', 'V': 'V', 'R': 'R', 'I': 'I', 'H': 'H', 'J': 'J', 'K': 'K',
              'u': 'u', 'b': 'b', 'v': 'v', 'y': 'y',
              'u\'': 'u,', 'g\'': 'g,', 'r\'': 'r,', 'i\'': 'i,', 'z\'': 'z,',
              'Astrodon ExoPlanet-BB': 'R',
              'UV': 'U', 'Rc': 'R', 'Ic': 'I', 'Re': 'R', 'Ie': 'I', 'Y': 'y,', 'r': 'r,', 'z': 'z,', 'i': 'i,',
              }

