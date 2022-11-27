
__all__ = ['get_planet', 'get_all_planets', 'locate_planet']


from ..errors import *
from ..__databases__ import plc_data
from ..models.exoplanet import Planet
from ..spacetime.angles import Degrees, Hours, _request_angle
from ..spacetime.targets import FixedTarget


def _flat_name(name):

    flat_name_list = [
        [' ', ''],
        ['-', ''],
        ['cancri', 'cnc'],
        ['hatp10', 'wasp11'],
        ['wasp40', 'hatp27'],
        ['wasp51', 'hatp30'],
        ['wasp86', 'kelt12'],
        ['kelt22', 'wasp173'],
    ]

    name = name.lower()

    for char in flat_name_list:
        name = name.replace(char[0], char[1])

    return name


def _search_by_planet(name):

    planets = plc_data.ecc()['planets']

    name_or = name

    if name in planets:
        return name

    else:
        for i in planets:
            if _flat_name(i) == _flat_name(name):
                return str(i)
            elif (_flat_name(i)[-1] == _flat_name(name)[-1] and _flat_name(i)[:-2] == _flat_name(name)[:-1] and
                  _flat_name(i)[-2] in ['a', 'b', 'n']):
                return str(i)

    raise PyLCInputError('No planet {0} found in the catalogue.'.format(name_or))


def get_all_planets():

    return list(plc_data.ecc()['planets'].keys())


def get_planet(name):

    name = _search_by_planet(name)

    planet_data = plc_data.ecc()['planets'][name]
    star_data = plc_data.ecc()['stars'][planet_data['star']]

    planet = Planet(
        name,
        Hours(star_data['ra']),
        Degrees(star_data['dec']),
        planet_data['logg'],
        planet_data['teff'],
        planet_data['meta'],
        planet_data['rp_over_rs'],
        planet_data['ephem_period'],
        planet_data['sma_over_rs'],
        planet_data['eccentricity'],
        planet_data['inclination'],
        planet_data['periastron'],
        planet_data['ephem_mid_time'],
        'BJD_TDB',
    )

    planet.all_data = {'planet': planet_data, 'star': star_data}

    return planet


def locate_planet(ra, dec, radius=Degrees(0, 1, 0)):

    # for compatibility reasons will be deprecated...
    try:
        ra = float(ra)
        dec = float(dec)
        radius = float(radius)
        ra = Degrees(ra)
        dec = Degrees(dec)
        radius = Degrees(radius)
    except:
        pass

    pointing = FixedTarget(ra, dec)
    _request_angle(radius)

    test_planets = []

    for test_planet_name in get_all_planets():
        test_planet = get_planet(test_planet_name).target
        test_planets.append([pointing.distance_on_sphere(test_planet).deg(), test_planet_name])

    test_planets.sort()

    if test_planets[0][0] < radius.deg():
        return get_planet(test_planets[0][1])
    else:
        raise PyLCLibraryError('Planet could not be located')
