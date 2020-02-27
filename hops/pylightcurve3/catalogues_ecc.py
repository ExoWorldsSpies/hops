from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
import time
import os
from .tools_files import open_dict
from ._1databases import databases
from .astro_angles_and_time import Hours, Degrees, Target


vizier = Vizier(columns=["*", "+_r"], catalog="II/246")
Simbad.add_votable_fields('ids', 'flux(V)', 'flux(R)', 'flux(I)', 'flux(J)', 'flux(H)', 'flux(K)', 'parallax', 'pmdec', 'pmra')


ecc_planets = open_dict(os.path.join(databases.oec, 'planets_catalogue.pickle'))
ecc_stars = open_dict(os.path.join(databases.oec, 'stars_catalogue.pickle'))

def flat_name(name):

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


def search_by_planet(name):

    if name in ecc_planets:
        return name

    else:
        for i in ecc_planets:
            if flat_name(i) == flat_name(name):
                return str(i)
            elif flat_name(i)[-1] == flat_name(name)[-1] and flat_name(i)[:-2] == flat_name(name)[:-1] and flat_name(i)[-2] in ['a', 'b', 'n']:
                return str(i)

    raise ValueError('No planet found in the catalogue')


def search_by_star(name):

    name_init = name
    found = False

    if name in ecc_stars:
        found = True

    if not found:

        for i in ecc_stars:
            if flat_name(i) == flat_name(name):
                name = i
                found = True
                break

    if not found:

        name = name.replace('HD-', 'HD ')
        name = name.replace('2MASS-', '2MASS ')
        name = name.replace('GJ-', 'GJ ')
        name = name.replace('HIP-', 'HIP ')
        name = name.replace('LHS-', 'LHS ')
        name = name.replace('TYC-', 'TYC ')
        name = name.replace('KIC-', 'KIC ')

        if 'BD' in name:
            name = name[:3] + name[3:].replace('-', ' ')

        result_table = None
        connected = False
        time.sleep(0.01)
        while not connected:
            try:
                result_table = Simbad.query_object(name)
                connected = True
            except:
                print('Retry to connect.')
                time.sleep(0.01)
                pass

        if result_table is not None:
            for star in ecc_stars:
                if ecc_stars[star]['simbad_id'] == result_table[0]['MAIN_ID'].decode("utf-8"):
                    name = star
                    found = True
                    break

    if not found:
        raise ValueError('No star found in the catalogue (search: {0})'.format(name_init))

    if len(ecc_stars[name]['planets']) == 1:
        return ecc_stars[name]['planets'][0]
    else:
        raise ValueError('More than one planets found in the catalogue for star {0} (Simbad name: {1}, '
                         'search: {2}): {3}'.format(name, ecc_stars[name]['simbad_id'], name_init,
                                                    ', '.join(ecc_stars[name]['planets'])))


class _Star:

    def __init__(self, planet):
        self.planet_name = planet
        star = planet[:-1]
        self.name = star
        self.planets = ecc_stars[star]['planets']
        self.simbad_id = ecc_stars[star]['simbad_id']
        self.simbad_link = ecc_stars[star]['simbad_link']
        self.gaia_id = ecc_stars[star]['gaia_id']
        self.twomass_id = ecc_stars[star]['twomass_id']
        self.ra = ecc_stars[star]['ra']
        self.dec = ecc_stars[star]['dec']
        self.v_mag = ecc_stars[star]['v_mag']
        self.r_mag = ecc_stars[star]['r_mag']
        self.i_mag = ecc_stars[star]['i_mag']
        self.v_mag_calc = ecc_stars[star]['v_mag_calc']
        self.r_mag_calc = ecc_stars[star]['r_mag_calc']
        self.i_mag_calc = ecc_stars[star]['i_mag_calc']
        self.gaia_ra = ecc_stars[star]['gaia_ra']
        self.gaia_dec = ecc_stars[star]['gaia_dec']
        self.gaia_g_mag = ecc_stars[star]['gaia_g_mag']
        self.gaia_bp_mag = ecc_stars[star]['gaia_bp_mag']
        self.gaia_rp_mag = ecc_stars[star]['gaia_rp_mag']
        self.gaia_parallax = ecc_stars[star]['gaia_parallax']
        self.gaia_pm_ra = ecc_stars[star]['gaia_pm_ra']
        self.gaia_pm_dec = ecc_stars[star]['gaia_pm_dec']
        self.j_mag = ecc_stars[star]['j_mag']
        self.h_mag = ecc_stars[star]['h_mag']
        self.k_mag = ecc_stars[star]['k_mag']
        self.teff = ecc_stars[star]['teff']
        self.logg = ecc_stars[star]['logg']
        self.met = ecc_stars[star]['met']
        self.mass = ecc_stars[star]['mass']
        self.rad = ecc_stars[star]['rad']

    def __str__(self):
        return 'plc.ExoPlanet.star(Simbad ID: {0}, host of {1})'.format(self.simbad_id, self.planet_name)

    def __repr__(self):
        return self.__str__()


class _Planet:
    def __init__(self, planet):
        self.name = planet
        self.rp_over_rs = ecc_planets[planet]['rp_over_rs']
        self.sma_over_rs = ecc_planets[planet]['sma_over_rs']
        self.inclination = ecc_planets[planet]['inclination']
        self.eccentricity = ecc_planets[planet]['eccentricity']
        self.periastron = ecc_planets[planet]['periastron']
        self.reference = None
        # self.mass =
        # self.radius =
        self.mid_time = ecc_planets[planet]['ephemeris']['mid_time']
        self.mid_time_error = ecc_planets[planet]['ephemeris']['mid_time_error']
        self.period = ecc_planets[planet]['ephemeris']['period']
        self.period_error = ecc_planets[planet]['ephemeris']['period_error']
        self.time_format = ecc_planets[planet]['ephemeris']['time_format']
        self.ephemeris_reference = ecc_planets[planet]['ephemeris']['reference']

    def __str__(self):
        return 'plc.ExoPlanet.planet(NEC ID: {0})'.format(self.name)

    def __repr__(self):
        return self.__str__()


class ExoPlanet:

    def __init__(self, planet=None, star=None):

        if planet:

            planet = search_by_planet(planet)

        elif star:

            planet = search_by_star(star)

        self.planet = _Planet(planet)
        self.star = _Star(planet)

    def __str__(self):
        return 'plc.ExoPlanet(planet NEC ID: {0}, star Simbad ID: {1})'.format(self.planet.name, self.star.simbad_id)

    def __repr__(self):
        return self.__str__()


def find_nearest(fov_coord):

    stars = [[ff, Target(Hours(ecc_stars[ff]['ra']), Degrees(ecc_stars[ff]['dec']))] for ff in ecc_stars]

    stars = sorted(stars, key=lambda x: fov_coord.distance_from_target(x[1]).deg)

    star = stars[0][0]

    return ExoPlanet(ecc_stars[star]['planets'][0])
