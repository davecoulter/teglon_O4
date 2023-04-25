import matplotlib

matplotlib.use("Agg")

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.pyplot import cm
from matplotlib.patches import CirclePolygon
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as path_effects


import sys

sys.path.append('../')

import os
import optparse

from configparser import RawConfigParser
import multiprocessing as mp
import mysql.connector

import mysql.connector as test

print(test.__version__)

from mysql.connector.constants import ClientFlag
from mysql.connector import Error
import csv
import time
import pickle
from collections import OrderedDict

import numpy as np
from scipy.special import erf
from scipy.optimize import minimize, minimize_scalar
import scipy.stats as st
from scipy.integrate import simps
from scipy.interpolate import interp2d

import astropy as aa
from astropy import cosmology
from astropy.cosmology import LambdaCDM # WMAP5, WMAP7,
from astropy.coordinates import Distance
from astropy.coordinates.angles import Angle
from astropy.cosmology import z_at_value
from astropy import units as u
import astropy.coordinates as coord
# from dustmaps.config import config
from dustmaps.config import config as dustmaps_config
from dustmaps.sfd import SFDQuery

import shapely as ss
from shapely.ops import transform as shapely_transform
from shapely.geometry import Point
from shapely.ops import linemerge, unary_union, polygonize, split
from shapely import geometry
from shapely.geometry import JOIN_STYLE

import healpy as hp
from ligo.skymap import distance

# from Detector import *
# from HEALPix_Helpers import *
# from Tile import *
# from SQL_Polygon import *
# from Pixel_Element import *
# from Completeness_Objects import *

from web.src.objects.Detector import *
from web.src.objects.Tile import *
# from src.objects.SQL_Polygon import *
from web.src.objects.Pixel_Element import *
from web.src.utilities.Database_Helpers import *


import psutil
import shutil
import urllib.request
import requests
from bs4 import BeautifulSoup
from dateutil.parser import parse

import glob
import gc
import json

import MySQLdb as my

# from spherical_geometry.polygon import SphericalPolygon
# import ephem
from datetime import datetime, timezone, timedelta

from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS
# from astropy import units as unit
from scipy import interpolate
from astropy.time import Time
from astropy.coordinates import get_sun
import matplotlib.animation as animation
from matplotlib.pyplot import cm

import ligo.skymap.plot
from ligo.skymap.postprocess.util import find_greedy_credible_levels

isDEBUG = False


# Set up dustmaps config
utilities_base_dir = "./web/src/utilities"
dustmaps_config["data_dir"] = utilities_base_dir

# Generate all pixel indices
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)


def get_healpix_pixel_id(galaxy_info):
    phi, theta = np.radians(float(galaxy_info[8])), 0.5 * np.pi - np.radians(float(galaxy_info[9]))

    # map NSIDE is last argument of galaxy_info
    # return the galaxy_id with the pixel index in the NSIDE of the healpix map
    return (galaxy_info[0], hp.ang2pix(int(galaxy_info[-1]), theta, phi))


class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--gw_id', default="", type="str",
                          help='LIGO superevent name, e.g. `S190425z` ')

        parser.add_option('--healpix_dir', default='./web/events/{GWID}', type="str",
                          help='Directory for where to look for the healpix file.')

        parser.add_option('--healpix_file', default="", type="str", help='Healpix filename.')

        parser.add_option('--tele', default="s", type="str",
                          help='''Telescope abbreviation for the telescope to extract files for. Default: `s`. 
                          Available: (s:Swope, t:Thacher)
                          ''')

        parser.add_option('--band', default="r", type="str",
                          help='Telescope abbreviation for the telescope to extract files for. Default: `r`. Available: (g, r, i)')

        parser.add_option('--extinct', default="0.5", type="float",
                          help='Extinction in mags in the specified band to be less than. Default: 0.5 mag. Must be > 0.0')

        parser.add_option('--tile_file', default="{FILENAME}", type="str", help='Filename of tiles to plot.')

        parser.add_option('--galaxies_file', default="{FILENAME}", type="str", help='Optional: Filename of galaxies to plot.')

        parser.add_option('--get_tiles_from_db', action="store_true", default=False,
                          help='''Ignore tile file and get all observed tiles for map from database''')

        parser.add_option('--prob_type', default="2D", type="str",
                          help='''Probability type to consider. Default `2D`. Available: (4D, 2D)''')

        parser.add_option('--cum_prob_outer', default="0.9", type="float",
                          help='''Cumulative prob to cover in tiles. Default 0.9. Must be > 0.2 and < 0.95 
                          and > `cum_prob_inner`''')

        parser.add_option('--cum_prob_inner', default="0.5", type="float",
                          help='''Cumulative prob to cover in tiles. Default 0.5. Must be > 0.2 and < 0.95 
                          and < `cum_prob_inner`''')

        return (parser)

    def main(self):

        fig = plt.figure()
        ax_h = plt.axes(projection='astro hours mollweide')
        ax_h.grid()

        # (orig_prob, orig_distmu, orig_distsigma, orig_distnorm), header_gen = hp.read_map("./web/events/S190425z/GW190425_PublicationSamples_flattened.fits.gz", field=(0, 1, 2, 3),
        #                                                                                 h=True)
        # _90_50_levels = find_greedy_credible_levels(orig_prob)

        healpix_map_select = "SELECT id, RescaledNSIDE, t_0 FROM HealpixMap WHERE GWID = '%s' and Filename = '%s'"
        healpix_map_result = query_db([healpix_map_select % (self.options.gw_id, self.options.healpix_file)])[0][0]
        healpix_map_id = int(healpix_map_result[0])
        healpix_map_event_time = Time(healpix_map_result[2], format="gps")
        healpix_map_nside = int(healpix_map_result[1])

        select_2D_pix = '''
            SELECT
                hp.id,
                hp.HealpixMap_id,
                hp.Pixel_Index,
                hp.Prob,
                hp.Distmu,
                hp.Distsigma,
                hp.Mean,
                hp.Stddev,
                hp.Norm,
                hp.N128_SkyPixel_id
            FROM HealpixPixel hp
            WHERE hp.HealpixMap_id = %s
        '''

        pixels_to_select = select_2D_pix % healpix_map_id
        map_pix_result = query_db([pixels_to_select])[0]

        map_pix_dict = {int(mpr[2]): float(mpr[3]) for mpr in map_pix_result }

        num_pix = hp.nside2npix(healpix_map_nside)
        # map_pix = list(np.arange(num_pix))
        map_pix = np.zeros(num_pix)
        for i, mp in enumerate(map_pix):
            if i in map_pix_dict:
                map_pix[i] = map_pix_dict[i]
        # for mpr in map_pix_result:
        #     map_pix[int(mpr[2])] = float(mpr[3])

        _90_50_levels = find_greedy_credible_levels(np.asarray(map_pix))

        # Probability
        ax_h.contourf_hpx(_90_50_levels, cmap='OrRd', levels=[0.0, 0.5, 0.9], linewidths=0.2, alpha=0.3)
        ax_h.contour_hpx(_90_50_levels, cmap='OrRd', levels=[0.0, 0.5, 0.9], linewidths=0.2, alpha=0.3)
        # ax_h.contourf_hpx(_90_50_levels, cmap='OrRd', linewidths=0.2, alpha=0.3)
        # ax_h.contour_hpx(_90_50_levels, cmap='OrRd', linewidths=0.2, alpha=0.3)

        # Dust
        select_mwe_pix = '''
            SELECT
                sp.Pixel_Index, 
                sp_ebv.EBV*(SELECT b.F99_Coefficient FROM Band b WHERE `Name`='SDSS r')
            FROM
                SkyPixel sp
            JOIN SkyPixel_EBV sp_ebv on sp_ebv.N128_SkyPixel_id = sp.id
            # WHERE 
            #     sp_ebv.EBV*(SELECT b.F99_Coefficient FROM Band b WHERE `Name`='SDSS r') > 0.5
            ORDER BY 
                sp.Pixel_Index;
        '''

        mwe_pix_result = query_db([select_mwe_pix])[0]
        mwe_hpx_arr = []
        for mp in mwe_pix_result:
            mwe_hpx_arr.append(float(mp[1]))
        ax_h.contour_hpx(np.asarray(mwe_hpx_arr), colors='dodgerblue', levels=[0.5, np.max(mwe_hpx_arr)],
                         linewidths=0.2, alpha=0.5)
        ax_h.contourf_hpx(np.asarray(mwe_hpx_arr), cmap='Blues', levels=[0.5, np.max(mwe_hpx_arr)], linewidths=0.2,
                          alpha=0.3)


        # Pointings
        # telescope_name = "SWOPE"
        telescope_name = "T80S_T80S-Cam"
        detector_select_by_name = "SELECT id, Name, Deg_width, Deg_height, Deg_radius, Area, MinDec, MaxDec, ST_AsText(Poly) FROM Detector WHERE Name='%s'"
        detector_result = query_db([detector_select_by_name % telescope_name])[0][0]
        detector_id = int(detector_result[0])
        detector_name = detector_result[1]
        detector_poly = detector_result[8]
        detector_vertices = Detector.get_detector_vertices_from_teglon_db(detector_poly)
        detector = Detector(detector_name, detector_vertices, detector_id=detector_id)
        #
        formatted_healpix_dir = self.options.healpix_dir
        formatted_healpix_dir = formatted_healpix_dir.replace("{GWID}", self.options.gw_id)
        tile_file_path = "%s/%s" % (formatted_healpix_dir, self.options.tile_file)
        tiles_to_plot = []



        with open('%s' % tile_file_path, 'r') as csvfile:
            csvreader = csv.reader(csvfile, delimiter=',', skipinitialspace=True)
            # Skip Header
            next(csvreader)

            for row in csvreader:
                name = row[0]
                c = coord.SkyCoord(row[1], row[2], unit=(u.hour, u.deg))
                t = Tile(c.ra.degree, c.dec.degree, detector=detector, nside=healpix_map_nside)

                t.field_name = name
                t.net_prob = float(row[6])
                tiles_to_plot.append(t)

        from matplotlib import colors
        tile_probs = [t.net_prob for t in tiles_to_plot[0:200]]
        min_prob = np.min(tile_probs)
        max_prob = np.max(tile_probs)
        print("min prob: %s" % min_prob)
        print("max prob: %s" % max_prob)

        clr_norm = colors.LogNorm(min_prob, max_prob)

        # This is an arbitrary # to plot... but they're sorted by highest prob first, and we'll never get 1000 pointings
        for i, t in enumerate(tiles_to_plot[0:200]):
            # HACK - protecting against plotting errs
            if t.ra_deg < 358 and t.ra_deg > 2:
                t.plot2(ax_h, edgecolor='k',
                        facecolor=plt.cm.Greens(clr_norm(t.net_prob)),
                        linewidth=0.1, alpha=1.0, zorder=9999) #


        # Sun
        from astroplan import Observer, AirmassConstraint
        from astropy.coordinates import AltAz, EarthLocation

        time_of_trigger = Time(healpix_map_event_time.to_datetime(), scale="utc")
        observatory = Observer.at_site("LCO")
        sun_set = observatory.sun_set_time(time_of_trigger, which='nearest', horizon=-18 * u.deg)
        sun_rise = observatory.sun_rise_time(sun_set, which='next', horizon=-18 * u.deg)
        sun_coord = get_sun(sun_set)

        hours_of_the_night = int((sun_rise - sun_set).to_value('hr'))
        observatory_location = EarthLocation(lat=-29.0146 * u.deg, lon=-70.6926 * u.deg, height=2380 * u.m)
        time_next = sun_set

        theta, phi = hp.pix2ang(nside=64, ipix=np.arange(hp.nside2npix(64)))
        ra = np.rad2deg(phi)
        dec = np.rad2deg(0.5 * np.pi - theta)
        all_sky_coords = coord.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))

        sun_separations = sun_coord.separation(all_sky_coords)
        deg_sep = np.asarray([sp.deg for sp in sun_separations])

        ax_h.contour_hpx(deg_sep, colors='darkorange', levels=[0, 60.0], alpha=0.5, linewidths=0.2)
        ax_h.contourf_hpx(deg_sep, colors='gold', levels=[0, 60.0], linewidths=0.2, alpha=0.3)
        # fmt = {}
        # strs = [r'$60\degree$']
        # for l, s in zip(sun_avoidance_contour.levels, strs):
        #     fmt[l] = s
        # c_labels = ax_h.clabel(sun_avoidance_contour, inline=True, fontsize=10, fmt=fmt,
        #                      levels=sun_avoidance_contour.levels, colors='k', use_clabeltext=True,
        #                      inline_spacing=60, zorder=9999)

        ax_h.plot(sun_coord.ra.degree, sun_coord.dec.degree, transform=ax_h.get_transform('world'), marker="o",
                  markersize=10, markeredgecolor="darkorange", markerfacecolor="gold", linewidth=0.05)


        # Swope Airmass constraints
        net_airmass = all_sky_coords.transform_to(AltAz(obstime=sun_set, location=observatory_location)).secz
        for i in np.arange(hours_of_the_night):
            time_next += timedelta(1/24.)
            temp = all_sky_coords.transform_to(AltAz(obstime=time_next, location=observatory_location)).secz
            temp_index = np.where((temp >= 1.0) & (temp <= 2.0))
            net_airmass[temp_index] = 1.5

        ax_h.contourf_hpx(net_airmass, colors='gray', levels=[np.min(net_airmass), 1.0], linewidths=0.2, alpha=0.3)
        ax_h.contourf_hpx(net_airmass, colors='gray', levels=[2.0, np.max(net_airmass)], linewidths=0.2, alpha=0.3)



        fig.savefig("./web/events/S190425z/test_fig.svg", bbox_inches='tight', format="svg")  # ,dpi=840
        plt.close('all')
        print("... Done.")




if __name__ == "__main__":
    useagestring = """python Plot_Teglon.py [options]

python Plot_Teglon.py --gw_id <gwid> --healpix_file <filename> --extinct 0.50 --prob_type 2D --cum_prob_outer 0.90 
--cum_prob_inner 0.50 --band r --tile_file <FILENAME> --tele s --galaxies_file <FILENAME> 
"""

    start = time.time()

    teglon = Teglon()
    parser = teglon.add_options(usage=useagestring)
    options, args = parser.parse_args()
    teglon.options = options

    teglon.main()

    end = time.time()
    duration = (end - start)
    print("\n********* start DEBUG ***********")
    print("Teglon `Plot_Teglon` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")


