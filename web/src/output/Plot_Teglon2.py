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

from src.objects.Detector import *
from src.objects.Tile import *
# from src.objects.SQL_Polygon import *
from src.objects.Pixel_Element import *
from src.utilities.Database_Helpers import *


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

from spherical_geometry.polygon import SphericalPolygon
import ephem
from datetime import datetime, timezone, timedelta

from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS
# from astropy import units as unit
from scipy import interpolate
from astropy.time import Time
from astropy.coordinates import get_sun
import matplotlib.animation as animation
from matplotlib.pyplot import cm

isDEBUG = False


# Set up dustmaps config
dustmaps_config["data_dir"] = "./"

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

        parser.add_option('--healpix_dir', default='../../Events/{GWID}', type="str",
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

        map_pix = query_db(["SELECT Pixel_Index, Prob FROM HealpixPixel Order By Pixel_Index"])[0]
        allsky = np.arange(0, 12*256**2)
        map_pix_arr = np.asarray([mp[0] for mp in map_pix])

        # Plot!!
        fig = plt.figure(figsize=(10, 10), dpi=800)
        ax = fig.add_subplot(111)
        m = Basemap(projection='moll', lon_0=180.0)

        hp.mollview(map_pix_arr, fig=plt.gcf().number, cmap=plt.cm.Greens)

        fig.savefig("test_fig.png", bbox_inches='tight')  # ,dpi=840
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


