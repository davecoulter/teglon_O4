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
from astropy.cosmology import WMAP5, WMAP7, LambdaCDM
from astropy.coordinates import Distance
from astropy.coordinates.angles import Angle
from astropy.cosmology import z_at_value
from astropy import units as u
import astropy.coordinates as coord
from dustmaps.config import config
from dustmaps.sfd import SFDQuery

import shapely as ss
from shapely.ops import transform as shapely_transform
from shapely.geometry import Point
from shapely.ops import linemerge, unary_union, polygonize, split
from shapely import geometry

import healpy as hp
from ligo.skymap import distance

from HEALPix_Helpers import *
from Tile import *
from SQL_Polygon import *
from Pixel_Element import *
from Completeness_Objects import *

import psutil
import shutil
import urllib.request
import requests
from bs4 import BeautifulSoup
from dateutil.parser import parse

import glob
import gc
import json

from astropy.table import Table



isDEBUG = False



# Generate all pixel indices
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)



class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--gw_id', default="", type="str", help='LIGO superevent name, e.g. `S190425z`.')

        parser.add_option('--healpix_file', default="", type="str",
                          help='Healpix filename. Used with `gw_id` to identify unique map.')

        parser.add_option('--healpix_dir', default='../Events/{GWID}', type="str",
                          help='Directory for where to look for the healpix file.')

        parser.add_option('--output_dir', default="../Events/{GWID}/ObservedTiles/TileStats", type="str",
                          help='Directory for where to look for observed tiles to import.')

        return (parser)

    def main(self):

        is_error = False

        # Parameter checks
        if self.options.gw_id == "":
            is_error = True
            print("GWID is required.")

        formatted_healpix_dir = self.options.healpix_dir
        if "{GWID}" in formatted_healpix_dir:
            formatted_healpix_dir = formatted_healpix_dir.replace("{GWID}", self.options.gw_id)

        formatted_output_dir = self.options.output_dir
        if "{GWID}" in formatted_output_dir:
            formatted_output_dir = formatted_output_dir.replace("{GWID}", self.options.gw_id)

        hpx_path = "%s/%s" % (formatted_healpix_dir, self.options.healpix_file)

        if self.options.healpix_file == "":
            is_error = True
            print("You must specify which healpix file to process.")

        if is_error:
            print("Exiting...")
            return 1

        # Get Map ID
        healpix_map_select = "SELECT id FROM HealpixMap WHERE GWID = '%s' and Filename = '%s'"
        healpix_map_id = int(query_db([healpix_map_select % (self.options.gw_id, self.options.healpix_file)])[0][0][0])

        # Get all conposed observed tiles for this map
        observed_tile_select = '''
            SELECT 
                ot.FileName,
                ot.id as ObservedTile_id, 
                ot.FieldName, 
                ot.RA, 
                ot._Dec, 
                ot.MJD, 
                ot.Exp_Time, 
                ot.Mag_Lim, 
                d.Name as Detector_Name, 
                ot.EBV, 
                b.Name as Band_Name, 
                ot.EBV*b.F99_Coefficient as A_lambda, 
                SUM(hp.Prob) as LIGO_2D_Tile_Prob, 
                SUM(hpc.Renorm2DProb) as Renorm_2D_Tile_Prob, 
                SUM(hpc.NetPixelProb) as Teglon_4D_Tile_Prob, # This should equal the contribution from the galaxies + redistributed prob
                AVG(hpc.PixelCompleteness) as Avg_Tile_Completeness, 
                AVG(hp.Mean) as Avg_LIGO_Distance 
            FROM ObservedTile ot 
            JOIN Detector d on d.id = ot.Detector_id 
            JOIN Band b on b.id = ot.Band_id 
            JOIN ObservedTile_HealpixPixel ot_hp on ot_hp.ObservedTile_id = ot.id 
            JOIN HealpixPixel hp on hp.id = ot_hp.HealpixPixel_id 
            JOIN HealpixPixel_Completeness hpc on hpc.HealpixPixel_id = hp.id 
            WHERE ot.HealpixMap_id = %s 
            GROUP BY 
                ot.FileName,
                ot.id, 
                ot.FieldName, 
                ot.RA, 
                ot._Dec, 
                ot.MJD, 
                ot.Exp_Time, 
                ot.Mag_Lim, 
                d.Name, 
                ot.EBV, 
                b.Name, 
                ot.EBV*b.F99_Coefficient 
            ORDER BY ot.id 
        '''
        observed_tile_result = query_db([observed_tile_select % healpix_map_id])[0]

        observed_tile_galaxy_select = '''
            SELECT
                ot.id as ObservedTile_id, 
                ot.FieldName, 
                ot.RA as ObservedTile_RA, 
                ot._Dec as ObservedTile_Dec, 
                hp.id as HealpixPixel_id, 
                gd2.id as Galaxy_id, 
                gd2.RA as Galaxy_RA, 
                gd2._Dec as Galaxy_Dec, 
                gd2.PGC, 
                gd2.Name_GWGC, 
                gd2.Name_HyperLEDA, 
                gd2.Name_2MASS, 
                gd2.Name_SDSS_DR12, 
                hp_gd2_w.GalaxyProb as Galaxy_4D_Prob, 
                gd2.z, 
                gd2.z_dist, 
                gd2.z_dist_err, 
                gd2.B, 
                gd2.K 
            FROM 
                GalaxyDistance2 gd2 
            JOIN 
                HealpixPixel_GalaxyDistance2 hp_gd2 on  hp_gd2.GalaxyDistance2_id = gd2.id 
            JOIN 
                HealpixPixel_GalaxyDistance2_Weight hp_gd2_w on hp_gd2_w.HealpixPixel_GalaxyDistance2_id = hp_gd2.id 
            JOIN 
                HealpixPixel hp on hp.id = hp_gd2.HealpixPixel_id 
            JOIN 
                HealpixMap hm on hm.id = hp.HealpixMap_id 
            JOIN 
                ObservedTile_HealpixPixel ot_hp on ot_hp.HealpixPixel_id = hp.id 
            JOIN 
                ObservedTile ot on ot.id = ot_hp.ObservedTile_id 
            WHERE 
                ot.id = %s AND hm.id = %s 
        '''

        observed_tile_pixel_select = '''
            SELECT 
                hp.Prob, 
                hpc.Renorm2DProb, 
                hpc.NetPixelProb, 
                hpc.PixelCompleteness, 
                hp.Mean 
            FROM ObservedTile ot 
            JOIN Detector d on d.id = ot.Detector_id 
            JOIN Band b on b.id = ot.Band_id 
            JOIN ObservedTile_HealpixPixel ot_hp on ot_hp.ObservedTile_id = ot.id 
            JOIN HealpixPixel hp on hp.id = ot_hp.HealpixPixel_id 
            JOIN HealpixPixel_Completeness hpc on hpc.HealpixPixel_id = hp.id 
            WHERE ot.id = %s 
        '''

        for ot in observed_tile_result:
            ot_id = ot[1]
            observed_tile_galaxy_result = query_db([observed_tile_galaxy_select % (ot_id, healpix_map_id)])[0]
            observed_tile_pixel_result = query_db([observed_tile_pixel_select % ot_id])[0]

            # Get weighted mean for distance and completeness based on the 4D prob
            pix_weights = np.asarray([float(otpr[2]) for otpr in observed_tile_pixel_result])
            pix_comp = np.asarray([float(otpr[3]) for otpr in observed_tile_pixel_result])
            pix_dist = np.asarray([float(otpr[4]) for otpr in observed_tile_pixel_result])

            avg_mean_dist = np.average(pix_dist)
            avg_pix_comp = np.average(pix_comp)

            # If we can weight the average, do so...
            if np.sum(pix_weights) > 0.0:
                avg_mean_dist = np.average(pix_dist, weights=pix_weights)
                avg_pix_comp = np.average(pix_comp, weights=pix_weights)

            field_name = ot[2]
            ascii_ecsv_fname = "%s_%s.txt" % (ot_id, field_name)
            ascii_ecsv_fpath = "%s/%s" % (formatted_output_dir, ascii_ecsv_fname)
            print("Creating `%s`" % ascii_ecsv_fpath)

            # Build ascii.ecsv formatted output
            cols = ['Galaxy_ID',
                    'Galaxy_RA',
                    'Galaxy_Dec',
                    'PGC',
                    'Name_GWGC',
                    'Name_HyperLEDA',
                    'Name_2MASS',
                    'Name_SDSS_DR12',
                    'Galaxy_4D_Prob',
                    'z',
                    'z_Dist',
                    'z_Dist_Err',
                    'B',
                    'K']
            dtype = ['i4', 'f8', 'f8', 'U64', 'U64', 'U64', 'U64', 'U64', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8']
            result_table = Table(dtype=dtype, names=cols)
            meta = ["{key}={value}".format(key="DCMP File", value=ot[0]),
                    "{key}={value}".format(key="DB ID", value=ot[1]),
                    "{key}={value}".format(key="Field Name", value=ot[2]),
                    "{key}={value}".format(key="RA (deg)", value=ot[3]),
                    "{key}={value}".format(key="Dec (deg)", value=ot[4]),
                    "{key}={value}".format(key="MJD", value=ot[5]),
                    "{key}={value}".format(key="Exp Time (sec)", value=ot[6]),
                    "{key}={value}".format(key="Lim Mag", value=ot[7]),
                    "{key}={value}".format(key="Instrument", value=ot[8]),
                    "{key}={value}".format(key="E(B-V)", value=ot[9]),
                    "{key}={value}".format(key="Filter", value=ot[10]),
                    "{key}={value}".format(key="A_lambda", value=ot[11]),
                    "{key}={value}".format(key="Net 2D Prob", value=ot[12]),
                    "{key}={value}".format(key="Redistributed 2D Prob", value=ot[13]),
                    "{key}={value}".format(key="Net 4D Prob", value=ot[14]),
                    "{key}={value}".format(key="Weighted Mean Pixel Comp", value=avg_pix_comp),
                    "{key}={value}".format(key="Weighted Mean Pixel Dist", value=avg_mean_dist)
                    ]
            result_table.meta['comment'] = meta
            for otgr in observed_tile_galaxy_result:
                result_table.add_row([otgr[5],
                                        otgr[6],
                                        otgr[7],
                                        otgr[8],
                                        otgr[9],
                                        otgr[10],
                                        otgr[11],
                                        otgr[12],
                                        otgr[13],
                                        otgr[14],
                                        otgr[15],
                                        otgr[16],
                                        otgr[17],
                                        otgr[18]])
            result_table.write(ascii_ecsv_fpath, overwrite=True, format='ascii.ecsv')


if __name__ == "__main__":
    useagestring = """python Generate_Tiles.py [options]

Example with healpix_dir defaulted to 'Events/<gwid>':
python LoadMap.py --gw_id <>gwid --healpix_file <filename>

Example with healpix_dir specified:
python LoadMap.py --gw_id <gwid> --healpix_dir Events/<directory name> --healpix_file <filename>
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
    print("Teglon `LoadObservedTiles` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")
