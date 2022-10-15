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

# from Detector import *
# from HEALPix_Helpers import *
# from Tile import *
# from SQL_Polygon import *
# from Pixel_Element import *
# from Completeness_Objects import *
from src.objects.Detector import *
from src.objects.Tile import *
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

isDEBUG = False


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

        parser.add_option('--tile_dir', default="../Events/{GWID}/ObservedTiles", type="str",
                          help='Directory for where to look for observed tiles to import.')

        parser.add_option('--tile_file', default="", type="str", help='File that contains the tile observations.')

        parser.add_option('--tele', default="", type="str",
                          help='Telescope abbreviation that `tile_file` corresponds to.')

        parser.add_option('--coord_format', default="hour", type="str",
                          help='What format the sky coordinates are in. Available: "hour" or "deg". Default: "hour".')

        return (parser)

    def main(self):

        HOUR_FORMAT = "hour"
        DEG_FORMAT = "deg"

        # Band abbreviation, band_id mapping
        band_mapping = {
            "g": "SDSS g",
            "r": "SDSS r",
            "i": "SDSS i",
            "Clear": "Clear",
            "open": "Clear",
            "J": "UKIRT J",
            "B":"Landolt B",
            "V": "Landolt V",
            "R": "Landolt R",
            "I": "Landolt I",
            "u": "SDSS u"
        }

        detector_mapping = {
            "s": "SWOPE",
            "t": "THACHER",
            "a": "ANDICAM",
            "n": "NICKEL",
            "m": "MOSFIRE",
            "k": "KAIT",
            "si": "SINISTRO",
            "w": "WISE",
            "gt4": "GOTO-4",
            "css": "MLS10KCCD-CSS",
            "sw": "Swift/UVOT",
            "mmt": "MMTCam",
            "z": "ZTF",
        }

        is_error = False

        # Parameter checks
        if self.options.gw_id == "":
            is_error = True
            print("GWID is required.")

        if self.options.healpix_file == "":
            is_error = True
            print("Healpix file is required.")

        if self.options.tile_file == "":
            is_error = True
            print("Tile file is required.")

        if self.options.tele == "":
            is_error = True
            print("Telescope abbreviation is required.")

        if self.options.coord_format != HOUR_FORMAT and self.options.coord_format != DEG_FORMAT:
            is_error = True
            print("Incorrect `coord_format`!")

        if is_error:
            print("Exiting...")
            return 1

        formatted_healpix_dir = self.options.healpix_dir
        if "{GWID}" in formatted_healpix_dir:
            formatted_healpix_dir = formatted_healpix_dir.replace("{GWID}", self.options.gw_id)

        formatted_tile_dir = self.options.tile_dir
        if "{GWID}" in formatted_tile_dir:
            formatted_tile_dir = formatted_tile_dir.replace("{GWID}", self.options.gw_id)

        hpx_path = "%s/%s" % (formatted_healpix_dir, self.options.healpix_file)
        tile_path = "%s/%s" % (formatted_tile_dir, self.options.tile_file)

        # Check if the above files exist...
        if not os.path.exists(hpx_path):
            is_error = True
            print("Healpix file `%s` does not exist." % hpx_path)

        if not os.path.exists(tile_path):
            is_error = True
            print("Tile file `%s` does not exist." % tile_path)

        if self.options.tele not in detector_mapping:
            is_error = True
            print("Unknown telescope abbreviation: %s " % self.options.tele)

        if is_error:
            print("Exiting...")
            return 1

        print("\tLoading NSIDE 128 pixels...")
        nside128 = 128
        N128_dict = None
        with open('N128_dict.pkl', 'rb') as handle:
            N128_dict = pickle.load(handle)
        del handle

        print("\tLoading existing EBV...")
        ebv = None
        with open('ebv.pkl', 'rb') as handle:
            ebv = pickle.load(handle)

        # Get Map ID
        healpix_map_select = "SELECT id, RescaledNSIDE FROM HealpixMap WHERE GWID = '%s' and Filename = '%s'"
        healpix_map_id = int(query_db([healpix_map_select % (self.options.gw_id, self.options.healpix_file)])[0][0][0])
        healpix_map_nside = int(
            query_db([healpix_map_select % (self.options.gw_id, self.options.healpix_file)])[0][0][1])

        print("Get map pixel")
        map_pixel_select = "SELECT id, HealpixMap_id, Pixel_Index, Prob, Distmu, Distsigma, Distnorm, Mean, Stddev, Norm, N128_SkyPixel_id FROM HealpixPixel WHERE HealpixMap_id = %s;"
        q = map_pixel_select % healpix_map_id
        map_pixels = query_db([q])[0]
        print("Retrieved %s map pixels..." % len(map_pixels))

        # Initialize map pix dict for later access
        map_pixel_dict = {}
        for p in map_pixels:
            map_pixel_dict[int(p[2])] = p

        band_select = "SELECT id, Name, F99_Coefficient FROM Band WHERE `Name`='%s'"
        detector_select_by_name = "SELECT id, Name, Deg_width, Deg_height, Deg_radius, Area, MinDec, MaxDec, ST_AsText(Poly) FROM Detector WHERE Name='%s'"

        obs_tile_insert_data = []
        detectors = {}

        tele_name = detector_mapping[self.options.tele]
        detector_result = query_db([detector_select_by_name % tele_name])[0][0]

        detector_name = detector_result[1]
        detector_poly = Detector.get_detector_vertices_from_teglon_db(detector_result[8])
        detector = Detector(detector_name, detector_poly, detector_id=int(detector_result[0]))
        # detector = Detector(detector_result[1], float(detector_result[2]), float(detector_result[2]))
        # detector.id = int(detector_result[0])
        # detector.area = float(detector_result[5])

        if detector.name not in detectors:
            detectors[detector.name] = detector
        print("Processing `%s` for %s" % (tile_path, detector.name))


        # Iterate over lines of a tile
        with open(tile_path, 'r') as csvfile:
            # Read CSV lines
            csvreader = csv.reader(csvfile, delimiter=' ', skipinitialspace=True)
            for row in csvreader:

                # default PA. Overwrite below if provided
                position_angle = 0.0

                # Format 0425  - Swope, KAIT
                file_name = row[0]
                field_name = row[1]
                ra = float(row[2])
                dec = float(row[3])
                mjd = float(row[4])
                band = row[5].strip()
                exp_time = float(row[6])
                mag_lim = None
                try:
                    mag_lim = float(row[7])
                except:
                    pass

                # # Format 0425  - Thacher, Nickel
                # file_name = row[0]
                # field_name = row[1]
                # ra = float(row[2])
                # dec = float(row[3])
                # exp_time = float(row[4])
                # mjd = float(row[5])
                # band = row[6].strip()
                # mag_lim = None
                # try:
                #     mag_lim = float(row[7])
                # except:
                #     pass

                # # Format 0425 Treasure Map
                # file_name = row[0]
                # field_name = row[1]
                # ra = float(row[2])
                # dec = float(row[3])
                # mjd = float(row[4])
                # band = row[5].strip()
                # exp_time = float(row[6])
                # mag_lim = float(row[7])
                # position_angle = float(row[8])

                # # Format 0814 Swope, Nickel, Wise
                # file_name = row[0]
                # field_name = row[0]
                # ra = float(row[1])
                # dec = float(row[2])
                # mjd = float(row[3])
                # band = row[4].strip()
                # exp_time = float(row[5])
                # mag_lim = float(row[6])

                # # Format 0814 Thacher
                # file_name = row[0]
                # field_name = row[1]
                # ra = float(row[2])
                # dec = float(row[3])
                # mjd = float(row[4])
                # band = row[5].strip()
                # exp_time = float(row[6])
                # mag_lim = float(row[7])

                # Format 0814 MOSFIRE
                # file_name = row[0]
                # field_name = row[1]
                # ra = float(row[2])
                # dec = float(row[3])
                # mjd = float(row[4])
                # band = row[5].strip()
                # exp_time = float(row[6])
                # mag_lim = None

                # # Format 0814 LCOGT, KAIT
                # file_name = row[0]
                # field_name = row[1]
                # ra = float(row[2])
                # dec = float(row[3])
                # mjd = float(row[4])
                # band = row[5].strip()
                # exp_time = float(row[6])
                # mag_lim = None
                # try:
                #     mag_lim = float(row[7])
                # except:
                #     pass


                # mag_lim = None


                # file_name = row[0]
                # field_name = row[1]
                # ra = float(row[2])
                # dec = float(row[3])
                # mjd = float(row[4])
                # band = row[5].strip()
                # exp_time = float(row[6])

                # exp_time = float(row[4])
                # mjd = float(row[5])
                # band = row[6].strip()

                # For ANDICAM, the filter is in the filename...
                # if self.options.tele == "a" and band == "___":
                #     band = file_name.split(".")[1]

                # # KAIT data format
                # ra = float(row[2])
                # dec = float(row[3])
                # mjd = float(row[4])
                # band = row[5].strip()

                # exp_time = None
                # try:
                #     exp_time = float(row[6])
                #     if exp_time <= 0.0:
                #         exp_time = None
                # except:
                #     pass


                # mag_lim = None
                # try:
                #     mag_lim = float(row[7])
                # except:
                #     pass

                # Get Band_id
                band_map = band_mapping[band]
                band_results = query_db([band_select % band_map])[0][0]

                band_id = band_results[0]
                band_name = band_results[1]
                band_F99 = float(band_results[2])

                dec_unit = u.hour
                if self.options.coord_format == DEG_FORMAT:
                    dec_unit = u.deg
                c = coord.SkyCoord(ra, dec, unit=(dec_unit, u.deg))

                n128_index = hp.ang2pix(nside128, 0.5 * np.pi - c.dec.radian, c.ra.radian)  # theta, phi
                n128_id = N128_dict[n128_index]

                # t = Tile(c.ra.degree, c.dec.degree, detector.deg_width, detector.deg_height, int(healpix_map_nside))
                t = Tile(c.ra.degree, c.dec.degree, detector, int(healpix_map_nside), position_angle_deg=position_angle)

                t.field_name = field_name
                t.N128_pixel_id = n128_id
                t.N128_pixel_index = n128_index
                t.mwe = ebv[n128_index] * band_F99
                t.mjd = mjd
                t.exp_time = exp_time
                t.mag_lim = mag_lim

                # observed_tiles.append(t)
                obs_tile_insert_data.append((
                    file_name,
                    detector.id,
                    t.field_name,
                    t.ra_deg,
                    t.dec_deg,
                    "POINT(%s %s)" % (t.dec_deg, t.ra_deg - 180.0),  # Dec, RA order due to MySQL convention for lat/lon
                    t.query_polygon_string,
                    str(t.mwe),
                    t.N128_pixel_id,
                    band_id,
                    t.mjd,
                    t.exp_time,
                    t.mag_lim,
                    healpix_map_id,
                    t.position_angle_deg))

        insert_observed_tile = '''
            INSERT INTO
                ObservedTile (FileName, Detector_id, FieldName, RA, _Dec, Coord, Poly, EBV, N128_SkyPixel_id, Band_id, MJD, Exp_Time, Mag_Lim, HealpixMap_id, PositionAngle)
            VALUES (%s, %s, %s, %s, %s, ST_PointFromText(%s, 4326), ST_GEOMFROMTEXT(%s, 4326), %s, %s, %s, %s, %s, %s, %s, %s)
        '''

        print("Inserting %s tiles..." % len(obs_tile_insert_data))
        batch_insert(insert_observed_tile, obs_tile_insert_data)
        print("Done...")

        # Associate map pixels with observed tiles
        print("Building observed tile-healpix map pixel relation...")

        obs_tile_select = '''
            SELECT FileName, id, Detector_id, FieldName, RA, _Dec, Coord, Poly, EBV, N128_SkyPixel_id, Band_id, MJD, Exp_Time, Mag_Lim, HealpixMap_id, PositionAngle 
            FROM ObservedTile 
            WHERE Detector_id = %s and HealpixMap_id = %s 
        '''

        # Obtain the tile with id's so they can be used in later INSERTs
        tile_pixel_data = []
        for d_name, d in detectors.items():
            obs_tiles_for_detector = query_db([obs_tile_select % (d.id, healpix_map_id)])[0]

            print("Observed Tiles for %s: %s" % (d.name, len(obs_tiles_for_detector)))
            for otfd in obs_tiles_for_detector:

                ot_id = int(otfd[1])
                ra = float(otfd[4])
                dec = float(otfd[5])
                pa = float(otfd[15])

                # t = Tile(ra, dec, d.deg_width, d.deg_height, healpix_map_nside)
                t = Tile(ra, dec, d, healpix_map_nside, position_angle_deg=pa)

                t.id = ot_id

                for p in t.enclosed_pixel_indices:
                    tile_pixel_data.append((t.id, map_pixel_dict[p][0]))

        print("Length of tile_pixel_data: %s" % len(tile_pixel_data))

        tile_pixel_upload_csv = "%s/ObservedTiles/%s_tile_pixel_upload.csv" % \
                                (formatted_healpix_dir, detector.name.replace("/", "_"))

        # Create CSV
        try:
            t1 = time.time()
            print("Creating `%s`" % tile_pixel_upload_csv)
            with open(tile_pixel_upload_csv, 'w') as csvfile:
                csvwriter = csv.writer(csvfile)
                for data in tile_pixel_data:
                    csvwriter.writerow(data)

            t2 = time.time()
            print("\n********* start DEBUG ***********")
            print("Tile-Pixel CSV creation execution time: %s" % (t2 - t1))
            print("********* end DEBUG ***********\n")
        except Error as e:
            print("Error in creating Tile-Pixel CSV:\n")
            print(e)
            print("\nExiting")
            return 1

        # clean up
        print("freeing `tile_pixel_data`...")
        del tile_pixel_data

        print("Bulk uploading `tile_pixel_data`...")
        ot_hp_upload_sql = """LOAD DATA LOCAL INFILE '%s' 
                    INTO TABLE ObservedTile_HealpixPixel 
                    FIELDS TERMINATED BY ',' 
                    LINES TERMINATED BY '\n' 
                    (ObservedTile_id, HealpixPixel_id);"""

        success = bulk_upload(ot_hp_upload_sql % tile_pixel_upload_csv)
        if not success:
            print("\nUnsuccessful bulk upload. Exiting...")
            return 1

        try:
            print("Removing `%s`..." % tile_pixel_upload_csv)
            os.remove(tile_pixel_upload_csv)

            print("... Done")
        except Error as e:
            print("Error in file removal")
            print(e)
            print("\nExiting")
            return 1


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
