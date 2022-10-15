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

import MySQLdb as my

# region config settings

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
configFile = os.path.join(__location__, 'Settings.ini')

db_config = RawConfigParser()
db_config.read(configFile)

db_name = db_config.get('database', 'DATABASE_NAME')
db_user = db_config.get('database', 'DATABASE_USER')
db_pwd = db_config.get('database', 'DATABASE_PASSWORD')
db_host = db_config.get('database', 'DATABASE_HOST')
db_port = db_config.get('database', 'DATABASE_PORT')

# Set up dustmaps config
config["data_dir"] = "./"

cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

# endregion

# region database methods

def bulk_upload(query):
    success = False
    try:

        conn = mysql.connector.connect(user=db_user, password=db_pwd, host=db_host, port=db_port, database=db_name)
        cursor = conn.cursor()
        cursor.execute(query)
        conn.commit()
        success = True

    except Error as e:
        print("Error in uploading CSV!")
        print(e)
    finally:
        cursor.close()
        conn.close()

    return success

def query_db(query_list, commit=False):
    # query_string = ";".join(query_list)

    results = []
    try:
        chunk_size = 1e+6

        db = my.connect(host=db_host, user=db_user, passwd=db_pwd, db=db_name, port=3306)
        cursor = db.cursor()

        for q in query_list:
            cursor.execute(q)

            if commit:  # used for updates, etc
                db.commit()

            streamed_results = []
            print("fetching results...")
            while True:
                r = cursor.fetchmany(1000000)
                count = len(r)
                streamed_results += r
                size_in_mb = sys.getsizeof(streamed_results) / 1.0e+6

                print("\tfetched: %s; current length: %s; running size: %0.3f MB" % (
                count, len(streamed_results), size_in_mb))

                if not r or count < chunk_size:
                    break

        results.append(streamed_results)

    # cnx = mysql.connector.connect(user=db_user, password=db_pwd, host=db_host, port=db_port, database=db_name)
    # cursor = cnx.cursor()
    # for result in cursor.execute(query_string, multi=True):

    # 	streamed_results = []
    # 	print("fetching results...")

    # 	i = 0
    # 	while True:
    # 		i += chunk_size
    # 		print("Fetching: %s records" % i)

    # 		partial_result = result.fetchmany(chunk_size)
    # 		count = len(partial_result)
    # 		streamed_results += partial_result
    # 		size_in_mb = sys.getsizeof(streamed_results)/1.0e+6

    # 		print("\tfetched: %s; current length: %s; running size: %0.3f MB" % (count, len(streamed_results), size_in_mb))

    # 		if not partial_result or count < chunk_size:
    # 			break

    # 	results.append(streamed_results)

    except Error as e:
        print('Error:', e)
    finally:
        cursor.close()
        # cnx.close()
        db.close()

    # fake = [[[(1)]]]
    return results

def batch_query(query_list):
    return_data = []
    batch_size = 500
    ii = 0
    jj = batch_size
    kk = len(query_list)

    print("\nLength of data to query: %s" % kk)
    print("Query batch size: %s" % batch_size)
    print("Starting loop...")

    number_of_queries = len(query_list) // batch_size
    if len(query_list) % batch_size > 0:
        number_of_queries += 1

    query_num = 1
    payload = []
    while jj < kk:
        t1 = time.time()

        print("%s:%s" % (ii, jj))
        payload = query_list[ii:jj]
        return_data += query_db(payload)

        ii = jj
        jj += batch_size
        t2 = time.time()

        print("\n********* start DEBUG ***********")
        print("Query %s/%s complete - execution time: %s" % (query_num, number_of_queries, (t2 - t1)))
        print("********* end DEBUG ***********\n")

        query_num += 1

    print("Out of loop...")

    t1 = time.time()

    print("\n%s:%s" % (ii, kk))

    payload = query_list[ii:kk]
    return_data += query_db(payload)

    t2 = time.time()

    print("\n********* start DEBUG ***********")
    print("Query %s/%s complete - execution time: %s" % (query_num, number_of_queries, (t2 - t1)))
    print("********* end DEBUG ***********\n")

    return return_data

def insert_records(query, data):
    _tstart = time.time()
    success = False
    try:
        conn = mysql.connector.connect(user=db_user, password=db_pwd, host=db_host, port=db_port, database=db_name)
        cursor = conn.cursor()
        cursor.executemany(query, data)

        conn.commit()
        success = True
    except Error as e:
        print('Error:', e)
    finally:
        cursor.close()
        conn.close()

    _tend = time.time()
    print("\n********* start DEBUG ***********")
    print("insert_records execution time: %s" % (_tend - _tstart))
    print("********* end DEBUG ***********\n")
    return success

def batch_insert(insert_statement, insert_data, batch_size=50000):
    _tstart = time.time()

    i = 0
    j = batch_size
    k = len(insert_data)

    print("\nLength of data to insert: %s" % len(insert_data))
    print("Insert batch size: %s" % batch_size)
    print("Starting loop...")

    number_of_inserts = len(insert_data) // batch_size
    if len(insert_data) % batch_size > 0:
        number_of_inserts += 1

    insert_num = 1
    payload = []
    while j < k:
        t1 = time.time()

        print("%s:%s" % (i, j))
        payload = insert_data[i:j]

        if insert_records(insert_statement, payload):
            i = j
            j += batch_size
        else:
            raise ("Error inserting batch! Exiting...")

        t2 = time.time()

        print("\n********* start DEBUG ***********")
        print("INSERT %s/%s complete - execution time: %s" % (insert_num, number_of_inserts, (t2 - t1)))
        print("********* end DEBUG ***********\n")

        insert_num += 1

    print("Out of loop...")

    t1 = time.time()

    print("\n%s:%s" % (i, k))

    payload = insert_data[i:k]
    if not insert_records(insert_statement, payload):
        raise ("Error inserting batch! Exiting...")

    t2 = time.time()

    print("\n********* start DEBUG ***********")
    print("INSERT %s/%s complete - execution time: %s" % (insert_num, number_of_inserts, (t2 - t1)))
    print("********* end DEBUG ***********\n")

    _tend = time.time()

    print("\n********* start DEBUG ***********")
    print("batch_insert execution time: %s" % (_tend - _tstart))
    print("********* end DEBUG ***********\n")

# endregion

class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--gw_id', default="", type="str",
                          help='LIGO superevent name, e.g. `S190425z` ')

        # parser.add_option('--orig_res', action="store_true", default=False,
        #                   help='''Upload the healpix file at the native resolution (Default = False)''')
        #
        # parser.add_option('--load_local_file', action="store_true", default=False,
        #                   help='''Circumvent downloading map and load map from disk (Default = False)''')

        return (parser)

    def main(self):

        api_token = "jOwrpZgnXU-nLVi7r8bhL2wliAI0ZXi5IZUudQ"
        base_url = 'http://treasuremap.space/api/v0'
        #
        # # get all photometric instruments
        # instrument_url = 'instruments'
        # instrument_synopsis_request = {
        #     "api_token": api_token,
        #     "type": "photometric"
        # }
        # instrument_list_target_url = "{}/{}?{}".format(base_url, instrument_url,
        #                                                urllib.parse.urlencode(instrument_synopsis_request))
        # instrument_response = requests.get(url=instrument_list_target_url)
        # instrument_result = [json.loads(r) for r in json.loads(instrument_response.text)]
        # instruments = {ir["id"]:ir["instrument_name"] for ir in instrument_result}
        # print("instruments:")
        # for id, name in instruments.items():
        #     print("{}:{}".format(name, id))
        #
        # get details for all instruments
        # footprint_url = 'footprints'
        # instrument_detail_request = {
        #     "api_token": api_token,
        #     # "ids": "[{}]".format(",".join([str(id) for id, name in instruments.items()]))
        # }
        # instrument_detail_target_url = "{}/{}?{}".format(base_url, footprint_url,
        #                                                urllib.parse.urlencode(instrument_detail_request))
        # instrument_detail_response = requests.get(url=instrument_detail_target_url)
        # instrument_detail_result = [json.loads(r) for r in json.loads(instrument_detail_response.text)]
        #
        # instrument_footprints = {}
        # for idr in instrument_detail_result:
        #
        #     if idr["instrumentid"] not in instrument_footprints:
        #         instrument_footprints[idr["instrumentid"]] = []
        #
        #     instrument_footprints[idr["instrumentid"]].append(idr["footprint"])
        #
        # # instrument_footprints = {idr["instrumentid"]:idr["footprint"] for idr in  instrument_detail_result}
        #
        # [print(p) for p in instrument_footprints[47]]
        #
        # test = 1

        pointing_url = 'pointings'
        pointing_request = {
            "api_token": api_token,
            "status": "completed",
            # "bands": [
            #     "U", "B", "V", "g", "r", "i"
            # ],
            "graceid": "GW190425"
        }
        pointing_target_url = "{}/{}?{}".format(base_url, pointing_url,
                                                         urllib.parse.urlencode(pointing_request))
        pointing_response = requests.get(url=pointing_target_url)
        pointing_result = [json.loads(r) for r in json.loads(pointing_response.text)]
        test = 1




        # t = Tile(359.85, 0.0, 0.49493333333333334, 0.4968666666666667, 1024)
        # # test1 = t.query_polygon_string
        # # test = 1
        #
        # t2 = Tile(359.8, 0.0, 0.49493333333333334, 0.4968666666666667, 1024)
        #
        # fig = plt.figure(figsize=(10, 10), dpi=800)
        # ax = fig.add_subplot(111)
        # m = Basemap(projection='ortho', lon_0=0.0, lat_0=0.0)
        # m.drawparallels(np.arange(-90., 120., 30.))
        # m.drawmeridians(np.arange(0., 420., 60.))
        #
        # t.plot(m, ax, edgecolor='blue', facecolor='None', linewidth=0.25, alpha=1.0,
        #        zorder=9900)
        # t2.plot(m, ax, edgecolor='red', facecolor='None', linewidth=0.25, alpha=1.0,
        #        zorder=9900)
        #
        # fig.savefig("test_tile.png", bbox_inches='tight')  # ,dpi=840
        # plt.close('all')


        # json_data = {
        #     "graceid": self.options.gw_id,
        #     "api_token": tm_api_key,
        #     "pointings": self.pointings
        # }

        # instruments_ids = {
        #     "ZTF":47,
        #
        #
        # }
        #
        # json_data = {
        #     "graceid": self.options.gw_id,
        #     "api_token": tm_api_key,
        #     "pointings": self.pointings
        # }


        # pass


if __name__ == "__main__":
    useagestring = """python GetTreasureMapData.py [options]

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
    print("Teglon `GetTreasureMapData` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")
