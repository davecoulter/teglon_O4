import argparse
import time
from configparser import RawConfigParser
import urllib.request
import requests
import json

from web.src.utilities.Database_Helpers import *
from web.src.objects.Detector import *

class Teglon_Detector_Util:
    def __int__(self):
        pass

    #
    def add_tm_detector(self, tm_detector_id, min_dec=-90, max_dec=90.0, detector_geometry=None,
                        detector_width=None, detector_height=None, detector_radius=None):

        configFile = './Settings.ini'
        config = RawConfigParser()
        config.read(configFile)

        treasure_map_detectors = None
        tm_api_token = None
        tm_base_url = None

        # Sanity
        is_error = False
        try:
            tm_api_token = config.get('treasuremap', 'TM_API_TOKEN')
            tm_base_url = config.get('treasuremap', 'TM_ENDPOINT')
        except:
            print("Error! No Treasure Map keys (`TM_API_TOKEN`, `TM_ENDPOINT`) in the Settings.ini file!")
            is_error = True

        if not any(tm_api_token) or not any(tm_base_url):
            print("Can't load TM detectors with API keys or API endpoints!")
            is_error = True

        if tm_detector_id is None:
            print("Error! You must provide a correct TM Detector Id!")
            is_error = True

        if is_error:
            print("Exiting...")
            return 1

        instrument_url = 'instruments'
        instrument_synopsis_request = {
            "api_token": tm_api_token,
            "id": tm_detector_id
        }
        instrument_list_target_url = "{}/{}?{}".format(tm_base_url, instrument_url,
                                                       urllib.parse.urlencode(instrument_synopsis_request))
        instrument_response = requests.get(url=instrument_list_target_url)
        instrument_list = json.loads(instrument_response.text)
        instrument_dict = None

        if not any(instrument_list):
            print("TM ID `%s` did not return any results. Exiting..." % tm_detector_id)
            return 1
        else:
            instrument_dict = instrument_list[0]

        instrument_full_name = instrument_dict["instrument_name"]
        instrument_nickname = instrument_dict["nickname"]

        # Check if detector is already in the db
        detector_select = '''
            SELECT id, LOWER(`Name`)
            FROM Detector 
            WHERE
                LOWER(`Name`) = '%s' OR 
                LOWER(`Name`)='%s';
        '''
        existing_detector_result = query_db([detector_select % (instrument_full_name, instrument_nickname)])[0]
        if any(existing_detector_result):
            print("\n*****\nDetector: `%s` is already in the database! Exiting...\n****" % instrument_full_name)
            return 1

        print("\n*****\nDetector: `%s` is new to Teglon. Adding...\n****" % instrument_full_name)

        footprint_url = 'footprints'
        instrument_detail_request = {
            "api_token": tm_api_token,
            "id": tm_detector_id
        }
        instrument_detail_target_url = "{}/{}?{}".format(tm_base_url, footprint_url,
                                                         urllib.parse.urlencode(
                                                             instrument_detail_request))
        instrument_detail_response = requests.get(url=instrument_detail_target_url)
        instrument_detail_dict = json.loads(instrument_detail_response.text)[0]
        instrument_vertices = Detector.\
            get_detector_vertices_from_treasuremap_footprint_string([instrument_detail_dict["footprint"]])
        tm_detector = Detector(instrument_full_name, instrument_vertices)

        is_detector_rectangular = False
        is_detector_circular = False
        if detector_geometry == 'rectangle' and detector_width is not None and detector_height is not None:
            tm_detector.deg_width = detector_width
            tm_detector.deg_height = detector_height
            is_detector_rectangular = True
        elif detector_geometry == 'circle' and detector_radius is not None:
            tm_detector.deg_radius = detector_radius
            is_detector_circular = True

        if is_detector_rectangular:
            insert_tm = '''INSERT INTO Detector (Name, Poly, Area, TM_id, MinDec, MaxDec, Deg_width, Deg_height) VALUES 
                                                        ('%s', ST_GEOMFROMTEXT('%s', 4326), %s, %s, %s, %s, %s, %s);'''
            q = insert_tm % (
                tm_detector.name, tm_detector.query_polygon_string, tm_detector.area, tm_detector_id, min_dec, max_dec,
                tm_detector.deg_width, tm_detector.deg_height)
            query_db([q], commit=True)
        elif is_detector_circular:
            insert_tm = '''INSERT INTO Detector (Name, Poly, Area, TM_id, MinDec, MaxDec, Deg_radius) VALUES 
                                                        ('%s', ST_GEOMFROMTEXT('%s', 4326), %s, %s, %s, %s, %s);'''
            q = insert_tm % (
                tm_detector.name, tm_detector.query_polygon_string, tm_detector.area, tm_detector_id, min_dec, max_dec,
                tm_detector.deg_radius)
            query_db([q], commit=True)
        else:
            insert_tm = '''INSERT INTO Detector (Name, Poly, Area, TM_id, MinDec, MaxDec) VALUES 
                                                        ('%s', ST_GEOMFROMTEXT('%s', 4326), %s, %s, %s, %s);'''
            q = insert_tm % (
                tm_detector.name, tm_detector.query_polygon_string, tm_detector.area, tm_detector_id, min_dec, max_dec)
            query_db([q], commit=True)

        print("Detector `%s` added!" % tm_detector.name)


if __name__ == "__main__":
    start = time.time()

    detector_util = Teglon_Detector_Util()

    parser = argparse.ArgumentParser()
    parser.add_argument('--tm_detector_id', type=int, help='Treasure Map detector ID')
    parser.add_argument('--min_dec', type=float, default=-90.0, help='Min dec in deg. Default = -90.0')
    parser.add_argument('--max_dec', type=float, default=90.0, help='Max dec in deg. Default = +90.0')
    parser.add_argument('--detector_geometry', type=str, default="rectangle", help='''Specify shape of detector. Valid 
    options: `rectangle` or `circle`. If geometry is an arbitrary polygon, leave blank.''')
    parser.add_argument('--detector_width', type=float, help='''Only used if `detector_geometry`==`rectangle`; width in 
    deg. Ignored if `detector_geometry` is not `rectangle`.''')
    parser.add_argument('--detector_height', type=float, help='''Only used if `detector_geometry`==`rectangle`; height in 
    deg. Ignored if `detector_geometry` is not `rectangle`.''')
    parser.add_argument('--detector_radius', type=float, help='''Only used if `detector_geometry`==`circle`; radius in 
            deg. Ignored if `detector_geometry` is not `circle`.''')

    args = parser.parse_args()
    detector_util.add_tm_detector(**vars(args))

    end = time.time()
    duration = (end - start)
    print("\n********* start DEBUG ***********")
    print("Teglon `add_detector` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")