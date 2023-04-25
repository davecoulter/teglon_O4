import matplotlib
matplotlib.use("Agg")

import matplotlib as mpl
import matplotlib.pyplot as plt
# from mpl_toolkits.basemap import Basemap

from matplotlib.patches import Polygon
from matplotlib.pyplot import cm
from matplotlib.patches import CirclePolygon
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

# import sys
# sys.path.append('../')

import os
from configparser import RawConfigParser
import multiprocessing as mp
import mysql.connector
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
from dustmaps.config import config as dustmaps_config
from dustmaps.sfd import SFDQuery

import shapely as ss
from shapely.ops import transform as shapely_transform
from shapely.geometry import Point
from shapely.ops import linemerge, unary_union, polygonize, split
from shapely import geometry

import healpy as hp
from ligo.skymap import distance

# import sys
# sys.path.append('../objects/')
# sys.path.append('./')


# Old
# from src.utilities.HEALPix_Helpers import *
# from src.utilities.Database_Helpers import *
#
# from src.objects.Detector import *
# from src.objects.Tile import *
# from src.objects.SQL_Polygon import *
# from src.objects.Pixel_Element import *
# from src.objects.Completeness_Objects import *

# NEW
from web.src.utilities.HEALPix_Helpers import *
from web.src.utilities.Database_Helpers import *

from web.src.objects.Detector import *
from web.src.objects.Tile import *
from web.src.objects.SQL_Polygon import *
from web.src.objects.Pixel_Element import *
from web.src.objects.Completeness_Objects import *

# from mpl_toolkits.basemap import Basemap

import multiprocessing as mp

class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--is_debug', action="store_true", default=False,
                          help='Output debug information during initialization (Default = False)')

        parser.add_option('--build_skydistances', action="store_true", default=False,
                          help='Build sky distance table (Default = False)')

        parser.add_option('--build_skypixels', action="store_true", default=False,
                          help='Build sky pixel table; HACK, remove table FK on parent pixel id first, then restore '
                               'after update (Default = False)')

        parser.add_option('--build_detectors', action="store_true", default=False,
                          help='Build detector table (Default = False)')

        parser.add_option('--build_TM_detectors', action="store_true", default=False,
                          help='Build detectors from Treasure Map (Default = False). Requires TM tokens in '
                               'settings.ini file')

        parser.add_option('--build_bands', action="store_true", default=False,
                          help='Build band table (Default = False)')

        parser.add_option('--build_MWE', action="store_true", default=False,
                          help='Build EBV table (Default = False)')

        parser.add_option('--build_galaxy_skypixel_associations', action="store_true", default=False,
                          help='Register galaxies to sky pixels (Default = False)')

        parser.add_option('--build_completeness', action="store_true", default=False,
                          help='Build completeness table (Default = False)')

        parser.add_option('--compose_completeness', action="store_true", default=False,
                          help='Compose completness records into function input (Default = False)')

        parser.add_option('--build_static_grids', action="store_true", default=False,
                          help='Build Swope and Thacher grids (Default = False)')

        return (parser)

    def generate_static_grid_inserts(self, detector_name_list, ebv_dict, skypixel_dict):

        detector_prefixes = {
            "SWOPE": "S",
            "THACHER": "T",
            "T80S_T80S-Cam": "T80_"
        }

        # Sanity - check if this is a supported detector, else return
        for detector_name in detector_name_list:
            if detector_name not in detector_prefixes:
                print("Unsupported detector! `%s`" % detector_name)
                return 1

        static_grid_nside = 128
        select_detect = '''
                        SELECT
                            Name, ST_AsText(Poly), Deg_width,
                            Deg_height, id, MinDec, MaxDec
                        FROM Detector
                        WHERE Name='%s';
                    '''
        insert_static_tile = '''INSERT INTO StaticTile (Detector_id, FieldName, RA, _Dec, Coord, Poly, EBV, N128_SkyPixel_id) 
                    VALUES (%s, %s, %s, %s, ST_PointFromText(%s, 4326), ST_GEOMFROMTEXT(%s, 4326), %s, %s)'''


        tiles = []
        for detector_name in detector_name_list:
            print("Building `%s` coordinate grid ..." % detector_name)
            detector_result = query_db([select_detect % detector_name])[0][0]
            detector_vertices = Detector.get_detector_vertices_from_teglon_db(detector_result[1])
            detector_id = int(detector_result[4])
            detector_min_dec = float(detector_result[5])
            detector_max_dec = float(detector_result[6])
            detector = Detector(detector_name=detector_name, detector_vertex_list_collection=detector_vertices,
                             detector_width_deg=detector_result[2], detector_height_deg=detector_result[3],
                             detector_id=detector_id)

            t1 = time.time()
            allsky_detector_coords = Cartographer.generate_all_sky_coords(detector)
            t2 = time.time()
            print("\n********* start DEBUG ***********")
            print("`%s` Coordinate grid creation - execution time: %s" % (detector_name, (t2 - t1)))
            print("********* end DEBUG ***********\n")

            print("Building `%s` tiles ..." % detector_name)
            t1 = time.time()
            for i, c in enumerate(allsky_detector_coords):
                # only include pointings that the telescope can reach
                if c[1] >= detector_min_dec and c[1] <= detector_max_dec:
                    t = Tile(c[0], c[1], detector, static_grid_nside)
                    t.field_name = "{detector_prefix}{field_number}".format(detector_prefix=detector_prefixes[detector_name],
                                                                            field_number=str(i).zfill(6))
                    # Sanity - due to discretization err, sometimes (0.5 * np.pi - t.dec_rad) < 0. If this is the case, set to 0
                    # theta = 0.5 * np.pi - t.dec_rad
                    # if theta < 0:
                    #     theta = 0.0
                    n128_index = hp.ang2pix(static_grid_nside, 0.5 * np.pi - t.dec_rad, t.ra_rad)  # theta, phi
                    t.mwe = ebv_dict[n128_index]
                    # HACK: DC - this is overloaded... the `id` is for the Tile id, but the N128_Pixel_id is put here for convenience
                    t.id = skypixel_dict[static_grid_nside][n128_index].id
                    tiles.append(t)

            t2 = time.time()
            print("\n********* start DEBUG ***********")
            print("`%s` tile creation - execution time: %s" % (detector_name, (t2 - t1)))
            print("********* end DEBUG ***********\n")

        print("Building INSERTS for all static tiles ...")
        static_tile_data = []
        for t in tiles:
            static_tile_data.append((t.detector.id,
                                     t.field_name,
                                     t.ra_deg,
                                     t.dec_deg,
                                     "POINT(%s %s)" % (t.dec_deg, t.ra_deg - 180.0),
                                     # Dec, RA order due to MySQL convention for lat/lon
                                     t.query_polygon_string,
                                     str(t.mwe),
                                     t.id) # HACK: DC - this is overloaded... the `id` is for the Tile id, but the N128_Pixel_id is put here for convenience
                                    )

        print("INSERTing all static tiles ...")
        batch_insert(insert_static_tile, static_tile_data)

    def main(self):

        configFile = './Settings.ini'
        config = RawConfigParser()
        config.read(configFile)

        utilities_base_dir = "./web/src/utilities"
        pickle_output_dir = "%s/pickles/" % utilities_base_dir
        if not os.path.exists(pickle_output_dir):
            os.makedirs(pickle_output_dir)
        # Set up dustmaps config
        # DC - Why is this directory set this way?
        dustmaps_config["data_dir"] = utilities_base_dir




        print("Initializing database named: %s" % db_name)


        # Generate all pixel indices
        h = 0.7
        cosmo = LambdaCDM(H0=100 * h, Om0=0.27, Ode0=0.73)

        nside2 = 2
        nside2_npix = hp.nside2npix(nside2)
        frac_sky_nside2 = (1/nside2_npix)

        nside4 = 4
        nside4_npix = hp.nside2npix(nside4)
        frac_sky_nside4 = (1/nside4_npix)

        nside8 = 8
        nside8_npix = hp.nside2npix(nside8)
        frac_sky_nside8 = (1/nside8_npix)

        nside16 = 16
        nside16_npix = hp.nside2npix(nside16)
        frac_sky_nside16 = (1/nside16_npix)

        nside32 = 32
        nside32_npix = hp.nside2npix(nside32)
        frac_sky_nside32 = (1/nside32_npix)

        nside64 = 64
        nside64_npix = hp.nside2npix(nside64)
        frac_sky_nside64 = (1/nside64_npix)

        nside128 = 128
        nside128_npix = hp.nside2npix(nside128)
        frac_sky_nside128 = (1/nside128_npix)

        nside_npix = [nside2_npix, nside4_npix, nside8_npix, nside16_npix, nside32_npix, nside64_npix, nside128_npix]
        nsides = [nside2, nside4, nside8, nside16, nside32, nside64, nside128]
        frac = [frac_sky_nside2, frac_sky_nside4, frac_sky_nside8, frac_sky_nside16, frac_sky_nside32, frac_sky_nside64, frac_sky_nside128]

        distance_at_resolution_change = [45,125,200,400,700,900,1220] # Mpc
        steps_in_resolution_bin = [6,24,19,37,38,13,8]


        if self.options.is_debug:
            print("NSIDE: %s, npix: %s, frac/pix: %0.3E" % (nside2, nside2_npix, frac_sky_nside2))
            print("NSIDE: %s, npix: %s, frac/pix: %0.3E" % (nside4, nside4_npix, frac_sky_nside4))
            print("NSIDE: %s, npix: %s, frac/pix: %0.3E" % (nside8, nside8_npix, frac_sky_nside8))
            print("NSIDE: %s, npix: %s, frac/pix: %0.3E" % (nside16, nside16_npix, frac_sky_nside16))
            print("NSIDE: %s, npix: %s, frac/pix: %0.3E" % (nside32, nside32_npix, frac_sky_nside32))
            print("NSIDE: %s, npix: %s, frac/pix: %0.3E" % (nside64, nside64_npix, frac_sky_nside64))
            print("NSIDE: %s, npix: %s, frac/pix: %0.3E" % (nside128, nside128_npix, frac_sky_nside128))
            print("\nDistances at resolution change: %s [Mpc]" % distance_at_resolution_change)

        ## Build Sky Distance Objects ##
        # start_comoving_volume is Quantity [Mpc^3]
        # end_comoving_volume is Quantity [Mpc^3]
        # start_distance is Quantity [Mpc]
        def compute_sky_distance(start_comoving_volume, end_comoving_volume, steps, start_distance, sky_fraction, nside):
            sky_distances = []
            volume_steps = np.linspace(start_comoving_volume, end_comoving_volume, steps)
            d1 = start_distance
            d2 = d1
            for i in range(len(volume_steps)):
                if i == 0:
                    continue
                z = z_at_value(cosmo.comoving_volume, volume_steps[i], zmin=2e-5)
                d2 = cosmo.luminosity_distance(z)
                partial_vol = (volume_steps[i]-volume_steps[i-1])*sky_fraction
                sky_distances.append(Sky_Distance(d1,d2,partial_vol,nside))
                d1 = d2
            return sky_distances

        ################################################

        if not self.options.build_skydistances:
            print("Skipping Sky Distances...")
        else:
            ## Build Sky Distance Objects ##
            print("\n\nBuilding Sky Distances...")

            sky_distances = []
            for i in range(len(distance_at_resolution_change)):

                start_vol = 0.0*u.Mpc**3
                end_vol = cosmo.comoving_volume((Distance(distance_at_resolution_change[i], u.Mpc)).z)
                start_distance = 0.0*u.Mpc

                if i > 0:
                    start_vol = cosmo.comoving_volume((Distance(distance_at_resolution_change[i-1], u.Mpc)).z)
                    start_distance = sky_distances[-1].D2*u.Mpc

                sky_distances += compute_sky_distance(start_vol, end_vol, steps_in_resolution_bin[i], start_distance, frac[i], nsides[i])

            if self.options.is_debug:
                print("\nGenerated distances:\n")
                for d in sky_distances:
                    print(d)

            sky_distance_insert = "INSERT INTO SkyDistance (D1, D2, dCoV, Sch_L10, NSIDE) VALUES (%s,%s,%s,%s,%s)"
            sky_distance_data = []

            for d in sky_distances:
                sky_distance_data.append((d.D1, d.D2, d.dCoV, d.Sch_L10, d.NSIDE))

            print("\nInserting %s sky distances..." % len(sky_distance_data))
            if insert_records(sky_distance_insert, sky_distance_data):
                print("Success!")
            else:
                raise("Error with INSERT! Exiting...")
            print("...Done")

        ################################################

        sky_pixels = None
        if not self.options.build_skypixels:
            print("Skipping Sky Pixels...")
            print("\tLoading existing pixels...")
            if os.path.exists(pickle_output_dir + "sky_pixels.pkl"):
                with open(pickle_output_dir + "sky_pixels.pkl", 'rb') as handle:
                    sky_pixels = pickle.load(handle)
            else:
                print('sky_pixels.pkl does not exist!')
        else:
            ## Build Sky Pixel Objects ##
            print("\n\nBuilding Sky Pixels...")
            sky_pixels = OrderedDict() # Ordered Dict because it will matter the order once we start iterating during the INSERT

            for i, nside in enumerate(nsides):

                print("Initializing NSIDE=%s pix..." % nside)
                sky_pixels[nside] = {}

                for j in range(nside_npix[i]):
                    pe = Pixel_Element(j, nside, 0.0)
                    pe.query_polygon_string # initialize
                    sky_pixels[nside][j] = pe

                if i > 0:
                    for j, (pi, pe) in enumerate(sky_pixels[nside].items()):
                        # Get Parent Pixel
                        theta, phi = hp.pix2ang(pe.nside, pe.index)
                        parent_index = hp.ang2pix(nsides[i-1], theta, phi)

                        if not parent_index in sky_pixels[nsides[i-1]]:
                            raise("Orphaned NSIDE=%s pixel! Orphaned index, theta, phi: (%s, %s, %s)" % (nside, pe.index, theta, phi))
                        else:
                            parent_pixel = sky_pixels[nsides[i-1]][parent_index]
                            sky_pixels[nside][j].parent_pixel = parent_pixel
                print("... created %s NSIDE=%s pix." % (len(sky_pixels[nside]), nside))

            sky_pixel_insert = "INSERT INTO SkyPixel (RA, _Dec, Coord, Poly, NSIDE, Pixel_Index, Parent_Pixel_id) VALUES (%s,%s,ST_PointFromText(%s, 4326),ST_GEOMFROMTEXT(%s, 4326),%s,%s,%s);"
            sky_pixel_select = "SELECT id, Pixel_Index FROM SkyPixel WHERE NSIDE = %s;"
            sky_pixel_select2 = "SELECT id, NSIDE, Pixel_Index FROM SkyPixel;"

            for i, (nside, pixel_dict) in enumerate(sky_pixels.items()):

                sky_pixel_data = []
                if i == 0: # NSIDE=2 has no parents to associate
                    for pi, pe in pixel_dict.items():

                        sky_pixel_data.append((pe.coord.ra.degree, pe.coord.dec.degree,
                        "POINT(%s %s)" % (pe.coord.dec.degree, pe.coord.ra.degree - 180.0), # Dec, RA order due to MySQL convention for lat/lon
                        pe.query_polygon_string, pe.nside, pe.index, None))

                else: # Associate parents

                    parent_ids = query_db([sky_pixel_select % (nsides[i-1])])

                    for row in parent_ids[0]:
                        sky_pixels[nsides[i-1]][int(row[1])].id = int(row[0]) # Set pixel id in parent...

                    for pi, pe in pixel_dict.items():

                        parent_id = sky_pixels[nsides[i-1]][pe.parent_pixel.index].id

                        sky_pixel_data.append((pe.coord.ra.degree, pe.coord.dec.degree,
                        "POINT(%s %s)" % (pe.coord.dec.degree, pe.coord.ra.degree - 180.0), # Dec, RA order due to MySQL convention for lat/lon
                        pe.query_polygon_string, pe.nside, pe.index, parent_id))

                print("\nInserting %s sky pixels (NSIDE=%s)..." % (len(sky_pixel_data), nside))
                batch_insert(sky_pixel_insert, sky_pixel_data)

            # # Get the n128 pixel ids
            # n128_pixel_tuples = query_db([sky_pixel_select % nside128])
            # for tup in n128_pixel_tuples[0]:
            #     sky_pixels[nside128][int(tup[1])].id = int(tup[0])
            # Get all pixel ids
            pixel_tuples = query_db([sky_pixel_select2])[0]
            for tup in pixel_tuples:
                db_id = int(tup[0])
                p_nside = int(tup[1])
                p_index = int(tup[2])

                sky_pixels[p_nside][p_index].id = db_id

            with open(pickle_output_dir + "sky_pixels.pkl", 'wb') as handle:
                pickle.dump(sky_pixels, handle, protocol=pickle.HIGHEST_PROTOCOL)

        ################################################

        # region manual detector definitions
        swope_deg_width = 4096*0.435/3600.
        swope_deg_height = 4112*0.435/3600.
        swope_max_dec = 30.0
        swope_min_dec = -90.0

        andicam_deg_width = 1024*0.371/3600.
        andicam_deg_height = 1024*0.371/3600.
        andicam_max_dec = 30.0
        andicam_min_dec = -90.0

        thacher_deg_width = 2048*0.609/3600.
        thacher_deg_height = 2048*0.609/3600.
        thacher_max_dec = 90.0
        thacher_min_dec = -30.0

        nickel_deg_width = 2048*0.368/3600.
        nickel_deg_height = 2048*0.368/3600.
        nickel_max_dec = 66.75
        nickel_min_dec = -30.0

        mosfire_deg_width = 0.1023333333
        mosfire_deg_height = 0.1023333333
        mosfire_max_dec = 90.0
        mosfire_min_dec = -70.0

        kait_deg_width = 0.1133333333
        kait_deg_height = 0.1133333333
        kait_max_dec = 90.0
        kait_min_dec = -30.0

        sinistro_deg_width = 0.4416666667
        sinistro_deg_height = 0.4416666667
        sinistro_max_dec = 90.0
        sinistro_min_dec = -90.0

        _2df_deg_radius = 1.05
        _2df_max_dec = 30.0
        _2df_min_dec = -90.0

        WISE_deg_width = 0.94833333333
        WISE_deg_height = 0.94833333333
        WISE_max_dec = 90.0
        WISE_min_dec = -30.0
        # endregion
        if not self.options.build_detectors:
            print("Skipping Detectors...")
        else:
            ## Build Detector Objects ##
            print("\n\nBuilding Detectors...")

            detector_data = [
                ("SWOPE", swope_deg_width, swope_deg_height, None, swope_deg_width*swope_deg_height, swope_min_dec, swope_max_dec),
                ("ANDICAM", andicam_deg_width, andicam_deg_height, None, andicam_deg_width*andicam_deg_height, andicam_min_dec, andicam_max_dec),
                ("THACHER", thacher_deg_width, thacher_deg_height, None, thacher_deg_width*thacher_deg_height, thacher_min_dec, thacher_max_dec),
                ("NICKEL", nickel_deg_width, nickel_deg_height, None, nickel_deg_width*nickel_deg_height, nickel_min_dec, nickel_max_dec),
                ("MOSFIRE", mosfire_deg_width, mosfire_deg_height, None, mosfire_deg_width * mosfire_deg_height, mosfire_min_dec, mosfire_max_dec),
                ("KAIT", kait_deg_width, kait_deg_height, None, kait_deg_width * kait_deg_height, kait_min_dec, kait_max_dec),
                ("SINISTRO", sinistro_deg_width, sinistro_deg_height, None, sinistro_deg_width * sinistro_deg_height, sinistro_min_dec, sinistro_max_dec),
                ("2dF", None, None, _2df_deg_radius, np.pi*_2df_deg_radius**2, _2df_min_dec, _2df_max_dec),
                ("WISE", WISE_deg_width, WISE_deg_height, None, WISE_deg_width * WISE_deg_height, WISE_min_dec, WISE_max_dec)
            ]

            detector_insert = "INSERT INTO Detector (Name, Deg_width, Deg_height, Deg_radius, Area, MinDec, MaxDec) VALUES (%s, %s, %s, %s, %s, %s, %s);"
            print("\nInserting %s detectors..." % len(detector_data))
            if insert_records(detector_insert, detector_data):
                print("Success!")
            else:
                raise("Error with INSERT! Exiting...")

            # Update base detectors with their polygon representation
            detector_select = '''
                SELECT id, Name, Deg_width, Deg_height, Deg_radius FROM Detector;
            '''
            detector_update = "UPDATE Detector SET Poly = ST_GEOMFROMTEXT('%s', 4326) WHERE id = %s"
            detectors_update_objs = []
            detector_result = query_db([detector_select])[0]
            for d in detector_result:
                d_id = int(d[0])
                d_name = d[1]
                d_width = d[2]
                d_height = d[3]
                d_radius = d[4]

                # Sanity
                if d_width is None and d_height is None and d_radius is None:
                    # skip these -- these are not "known detectors" -- these are Treasuremap Detectors
                    continue

                is_rectangular = d[2] is not None
                if is_rectangular:
                    detector_deg_width = float(d_width)
                    detector_deg_height = float(d_height)
                    detector_vertex_collection = Detector.get_detector_vertices_rectangular(detector_deg_width,
                                                                                            detector_deg_height)
                else:
                    detector_radius = float(d_radius)
                    detector_vertex_collection = Detector.get_detector_vertices_circular(detector_radius)

                output_detector = Detector(d_name, detector_vertex_collection)
                detectors_update_objs.append((d_id, output_detector.query_polygon_string))

            for d in detectors_update_objs:
                q = detector_update % (d[1], d[0])
                query_db([q], commit=True)

            if self.options.build_TM_detectors:

                def get_skycoord_from_tm_point(tm_point_string):
                    ra_dec_list = tm_point_string.split("(")[1].split(")")[0].split(" ")

                    ra = float(ra_dec_list[0])
                    dec = float(ra_dec_list[1])

                    return ra, dec

                treasure_map_detectors = None
                tm_api_token = config.get('treasuremap', 'TM_API_TOKEN')
                tm_base_url = config.get('treasuremap', 'TM_ENDPOINT')

                # Sanity - only proceed if these values are present
                if not (tm_api_token == '' or tm_api_token is None or tm_base_url == '' or tm_base_url is None):
                    # get all photometric instruments
                    instrument_url = 'instruments'
                    instrument_synopsis_request = {
                        "api_token": tm_api_token,
                        "type": "photometric"
                    }
                    instrument_list_target_url = "{}/{}?{}".format(tm_base_url, instrument_url,
                                                                   urllib.parse.urlencode(instrument_synopsis_request))
                    instrument_response = requests.get(url=instrument_list_target_url)
                    instrument_result = [json.loads(r) for r in json.loads(instrument_response.text)]
                    instrument_names_full = {ir["id"]: (ir["instrument_name"], ir['nickname']) for ir in
                                             instrument_result}
                    instrument_names = {}
                    for inst_id, inst_name_tuple in instrument_names_full.items():
                        full_name = inst_name_tuple[0]
                        short_name = inst_name_tuple[1]

                        # Sanity - we already have some of these detectors loaded (e.g. Swope, Sinistro, etc).
                        # Omit duplicate detectors
                        existing_detector_names = [d[0].strip().lower() for d in detector_data]
                        if full_name.strip().lower() in existing_detector_names or short_name.strip().lower() in existing_detector_names:
                            print("\n\t**** Skipping `%s`. Detector already in db! ****\n" % full_name)
                            continue

                        if short_name.strip() == '':
                            instrument_names[inst_id] = full_name
                        else:
                            if len(full_name) > 20:
                                instrument_names[inst_id] = short_name
                            else:
                                instrument_names[inst_id] = full_name

                    # get all detector footprints
                    footprint_url = 'footprints'
                    instrument_detail_request = {
                        "api_token": tm_api_token
                    }
                    instrument_detail_target_url = "{}/{}?{}".format(tm_base_url, footprint_url,
                                                                     urllib.parse.urlencode(
                                                                         instrument_detail_request))
                    instrument_detail_response = requests.get(url=instrument_detail_target_url)
                    instrument_detail_result = [json.loads(r) for r in
                                                json.loads(instrument_detail_response.text)]

                    instrument_footprints = {}
                    for idr in instrument_detail_result:
                        if idr["instrumentid"] not in instrument_footprints:
                            instrument_footprints[idr["instrumentid"]] = []
                        instrument_footprints[idr["instrumentid"]].append(idr["footprint"])

                    instrument_vertices = {}
                    for inst_id, inst_fp in instrument_footprints.items():
                        vertext_collection_list = Detector.get_detector_vertices_from_treasuremap_footprint_string(
                            inst_fp)
                        instrument_vertices[inst_id] = vertext_collection_list

                    treasure_map_detectors = {}
                    for inst_id, inst_vert in instrument_vertices.items():

                        # Sanity -- TM returns instruments that don't exist on the detail pages...
                        # ex: instrument id = 21 is returned from the footprint query, however, navigating to:
                        #   http://gwtmlb-1777941625.us-east-2.elb.amazonaws.com/instrument_info?id=21
                        # results in a 404. Omit these results.
                        if inst_id not in instrument_names:
                            continue

                        # Don't handle Swift XRT yet...
                        xrt_id = 13
                        if inst_id == xrt_id:
                            continue

                        name = instrument_names[inst_id]
                        tm_detector = Detector(name, inst_vert)
                        treasure_map_detectors[inst_id] = tm_detector

                    insert_tm = '''INSERT INTO Detector (Name, Poly, Area, TM_id) VALUES 
                                    ('%s', ST_GEOMFROMTEXT('%s', 4326), %s, %s);'''
                    for inst_id, tm_detect in treasure_map_detectors.items():
                        q = insert_tm % (
                        tm_detect.name, tm_detect.query_polygon_string, tm_detect.area, inst_id)
                        query_db([q], commit=True)

                    # HACK: DC - TM doesn't store declination limits natively, so I am inserting this here since it will
                    # be used to construct a static grid for Charlie...
                    update_T80_tm = '''
                        WITH T80S_Cam AS (
                            SELECT id FROM Detector WHERE TM_id = 80
                        )
                        UPDATE 
                            Detector 
                        SET
                            MaxDec = 10.0,
                            Deg_width = 1.4216,
                            Deg_height = 1.4238
                        WHERE id = (SELECT id FROM T80S_Cam);
                    '''
                    query_db([update_T80_tm], commit=True)
                else:
                    print("Can't load TM detectors with API keys or API endpoints!!")
            print("...Done")

        ################################################

        if not self.options.build_bands:
            print("Skipping Bands...")
        else:
            ## Build Band Objects ##
            print("\n\nBuilding Photometric Bands...")
            # Schlafly & Finkbeiner, 2018: https://arxiv.org/pdf/1012.4804.pdf
            # Name, Effective_Wavelength, F99 coefficient (Rv = 3.1); See Table 6.
            band_data = [
                ("SDSS u", 3586.8, 4.239),
                ("SDSS g", 4716.7, 3.303),
                ("SDSS r", 6165.1, 2.285),
                ("SDSS i", 7475.9, 1.698),
                ("SDSS z", 8922.9, 1.263),
                ("Landolt B", 4329.0,  3.626),
                ("Landolt V", 5421.7, 2.742),
                ("Landolt R", 6427.8,  2.169),
                ("Landolt I", 8048.4, 1.505),
                ("UKIRT J", 12482.9, 0.709),
                ("UKIRT H", 16588.4, 0.449),
                ("UKIRT K", 21897.7, 0.302),
                ("Clear", 5618.41, 0.91)
            ]

            band_insert = "INSERT INTO Band (Name, Effective_Wavelength, F99_Coefficient) VALUES (%s, %s, %s);"

            print("\nInserting %s bands..." % len(band_data))
            if insert_records(band_insert, band_data):
                print("Success!")
            else:
                raise("Error with INSERT! Exiting...")
            print("...Done")

        ################################################

        ebv = None
        if not self.options.build_MWE:
            print("Skipping MWE...")
            print("\tLoading existing EBV...")
            if os.path.exists(pickle_output_dir + 'ebv.pkl'):
                with open(pickle_output_dir + 'ebv.pkl', 'rb') as handle:
                    ebv = pickle.load(handle)
            else:
                print('ebv.pkl does not exist!')
        else:
            ## Build E(B-V) Map ##
            print("\n\nBuilding MWE...")

            nside128_RA = []
            nside128_DEC = []
            nside128_ids = []

            t1 = time.time()

            for n128, pe in sky_pixels[nside128].items():
                nside128_RA.append(pe.coord.ra.degree)
                nside128_DEC.append(pe.coord.dec.degree)
                nside128_ids.append(pe.id)

            c = coord.SkyCoord(nside128_RA, nside128_DEC, unit=(u.deg,u.deg))
            sfd = SFDQuery()
            ebv = sfd(c)

            t2 = time.time()

            print("\n********* start DEBUG ***********")
            print("EBV query - execution time: %s" % (t2 - t1))
            print("********* end DEBUG ***********\n")

            print("E(B-V) data length: %s" % len(ebv))
            min_ebv = np.min(ebv)
            max_ebv = np.max(ebv)
            print("Max E(B-V): %s" % min_ebv)
            print("Min E(B-V): %s" % max_ebv)

            mwe_insert = "INSERT INTO SkyPixel_EBV (N128_SkyPixel_id, EBV) VALUES (%s, %s);"
            mwe_data = [(nside128_ids[i], str(e)) for i,e in enumerate(ebv)]

            batch_insert(mwe_insert, mwe_data)

            with open(pickle_output_dir + 'ebv.pkl', 'wb') as handle:
                pickle.dump(ebv, handle, protocol=pickle.HIGHEST_PROTOCOL)

        ################################################

        ## IS this necessary??
        if not self.options.build_galaxy_skypixel_associations:
            print("Skipping SkyPixel - Galaxy Associations...")
        else:
            ## Build SkyPixel - Galaxy associations
            print("\n\nBuilding SkyPixel - Galaxy Associations...")

            insert_into_SP_G_query = '''
                INSERT INTO SkyPixel_Galaxy (SkyPixel_id, Galaxy_id)
                SELECT sp.id, g.id
                FROM Galaxy g
                JOIN SkyPixel sp on sp.Pixel_Index = g.N%s_Pixel_Index
                WHERE sp.NSIDE = %s AND
                      g.z_dist IS NOT NULL AND
                      g.B IS NOT NULL;
            '''

            queries_to_exe = []
            for i, (nside, pixel_dict) in enumerate(sky_pixels.items()):
                queries_to_exe.append(insert_into_SP_G_query % (nside, nside))

            query_db(queries_to_exe, commit=True)

        ################################################

        if not self.options.build_completeness:
            print("Skipping Completenesses...")
        else:
            print("Building Completenesses...")
            solar_B_abs = 5.48 # mag

            sky_completeness_select = '''
                SELECT 
                    d.id as SkyDistance_id, 
                    sp.id as SkyPixel_id,
                    d.D1, 
                    d.D2, 
                    SUM(POW(10,-0.4*((g.B - (5*log10(g.z_dist*1e+6)-5)) - %s))/1e+10) as L10, 
                    SUM(POW(10,-0.4*((g.B - (5*log10(g.z_dist*1e+6)-5)) - %s))/1e+10)/d.Sch_L10 as Completeness,
                    sp.Pixel_Index 
                FROM teglon.Galaxy g 
                JOIN teglon.SkyPixel_Galaxy sp_g on sp_g.Galaxy_id = g.id 
                JOIN teglon.SkyPixel sp on sp.id = sp_g.SkyPixel_id 
                JOIN teglon.SkyDistance d on d.NSIDE = sp.NSIDE 
                WHERE g.z_dist BETWEEN d.D1 and d.D2 and d.id = %s
                GROUP BY SkyDistance_id, SkyPixel_id, d.D1, d.D2 
                ORDER BY sp.id, d.D1;
            '''

            sky_distance_select = '''
                SELECT 
                    id, 
                    D1, 
                    D2, 
                    dCoV, 
                    Sch_L10, 
                    NSIDE 
                FROM SkyDistance 
                ORDER BY D1
            '''

            sky_completeness_insert = '''
                INSERT INTO SkyCompleteness (SkyPixel_id, SkyDistance_id, L10, Completeness, SmoothedCompleteness)  
                VALUES (%s, %s, %s, %s, %s) 
            '''

            sky_distances = batch_query([sky_distance_select])[0]
            smoothing_radius_Mpc = 30.0  # Mpc
            for sd_index, d in enumerate(sky_distances):

                curr_iter = sd_index + 1
                print("\n**** Creating Completeness Slice %s/%s ****\n" % (curr_iter, len(sky_distances)))

                d_ID = int(d[0])
                d_d1 = float(d[1])
                d_d2 = float(d[2])
                d_NSIDE = int(d[5])

                avg_dist = (d_d1 + d_d2) / 2.
                z = z_at_value(cosmo.luminosity_distance, avg_dist * u.Mpc)
                mpc_radian = cosmo.angular_diameter_distance(z)  # Mpc/radian
                radian_radius = smoothing_radius_Mpc / mpc_radian.value

                completeness_result = query_db([sky_completeness_select % (solar_B_abs, solar_B_abs, d_ID)])[0]

                # retrieve the pixels for a given sky distance
                current_pix = sky_pixels[d_NSIDE]
                current_completenesses = OrderedDict()

                # Initialize all pixels
                for pix_key, pix_val in current_pix.items():
                    current_completenesses[pix_val.index] = [pix_val.id, d_ID, 0.0, 0.0, 0.0, pix_val.index]

                # Update pixels that have completeness records
                for cr in completeness_result:
                    current_pix_index = int(cr[6])
                    cr_L10 = float(cr[4])
                    cr_complete = float(cr[5])

                    current_completenesses[current_pix_index][2] = cr_L10
                    current_completenesses[current_pix_index][3] = cr_complete

                pix_completenesses = []
                for pix_key, pix_val in current_completenesses.items():
                    pix_completenesses.append(pix_val[3])

                smoothed_completeness = hp.sphtfunc.smoothing(np.asarray(pix_completenesses), fwhm=radian_radius, iter=1)

                for sc_index, sc in enumerate(smoothed_completeness):

                    # Smoothing can introduce unphysical completeness values; fix this.
                    recitified_smoothed_completeness = sc
                    if sc < 0.0:
                        recitified_smoothed_completeness = 0.0
                    elif sc > 1.0:
                        recitified_smoothed_completeness = 1.0

                    current_completenesses[sc_index][4] = recitified_smoothed_completeness

                insert_completeness_data = []
                for sp_key, sp_value in current_completenesses.items():
                    insert_completeness_data.append((sp_value[0], sp_value[1], sp_value[2], sp_value[3], sp_value[4]))

                batch_insert(sky_completeness_insert, insert_completeness_data)


        composed_completeness_dict = None
        if not self.options.compose_completeness:
            print("Skipping Compose Completeness...")

            if self.options.is_debug:
                # Check that you can deserialize this...
                if os.path.exists(pickle_output_dir + 'composed_completeness_dict.pkl'):
                    with open(pickle_output_dir + 'composed_completeness_dict.pkl', 'rb') as handle:
                        composed_completeness_dict = pickle.load(handle)
                else:
                    print("composed_completeness_dict.pkl does not exist!")
        else:
            print("Composing Completeness...")
            composed_completeness_dict = OrderedDict()

            n128_ids = []
            for pix_index, pix_ele in sky_pixels[nside128].items():
                n128_ids.append(pix_ele.id)

            select_insert_completeness = '''
                        INSERT INTO CompletenessGrid (MeanDistance, SmoothedCompleteness, N128_SkyPixel_id)
                        SELECT
                            0.5*(sd.D1+sd.D2) as Dist, sc.SmoothedCompleteness, sp1.id as N128_id
                        FROM SkyPixel sp1
                        JOIN SkyPixel sp2 on sp2.id = sp1.Parent_Pixel_id
                        JOIN SkyPixel sp3 on sp3.id = sp2.Parent_Pixel_id
                        JOIN SkyPixel sp4 on sp4.id = sp3.Parent_Pixel_id
                        JOIN SkyPixel sp5 on sp5.id = sp4.Parent_Pixel_id
                        JOIN SkyPixel sp6 on sp6.id = sp5.Parent_Pixel_id
                        JOIN SkyPixel sp7 on sp7.id = sp6.Parent_Pixel_id
                        JOIN SkyCompleteness sc on sc.SkyPixel_id in (sp1.id, sp2.id, sp3.id, sp4.id, sp5.id, sp6.id, sp7.id)
                        JOIN SkyDistance sd on sd.id = sc.SkyDistance_id
                        WHERE
                            sp1.id = %s AND 
                            sp2.NSIDE = 64 AND
                            sp3.NSIDE = 32 AND
                            sp4.NSIDE = 16 AND
                            sp5.NSIDE = 8 AND
                            sp6.NSIDE = 4 AND
                            sp7.NSIDE = 2;
                        '''


            print("Executing INSERT to Completeness Grid...")
            queries_to_exe = []
            for n128_id in n128_ids:
                queries_to_exe.append(select_insert_completeness % n128_id)

            print("\t...sending %s queries" % len(queries_to_exe))
            query_db(queries_to_exe, commit=True)

            select_composed_completeness = '''
                SELECT id, MeanDistance, SmoothedCompleteness, N128_SkyPixel_id 
                FROM CompletenessGrid 
                ORDER BY N128_SkyPixel_id, MeanDistance;
            '''

            composed_completeness_result = query_db([select_composed_completeness])[0]
            for ccr in composed_completeness_result:
                mean_dist = float(ccr[1])
                smoothed_completeness = float(ccr[2])
                n128_pix_id = int(ccr[3])

                if n128_pix_id not in composed_completeness_dict:
                    composed_completeness_dict[n128_pix_id] = [[], []]

                composed_completeness_dict[n128_pix_id][0].append(mean_dist)
                composed_completeness_dict[n128_pix_id][1].append(smoothed_completeness)

            # append a point at infinity
            for key, val in composed_completeness_dict.items():
                composed_completeness_dict[key][0].append(np.inf)
                composed_completeness_dict[key][1].append(0.0)

            with open(pickle_output_dir + 'composed_completeness_dict.pkl', 'wb') as handle:
                pickle.dump(composed_completeness_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

        ################################################

        if not self.options.build_static_grids:
            print("Skipping Static Grids...")
        else:
            print("Building Static Grids...")

            grids_to_build = [
                "SWOPE",
                "THACHER",
                "T80S_T80S-Cam"
            ]

            t1 = time.time()
            self.generate_static_grid_inserts(detector_name_list=grids_to_build, ebv_dict=ebv, skypixel_dict=sky_pixels)
            t2 = time.time()

            print("\n********* start DEBUG ***********")
            print("Static Tile Creation - execution time: %s" % (t2 - t1))
            print("********* end DEBUG ***********\n")

            # select_detect = '''
            #     SELECT Name, ST_AsText(Poly), Deg_width, Deg_height, id, MinDec, MaxDec FROM Detector WHERE Name='%s';
            # '''
            # insert_static_tile = '''INSERT INTO StaticTile (Detector_id, FieldName, RA, _Dec, Coord, Poly, EBV, N128_SkyPixel_id)
            # VALUES (%s, %s, %s, %s, ST_PointFromText(%s, 4326), ST_GEOMFROMTEXT(%s, 4326), %s, %s)'''
            #
            # # Only build out static grids for Swope and Thacher to begin with
            # t1 = time.time()
            # swope_result = query_db([select_detect % 'SWOPE'])[0][0]
            # swope_vertices = Detector.get_detector_vertices_from_teglon_db(swope_result[1])
            # swope_id = int(swope_result[4])
            # swope_min_dec = float(swope_result[5])
            # swope_max_dec = float(swope_result[6])
            # swope = Detector(detector_name="SWOPE", detector_vertex_list_collection=swope_vertices,
            #                  detector_width_deg=swope_result[2], detector_height_deg=swope_result[2],
            #                  detector_id=swope_id)
            # swope_coords = Cartographer.generate_all_sky_coords(swope)
            # t2 = time.time()
            # print("\n********* start DEBUG ***********")
            # print("Swope SkyCoord Creation - execution time: %s" % (t2 - t1))
            # print("********* end DEBUG ***********\n")
            #
            # t1 = time.time()
            # thacher_result = query_db([select_detect % 'THACHER'])[0][0]
            # thacher_vertices = Detector.get_detector_vertices_from_teglon_db(swope_result[1])
            # thacher_id = int(thacher_result[4])
            # thacher_min_dec = float(thacher_result[5])
            # thacher_max_dec = float(thacher_result[6])
            # thacher = Detector(detector_name="THACHER", detector_vertex_list_collection=thacher_vertices,
            #                  detector_width_deg=thacher_result[2], detector_height_deg=thacher_result[2],
            #                  detector_id=thacher_id)
            # thacher_coords = Cartographer.generate_all_sky_coords(thacher)
            # t2 = time.time()
            # print("\n********* start DEBUG ***********")
            # print("Thacher SkyCoord Creation - execution time: %s" % (t2 - t1))
            # print("********* end DEBUG ***********\n")
            #
            # swope_tiles = []
            # t1 = time.time()
            # for i, c in enumerate(swope_coords):
            #     ra_deg, dec_deg = c[0], c[1]
            #     if dec_deg >= swope_min_dec and dec_deg <= swope_max_dec:
            #         t = Tile(c[0], c[1], swope, nside128)
            #         t.field_name = "S%s" % str(i).zfill(6)
            #         n128_index = hp.ang2pix(nside128, 0.5*np.pi - t.dec_rad, t.ra_rad) # theta, phi
            #         t.mwe = ebv[n128_index]
            #         t.id = sky_pixels[nside128][n128_index].id
            #         swope_tiles.append(t)
            #
            # t2 = time.time()
            # print("\n********* start DEBUG ***********")
            # print("Swope Static Tile Creation - execution time: %s" % (t2 - t1))
            # print("********* end DEBUG ***********\n")
            #
            # thacher_tiles = []
            # t1 = time.time()
            # for i, c in enumerate(thacher_coords):
            #     ra_deg, dec_deg = c[0], c[1]
            #     if dec_deg >= thacher_min_dec and dec_deg <= thacher_max_dec:
            #         t = Tile(c[0], c[1], thacher, nside128)
            #         t.field_name = "T%s" % str(i).zfill(6)
            #         n128_index = hp.ang2pix(nside128, 0.5*np.pi - t.dec_rad, t.ra_rad) # theta, phi
            #         t.mwe = ebv[n128_index]
            #         t.id = sky_pixels[nside128][n128_index].id
            #         thacher_tiles.append(t)
            #
            # t2 = time.time()
            # print("\n********* start DEBUG ***********")
            # print("Thacher Static Tile Creation - execution time: %s" % (t2 - t1))
            # print("********* end DEBUG ***********\n")
            #
            # t1 = time.time()
            # static_tile_data = []
            # for t in swope_tiles:
            #     static_tile_data.append((swope_id,
            #         t.field_name,
            #         t.ra_deg,
            #         t.dec_deg,
            #         "POINT(%s %s)" % (t.dec_deg, t.ra_deg - 180.0),  # Dec, RA order due to MySQL convention for lat/lon
            #         t.query_polygon_string,
            #         str(t.mwe),
            #         t.id))
            #
            # for t in thacher_tiles:
            #     static_tile_data.append((thacher_id,
            #         t.field_name,
            #         t.ra_deg,
            #         t.dec_deg,
            #         "POINT(%s %s)" % (t.dec_deg, t.ra_deg - 180.0),  # Dec, RA order due to MySQL convention for lat/lon
            #         t.query_polygon_string,
            #         str(t.mwe),
            #         t.id))
            # t2 = time.time()
            # print("\n********* start DEBUG ***********")
            # print("Building Static Tile Insert Data - execution time: %s" % (t2 - t1))
            # print("********* end DEBUG ***********\n")
            #
            # batch_insert(insert_static_tile, static_tile_data)

if __name__ == "__main__":
    useagestring = """python InitializeTeglon.py [options]

Example non-box command line:
python InitializeTeglon.py  --is_debug True --build_skydistances True --build_skypixels True --build_detectors True
                            --build_bands True --build_MWE True --plot_MWE True 
                            --build_galaxy_skypixel_associations True --build_completeness True 
                            --plot_completeness True --build_static_grids True

All args default to true, to only install/restore certain tables, set non-used args to False
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
    print("Teglon `InitializeTeglon` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")