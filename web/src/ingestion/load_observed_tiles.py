import os
import csv
import time
import pickle
from collections import OrderedDict

import numpy as np
from astropy import units as u
import astropy.coordinates as coord
import healpy as hp
import pprint
from astropy.table import Table
from astropy.io import ascii

from web.src.objects.Detector import *
from web.src.objects.Tile import *
from web.src.objects.Pixel_Element import *
from web.src.utilities.Database_Helpers import *


isDEBUG = False
class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--gw_id', default="", type="str", help='LIGO superevent name, e.g. `S190425z`.')

        parser.add_option('--healpix_file', default="", type="str",
                          help='Healpix filename. Used with `gw_id` to identify unique map.')

        parser.add_option('--healpix_dir', default='./web/events/{GWID}', type="str",
                          help='Directory for where to look for the healpix file.')

        parser.add_option('--tile_dir', default="./web/events/{GWID}/observed_tiles", type="str",
                          help='Directory for where to look for observed tiles to import.')

        parser.add_option('--tile_file', default="", type="str", help='File that contains the tile observations.')

        return (parser)

    def main(self):

        utilities_base_dir = "/app/web/src/utilities"
        pickle_output_dir = "%s/pickles/" % utilities_base_dir
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

        if is_error:
            print("Exiting...")
            return 1

        print("\tLoading NSIDE 128 pixels...")
        nside128 = 128
        N128_dict = None
        print("\tLoading NSIDE 128 pixels...")
        with open(pickle_output_dir + 'N128_dict.pkl', 'rb') as handle:
            N128_dict = pickle.load(handle)
        del handle

        print("\tLoading existing mwe...")
        ebv = None
        if os.path.exists(pickle_output_dir + "ebv.pkl"):
            with open(pickle_output_dir + "ebv.pkl", 'rb') as handle:
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

        duplicate_tile_data = []
        obs_tile_insert_data = {}
        # detectors = {}

        input_table = ascii.read(tile_path, delimiter=' ')
        input_detector_name = input_table.meta["comments"]["detector_name"]
        detector_result = query_db([detector_select_by_name % input_detector_name])[0][0]

        detector_name = detector_result[1]
        detector_poly = Detector.get_detector_vertices_from_teglon_db(detector_result[8])
        detector = Detector(detector_name, detector_poly, detector_id=int(detector_result[0]))

        filter_select = '''
            SELECT `Name`, id FROM Band;
        '''
        filter_result = query_db([filter_select])[0]
        filter_lookup = {br[0]: int(br[1]) for br in filter_result}

        print("Processing `%s` for %s" % (tile_path, detector.name))


        # Sanity checks:
        #   0. Check for duplication: What defines the tile being the same:
        #       (Detector_id, MJD, Filter, N128_SkyPixel_id)
        #       While exp time, mag lim, etc might be different, there's no legitimate way for the same detector to
        #       take > 1 "r-band" image at the exact same.
        #   1. Check for duplication within the db. This preempts duplication within the payload.
        #   2. Check for duplication within the payload.
        #   3. Only upload tile-pixel relations for uploaded tiles (vs all tiles in the db)

        select_existing = '''
            SELECT 
                COUNT(*)
            FROM ObservedTile 
            WHERE
                Detector_id = %s AND
                MJD = %s AND
                Band_Id = %s AND
                N128_SkyPixel_id = %s 
        '''

        for row in input_table:

            source = row['source']
            field_name = row['field_name']
            ra = row['ra']
            dec = row['dec']

            n128_pix_index = hp.ang2pix(128, ra, dec, lonlat=True)
            tile_ebv = ebv[n128_pix_index]
            n128_id = N128_dict[n128_pix_index]

            filter_name = row['filter']
            filter_id = filter_lookup[filter_name]

            mjd = row['mjd']
            exp_time = row['exp_time']
            mag_lim = row['mag_lim']
            position_angle = row['position_angle']

            # 1. db check
            error_msg = "\tDuplicate Tile! (Detector: %s, MJD: %s, Filter_id: %s, N128_id: %s)" % (detector.id, mjd,
                                                                                                 filter_id, n128_id) + \
                        "\nSkipping! ***********\n"
            unique_key = (detector.id, mjd, filter_id, n128_id)
            already_exists = query_db([select_existing % unique_key])[0]
            if len(already_exists) > 0 and already_exists[0][0] >= 1:
                print("\n *********** Duplicated in DB: *********** ")
                print(error_msg)
                duplicate_tile_data.append(unique_key)
                continue

            # 2. payload check
            if unique_key not in obs_tile_insert_data:
                obs_tile_insert_data[unique_key] = None
            else:
                print("\n*********** Duplicate in Payload: *********** ")
                print(error_msg)
                duplicate_tile_data.append(unique_key)
                continue

            x0 = row['x0']
            if np.isnan(x0):
                x0 = None

            x0_err = row['x0_err']
            if np.isnan(x0_err):
                x0_err = None

            a = row['a']
            if np.isnan(a):
                a = None

            a_err = row['a_err']
            if np.isnan(a_err):
                a_err = None

            n = row['n']
            if np.isnan(n):
                n = None

            n_err = row['n_err']
            if np.isnan(n_err):
                n_err = None

            t = Tile(ra, dec, detector, int(healpix_map_nside), position_angle_deg=position_angle)

            obs_tile_insert_data[unique_key] = (
                source,
                detector.id,
                field_name,
                ra,
                dec,
                "POINT(%s %s)" % (dec, ra - 180.0),  # Dec, RA order due to MySQL convention for lat/lon
                t.query_polygon_string,
                float(tile_ebv),
                n128_id,
                filter_id,
                mjd,
                exp_time,
                mag_lim,
                healpix_map_id,
                position_angle,
                x0,
                x0_err,
                a,
                a_err,
                n,
                n_err
            )

        # Convert dictionary to list for db insert
        insert_observed_tile_list = [value for key, value in obs_tile_insert_data.items()]

        insert_observed_tile = '''
            INSERT INTO
                ObservedTile (Source, Detector_id, FieldName, RA, _Dec, Coord, Poly, EBV, N128_SkyPixel_id, Band_id, 
                MJD, Exp_Time, Mag_Lim, HealpixMap_id, PositionAngle, x0, x0_err, a, a_err, n, n_err)
            VALUES (%s, %s, %s, %s, %s, ST_PointFromText(%s, 4326), ST_GEOMFROMTEXT(%s, 4326), %s, %s, %s, %s, %s, %s, 
            %s, %s, %s, %s, %s, %s, %s, %s)
        '''

        print("Inserting %s tiles..." % len(insert_observed_tile_list))
        batch_insert(insert_observed_tile, insert_observed_tile_list)
        print("Done...")

        # Only insert tile-pixel data for things just inserted!!!
        select_inserted = '''
            WITH NewINSERT AS (
                SELECT id FROM ObservedTile ORDER BY id DESC LIMIT %s
            )
            SELECT id FROM NewINSERT ORDER BY id ASC;
        '''
        inserted_ids = query_db([select_inserted % len(insert_observed_tile_list)])[0]
        # generate id list
        id_list = ",".join([str(id[0]) for id in inserted_ids])

        # 3. Only insert tile-pixel data if there were tiles inserted (previous bug would regen tile-pixel relations
        #       for ALL existing tiles...
        if len(insert_observed_tile_list) <= 0:
            print("No tiles inserted! Exiting...")
            return 0;

        # Associate map pixels with observed tiles
        print("Tiles inserted: Building observed tile-healpix map pixel relation...")

        obs_tile_select = '''
            SELECT `Source`, id, Detector_id, FieldName, RA, _Dec, Coord, Poly, EBV, N128_SkyPixel_id, Band_id, MJD, Exp_Time, Mag_Lim, HealpixMap_id, PositionAngle 
            FROM ObservedTile 
            WHERE id IN (%s)
        '''

        # Obtain the tile with id's so they can be used in later INSERTs
        tile_pixel_data = []
        obs_tiles_for_detector = query_db([obs_tile_select % id_list])[0]
        for otfd in obs_tiles_for_detector:

            ot_id = int(otfd[1])
            ra = float(otfd[4])
            dec = float(otfd[5])
            pa = float(otfd[15])

            t = Tile(ra, dec, detector, healpix_map_nside, position_angle_deg=pa)
            t.id = ot_id

            for p in t.enclosed_pixel_indices:
                tile_pixel_data.append((t.id, map_pixel_dict[p][0], healpix_map_id))

        print("Length of tile_pixel_data: %s" % len(tile_pixel_data))

        tile_pixel_upload_csv = "%s/%s_tile_pixel_upload.csv" % (formatted_tile_dir,
                                                                 detector.name.replace("/", "_"))

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
                    (ObservedTile_id, HealpixPixel_id, HealpixMap_id);"""

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

        print("Duplicates removed from input file:")
        pprint.pprint(duplicate_tile_data)


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
