import csv
import os.path
import sys

import astropy.coordinates as coord
import numpy as np
from astropy import units as u
from web.src.utilities.Database_Helpers import *
from web.src.utilities.filesystem_utilties import *
from astropy.table import Table
from astropy.io import fits

from collections import OrderedDict


class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--gw_id', default="", type="str",
                          help='LIGO superevent name, e.g. `S190425z` ')

        parser.add_option('--healpix_dir', default='./web/events/{GWID}', type="str",
                          help='Directory for where to look for the healpix file.')

        parser.add_option('--healpix_file', default="bayestar.fits.gz", type="str", help='healpix filename.')

        parser.add_option('--tele', default="a", type="str",
                          help='''Telescope abbreviation for the telescope to extract files for. Default: `a`. 
                                  Available: (s:Swope, t:Thacher, n:Nickel, t80: T80S_T80S-Cam, a:All). If `All` 
                                  selected,combination plot is generated (vs single plot).''')

        parser.add_option('--band', default="r", type="str",
                          help='Only used if `tele` specified. Default: `r`. Available: (g, r, i, z, I)')

        parser.add_option('--extinct', default="0.5", type="float",
                          help='''Extinction in mags in the specified band(s) to be less than. Default: 0.5 mag. 
                          Must be > 0.0''')

        parser.add_option('--prob_type', default="4D", type="str",
                          help='''Probability type to consider. Default `4D`. Available: (4D, 2D)''')

        parser.add_option('--cum_prob', default="0.9", type="float",
                          help='Cumulative prob to cover. Default 0.9. Must be > 0.2 and < 0.95')

        parser.add_option('--num_tiles', default="1000", type="int",
                          help='''Return the top `num_tiles` tiles per telescope. Must be > 1. Default 1000.''')

        parser.add_option('--min_ra', default="-1.0", type="float",
                          help='''Optional, decimal format. If specified, `min_ra`, `max_ra`, `min_dec`, 
                      and `max_dec` are all required. `min_ra` must be >= 0.0 and < `max_ra`. Within the defined box, 
                      tiles will be extracted within a given telescope's declination limits.''')

        parser.add_option('--max_ra', default="-1.0", type="float",
                          help='''Optional, decimal format. If specified, `min_ra`, `max_ra`, `min_dec`, 
                          and `max_dec` are all required. `max_ra` must be <= 360.0 and > `min_ra`. Within the defined 
                          box, tiles will be extracted within a given telescope's declination limits.''')

        parser.add_option('--max_dec', default="-1.0", type="float",
                          help='''Optional, decimal format. Must be If specified, `min_ra`, `max_ra`, `min_dec`, 
                          and `max_dec` are all required. `max_dec` must be <= 90.0 and > `min_dec`. Within the defined 
                          box, tiles will be extracted within a given  telescope's declination limits.''')

        parser.add_option('--min_dec', default="-1.0", type="float",
                          help='''Optional, decimal format. Must be If specified, `min_ra`, `max_ra`, `min_dec`, 
                          and `max_dec` are all required. `min_dec` must be >= -90.0 and < `max_dec`. Within the defined 
                          box, tiles will be extracted within a given telescope's declination limits.''')

        return (parser)

    def compose_tile_query(self, detector_id, band_F99, extinct, healpix_map_id, tile_where, tile_aggregate,
                           tile_composed, num_tiles):

        w = tile_where.format(f99_coefficient=band_F99, extinct=extinct)
        a = tile_aggregate.format(detector_id=detector_id, healpix_map_id=healpix_map_id)
        tile_composed = tile_composed.format(aggregate_tile_cte=a, f99_coefficient=band_F99,
                                             where_filter=w,
                                             num_tiles=num_tiles)
        return tile_composed

    def compose_galaxy_query(self, band_F99, extinct, galaxy_select, galaxy_where, num_tiles):

        w = galaxy_where.format(f99_coefficient=band_F99, extinct=extinct)
        galaxy_composed = galaxy_select.format(where_filter=w, f99_coefficient=band_F99, num_tiles=num_tiles)
        return galaxy_composed

    def GetSexigesimalString(self, c):
        ra = c.ra.hms
        dec = c.dec.dms

        ra_string = "%02d:%02d:%05.2f" % (ra[0], ra[1], ra[2])
        if dec[0] >= 0:
            dec_string = "+%02d:%02d:%05.2f" % (dec[0], np.abs(dec[1]), np.abs(dec[2]))
        else:
            dec_string = "%03d:%02d:%05.2f" % (dec[0], np.abs(dec[1]), np.abs(dec[2]))

        # Python has a -0.0 object. If the deg is this (because object lies < 60 min south), the string formatter will drop the negative sign
        if c.dec < 0.0 and dec[0] == 0.0:
            dec_string = "-00:%02d:%05.2f" % (np.abs(dec[1]), np.abs(dec[2]))
        return (ra_string, dec_string)

    def main(self):

        print(self.options)

        t1 = time.time()

        # Valid prob types
        _4D = "4D"
        _2D = "2D"

        detector_mapping = {
            "s": "SWOPE",
            "t": "THACHER",
            "n": "NICKEL",
            "t80": "T80S_T80S-Cam",
            "a": "All"
        }

        band_mapping = {
            "g": "SDSS g",
            "r": "SDSS r",
            "i": "SDSS i",
            "z": "SDSS z",
            "I": "Landolt I"
        }

        default_settings = {
            "SWOPE": "r",
            "THACHER": "r",
            "NICKEL": "r",
            "T80S_T80S-Cam": "I",
        }

        is_error = False
        # region Parameter sanity checks
        formatted_healpix_dir = self.options.healpix_dir

        if self.options.gw_id == "":
            is_error = True
            print("GWID is required.")
        else:
            if "{GWID}" in formatted_healpix_dir:
                formatted_healpix_dir = formatted_healpix_dir.replace("{GWID}", self.options.gw_id)

        healpix_file_path = "{base_dir}/{file_name}".format(base_dir=formatted_healpix_dir,
                                                            file_name=self.options.healpix_file)
        if self.options.healpix_file == "":
            is_error = True
            print("You must specify which healpix file to process.")
        else:
            if not os.path.exists(healpix_file_path):
                is_error = True
                print("Healpix file does not exist!")

        if self.options.tele not in detector_mapping:
            is_error = True
            print("Invalid telescope selection. Available telescopes: %s" % detector_mapping.keys())

        if self.options.band not in band_mapping:
            is_error = True
            print("Invalid band selection. Available bands: %s" % band_mapping.keys())

        if self.options.extinct <= 0.0:
            is_error = True
            print("Extinction must be a valid float > 0.0")

        if self.options.num_tiles < 1:
            is_error = True
            print("Number of tiles  must > 1")

        if not (self.options.prob_type == _4D or self.options.prob_type == _2D):
            is_error = True
            print("Prob type must either be `4D` or `2D`")

        if self.options.cum_prob > 0.95 or self.options.cum_prob < 0.20:
            is_error = True
            print("Cumulative prob must be between 0.2 and 0.95")

        is_box_query = False
        if self.options.min_ra != -1 and \
                self.options.max_ra != -1 and \
                self.options.min_dec != -1 and \
                self.options.max_dec != -1:
            # if any of these are non-default values, then we will validate them all
            is_box_query = True

            if self.options.min_ra < 0.0 or self.options.min_ra >= self.options.max_ra:
                is_error = True
                print("Min RA must be >= 0.0 and be < Max RA")

            if self.options.max_ra > 360.0 or self.options.max_ra <= self.options.min_ra:
                is_error = True
                print("Max RA must be <= 360.0 and be > Min RA")

            if self.options.min_dec < -90.0 or self.options.min_dec >= self.options.max_dec:
                is_error = True
                print("Min Dec must be >= -90.0 and be < Max Dec")

            if self.options.max_dec > 90.0 or self.options.max_dec <= self.options.min_dec:
                is_error = True
                print("Max Dec must be <= 90.0 and be > Min Dec")

        if is_error:
            print("Exiting...")
            return 1
        # endregion

        extract_all = True
        if self.options.tele != 'a':
            extract_all = False

        # region SQL Statements
        healpix_map_select = "SELECT id FROM HealpixMap WHERE GWID = '%s' and Filename = '%s'"
        band_select = "SELECT id, Name, F99_Coefficient FROM Band WHERE `Name`='%s'"
        detector_select_by_name = "SELECT id, Name, Deg_width, Deg_height, Deg_radius, Area, MinDec, MaxDec, ST_AsText(Poly) FROM Detector WHERE Name='%s'"

        composed_tile_query = '''
            WITH
            {aggregate_tile_cte},
            running_tile_prob AS (
                SELECT 
                    aggregate_tile.id, 
                    aggregate_tile.Detector_id, 
                    aggregate_tile.net_prob, 
                    aggregate_tile.mean_pixel_dist, 
                    aggregate_tile.EBV, 
                    SUM(aggregate_tile.net_prob) OVER(ORDER BY SUM(aggregate_tile.net_prob) DESC) AS cum_prob 
                 FROM 
                    aggregate_tile 
                GROUP BY 
                    aggregate_tile.id, 
                    aggregate_tile.Detector_id, 
                    aggregate_tile.net_prob, 
                    aggregate_tile.mean_pixel_dist, 
                    aggregate_tile.EBV 
                ORDER BY 
                    aggregate_tile.net_prob DESC
            )
            SELECT 
                running_tile_prob.id as static_tile_id,  
                st.FieldName, 
                st.RA, 
                st._Dec, 
                running_tile_prob.net_prob,
                running_tile_prob.EBV,
                running_tile_prob.EBV*{f99_coefficient} as A_lambda, 
                running_tile_prob.mean_pixel_dist, 
                running_tile_prob.cum_prob 
            FROM running_tile_prob 
            JOIN StaticTile st on st.id = running_tile_prob.id 
            JOIN Detector d on d.id = running_tile_prob.Detector_id 
            WHERE 
                {where_filter}
            ORDER BY 
                running_tile_prob.net_prob DESC 
            LIMIT 1,{num_tiles};
        '''

        aggregate_tile_4D_cte = '''
            aggregate_tile AS (
                SELECT 
                    st.id, 
                    st.Detector_id, 
                    SUM(hpc.NetPixelProb) as net_prob, 
                    AVG(hp.Mean) as mean_pixel_dist, 
                    sp_ebv.EBV 
                FROM 
                    StaticTile st 
                JOIN StaticTile_HealpixPixel st_hp on st_hp.StaticTile_id = st.id 
                JOIN HealpixPixel_Completeness hpc on hpc.HealpixPixel_id = st_hp.HealpixPixel_id 
                JOIN HealpixPixel hp on hp.id = hpc.HealpixPixel_id 
                JOIN SkyPixel_EBV sp_ebv on sp_ebv.N128_SkyPixel_id = st.N128_SkyPixel_id 
                JOIN Detector d on d.id = st.Detector_id 
                WHERE d.id = {detector_id} and hp.HealpixMap_id = {healpix_map_id} 
                GROUP BY 
                    st.id, 
                    st.Detector_id, 
                    sp_ebv.EBV
            )
        '''

        aggregate_tile_2D_cte = '''
            aggregate_tile AS (
                SELECT 
                    st.id, 
                    st.Detector_id, 
                    SUM(hp.prob) as net_prob, 
                    AVG(hp.Mean) as mean_pixel_dist, 
                    sp_ebv.EBV 
                FROM 
                    StaticTile st 
                JOIN StaticTile_HealpixPixel st_hp on st_hp.StaticTile_id = st.id 
                JOIN HealpixPixel_Completeness hpc on hpc.HealpixPixel_id = st_hp.HealpixPixel_id 
                JOIN HealpixPixel hp on hp.id = hpc.HealpixPixel_id 
                JOIN SkyPixel_EBV sp_ebv on sp_ebv.N128_SkyPixel_id = st.N128_SkyPixel_id 
                JOIN Detector d on d.id = st.Detector_id 
                WHERE d.id = {detector_id} and hp.HealpixMap_id = {healpix_map_id} 
                GROUP BY 
                    st.id, 
                    st.Detector_id, 
                    sp_ebv.EBV
            )
        '''

        tile_where_filter = '''
            running_tile_prob.net_prob > 0.0 AND
            running_tile_prob.cum_prob <= {cum_prob} AND
            running_tile_prob.EBV*{f99_coefficient} <= {extinct} 
            {box_addendum}
        '''

        tile_box_addendum = '''
            AND st.RA BETWEEN {min_ra} AND {max_ra} 
            AND st._Dec BETWEEN {min_dec} AND {max_dec} 
        '''

        galaxies_select = '''
            SELECT 
                g.id, 
                CASE
                    WHEN NOT (ISNULL(g.Name_GWGC) OR g.Name_GWGC = '' OR g.Name_GWGC = 'null') THEN g.Name_GWGC
                    WHEN NOT (ISNULL(g.PGC) OR g.PGC = '' OR g.PGC = 'null') THEN g.PGC
                    WHEN NOT (ISNULL(g.Name_HyperLEDA) OR g.Name_HyperLEDA = '' OR g.Name_HyperLEDA = 'null') THEN g.Name_HyperLEDA
                    WHEN NOT (ISNULL(g.Name_2MASS) OR g.Name_2MASS = '' OR g.Name_2MASS = 'null') THEN g.Name_2MASS
                    WHEN NOT (ISNULL(g.NAME_SDSS_DR12) OR g.NAME_SDSS_DR12 = '' OR g.NAME_SDSS_DR12 = 'null') THEN g.NAME_SDSS_DR12
                  END as FieldName,
                g.RA, 
                g._Dec,
                hp_g_w.GalaxyProb,
                sp_ebv.EBV, 
                sp_ebv.EBV*{f99_coefficient} AS A_lambda,
                g.z_dist, 
                g.B, 
                g.K
            FROM 
                Galaxy g 
            JOIN HealpixPixel_Galaxy_Weight hp_g_w on hp_g_w.Galaxy_id = g.id 
            JOIN HealpixPixel hp on hp.id = hp_g_w.HealpixPixel_id 
            JOIN SkyPixel_EBV sp_ebv on sp_ebv.N128_SkyPixel_id = hp.N128_SkyPixel_id 
            WHERE 
                {where_filter}
            ORDER BY hp_g_w.GalaxyProb DESC
            LIMIT 1,{num_tiles};
        '''

        galaxies_where_filter = '''
            hp_g_w.GalaxyProb > 0.0 AND
            hp.HealpixMap_id = {healpix_map_id} AND 
            sp_ebv.EBV*{f99_coefficient} <= {extinct} AND 
            g._Dec BETWEEN {min_dec} AND {max_dec} 
            {box_addendum}
        '''

        galaxies_box_addendum = '''
            AND g.RA BETWEEN {min_ra} AND {max_ra}
        '''
        # endregion

        healpix_map_id = int(query_db([healpix_map_select % (self.options.gw_id, self.options.healpix_file)])[0][0][0])

        # region Query Builder
        tile_where = tile_where_filter.format(box_addendum="", cum_prob=self.options.cum_prob,
                                              f99_coefficient="{f99_coefficient}", extinct="{extinct}")
        galaxy_where = galaxies_where_filter.format(box_addendum="", f99_coefficient="{f99_coefficient}",
                                                    healpix_map_id=healpix_map_id, extinct="{extinct}",
                                                    min_dec="{min_dec}", max_dec="{max_dec}")
        if is_box_query:
            tile_box = tile_box_addendum.format(min_ra=self.options.min_ra, max_ra=self.options.max_ra,
                                                min_dec=self.options.min_dec, max_dec=self.options.max_dec)
            tile_where = tile_where_filter.format(box_addendum=tile_box, cum_prob=self.options.cum_prob,
                                                  f99_coefficient="{f99_coefficient}", extinct="{extinct}")

            galaxy_box = galaxies_box_addendum.format(min_ra=self.options.min_ra, max_ra=self.options.max_ra)
            galaxy_where = galaxies_where_filter.format(box_addendum=galaxy_box, f99_coefficient="{f99_coefficient}",
                                                        healpix_map_id=healpix_map_id, extinct="{extinct}",
                                                        min_dec=self.options.min_dec, max_dec=self.options.max_dec)

        tile_aggregate = aggregate_tile_4D_cte.format(detector_id="{detector_id}", healpix_map_id=healpix_map_id)
        if self.options.prob_type == "2D":
            tile_aggregate = aggregate_tile_2D_cte.format(detector_id="{detector_id}", healpix_map_id=healpix_map_id)

        queries_to_execute = {}
        if extract_all:
            # for key, detector_name in detector_mapping.items():
            for detector_name, band_shortname in default_settings.items():

                detector_result = query_db([detector_select_by_name % detector_name])[0][0]
                detector_id = int(detector_result[0])
                detector_min_dec = float(detector_result[6])
                detector_max_dec = float(detector_result[7])

                band_name = band_mapping[band_shortname]
                band_result = query_db([band_select % band_name])[0][0]
                band_F99 = float(band_result[2])

                if detector_name != 'NICKEL':
                    queries_to_execute[(detector_name, band_name, band_F99)] = self.compose_tile_query(detector_id,
                                                                                                       band_F99,
                                                                                                       self.options.extinct,
                                                                                                       healpix_map_id,
                                                                                                       tile_where,
                                                                                                       tile_aggregate,
                                                                                                       composed_tile_query,
                                                                                                       self.options.num_tiles)
                else:

                    w = galaxy_where
                    if not is_box_query:
                        w = galaxy_where.format(f99_coefficient="{f99_coefficient}", extinct="{extinct}",
                                                min_dec=detector_min_dec, max_dec=detector_max_dec)

                    queries_to_execute[(detector_name, band_name, band_F99)] = self.compose_galaxy_query(band_F99,
                                                                                                         self.options.extinct,
                                                                                                         galaxies_select,
                                                                                                         w,
                                                                                                         self.options.num_tiles)
        else:
            detector_name = detector_mapping[self.options.tele]
            detector_result = query_db([detector_select_by_name % detector_name])[0][0]
            detector_id = int(detector_result[0])
            detector_min_dec = float(detector_result[6])
            detector_max_dec = float(detector_result[7])

            band_name = band_mapping[self.options.band]
            band_result = query_db([band_select % band_name])[0][0]
            band_F99 = float(band_result[2])

            if detector_name != 'NICKEL':

                queries_to_execute[(detector_name, band_name, band_F99)] = self.compose_tile_query(detector_id,
                                                                                                   band_F99,
                                                                                                   self.options.extinct,
                                                                                                   healpix_map_id,
                                                                                                   tile_where,
                                                                                                   tile_aggregate,
                                                                                                   composed_tile_query,
                                                                                                   self.options.num_tiles)
            else:
                w = galaxy_where
                if not is_box_query:
                    w = galaxy_where.format(f99_coefficient="{f99_coefficient}", extinct="{extinct}",
                                            min_dec=detector_min_dec, max_dec=detector_max_dec)

                queries_to_execute[(detector_name, band_name, band_F99)] = self.compose_galaxy_query(band_F99,
                                                                                                     self.options.extinct,
                                                                                                     galaxies_select, w,
                                                                                                     self.options.num_tiles)
        # endregion

        # Execute queries and build ECSV output
        output_file_formatter = "{working_dir}/{gw_id}_{detector_name}_{prob_type}_{cum_prob}_{hpx_file}{box_query}.txt"
        box_query_addendum = ""
        if is_box_query:
            box_query_addendum = "_box"
        for (detector_name, band_name, band_F99), query in queries_to_execute.items():
            query_results = query_db([query])[0]
            file_output = output_file_formatter.format(working_dir=formatted_healpix_dir, gw_id=self.options.gw_id,
                                                       detector_name=detector_name, prob_type=self.options.prob_type,
                                                       cum_prob=self.options.cum_prob,
                                                       hpx_file=self.options.healpix_file, box_query=box_query_addendum)

            cols = ['Field_Name', 'RA', 'Dec', 'Prob', 'EBV', 'A_lambda', 'Lum_Dist', 'B_mag', 'K_mag']
            dtype = ['U64', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8']

            output_table = Table(dtype=dtype, names=cols)

            for qr in query_results:
                row_to_add = [qr[1], qr[2], qr[3], qr[4], qr[5], qr[6], qr[7], None, None]
                if detector_name == "NICKEL":
                    row_to_add[-2] = qr[8]
                    row_to_add[-1] = qr[9]
                output_table.add_row(row_to_add)

            hdul = fits.open(healpix_file_path)
            hdr = hdul[1].header
            del hdr['HISTORY']
            hdr_dict = dict(hdr)
            output_table.meta['keywords'] = OrderedDict(hdr_dict)

            comments = OrderedDict()
            comments["GW_ID"] = self.options.gw_id
            comments["HPX_DIR"] = formatted_healpix_dir
            comments["HPX_FILE"] = self.options.healpix_file
            comments["TELESCOPE"] = detector_name
            comments["A_LAMBDA_BAND"] = band_name
            comments["F99_BAND_COEF"] = band_F99
            comments["EXTINCT_THRESH"] = self.options.extinct
            comments["PROB_TYPE"] = self.options.prob_type
            comments["CUM_PROB"] = self.options.cum_prob
            comments["NUM_TILES"] = self.options.num_tiles
            comments["BOX_QUERY"] = is_box_query
            comments["MIN_RA"] = self.options.min_ra
            comments["MAX_RA"] = self.options.max_ra
            comments["MIN_DEC"] = self.options.min_dec
            comments["MAX_DEC"] = self.options.max_dec
            s = ""
            for arg in sys.argv:
                s = s + arg + " "
            cmd = os.path.basename(sys.executable) + " " + s.strip()
            comments["CMD"] = cmd

            output_table.meta['comments'] = comments
            output_table.write(file_output, overwrite=True, format='ascii.ecsv')
            chmod_outputfile(file_output)

        t2 = time.time()
        print("\n********* start DEBUG ***********")
        print("Extract Tiles execution time: %s" % (t2 - t1))
        print("********* end DEBUG ***********\n")


if __name__ == "__main__":
    useagestring = """python ExtractTiles.py [options]

Example non-box command line:
python ExtractTiles.py   --gw_id <gw id> --healpix_file <file name> 
                    --s_band r --s_exp_time 60 --s_cum_prob 0.90 
                    --t_band r --t_exp_time 120 --t_cum_prob 0.75 
                    --extinct 0.5 --prob_type 4D

Example box command line:
python ExtractTiles.py   --gw_id <gw id> --healpix_file <file name> 
                    --s_band r --s_exp_time 60 --s_cum_prob 0.90 
                    --s_min_ra 180.0 --s_max_ra 270.0 --s_min_dec -65.0 --s_max_dec 25.0 
                    --t_band r --t_exp_time 120 --t_cum_prob 0.75 
                    --t_min_ra 90.0 --t_max_ra 115.0 --t_min_dec -20.0 --t_max_dec 90.0 
                    --extinct 0.5 --prob_type 4D

NOTE: Passing s_/t_/g_exp_time <= 0 will SKIP that extraction!
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
    print("Teglon `ExtractTiles` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")
