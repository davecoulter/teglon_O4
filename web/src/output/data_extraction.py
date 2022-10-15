import csv
import astropy.coordinates as coord
import numpy as np
from astropy import units as u
import time
import os


from web.src.utilities.Database_Helpers import query_db


# Function used in serializing output
def get_sexigesimal_string(c):
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


def extract_tiles(gw_id, healpix_file, s_band='r', t_band='r', g_band='r', s_exp_time=60.0, t_exp_time=120.0,
                  g_exp_time=120.0, extinct=0.5, prob_type='4D', s_cum_prob=0.9, t_cum_prob=0.9, s_min_ra=-1,
                  s_max_ra=-1, s_min_dec=-1, s_max_dec=-1, t_min_ra=-1, t_max_ra=-1, t_min_dec=-1, t_max_dec=-1,
                  g_min_ra=-1, g_max_ra=-1, g_min_dec=-1, g_max_dec=-1, galaxies_detector='NICKEL', num_gal=10000):


    # HACK - construct file path relative to this script
    path_tokens = os.path.dirname(__file__).split("/")
    healpix_dir = "/".join(path_tokens[0:-2]) + "/events/{gw_id}".format(gw_id=gw_id)

    t1 = time.time()

    # region sanity_checks
    # Valid prob types
    _4D = "4D"
    _2D = "2D"

    # Valid detectors names
    valid_detector_names = {
        "THACHER",
        "SWOPE",
        "NICKEL",
        "SINISTRO"
    }

    band_mapping = {
        "g": "SDSS g",
        "r": "SDSS r",
        "i": "SDSS i"
    }

    is_error = False

    # Parameter sanity checks
    if gw_id == "":
        is_error = True
        print("GWID is required.")

    formatted_healpix_dir = healpix_dir
    if "{GWID}" in formatted_healpix_dir:
        formatted_healpix_dir = formatted_healpix_dir.replace("{GWID}", gw_id)

    if healpix_file == "":
        is_error = True
        print("You must specify which healpix file to process.")

    if s_band not in band_mapping:
        is_error = True
        print("Invalid Swope band selection. Available bands: %s" % band_mapping.keys())

    if t_band not in band_mapping:
        is_error = True
        print("Invalid Thacher band selection. Available bands: %s" % band_mapping.keys())

    if g_band not in band_mapping:
        is_error = True
        print("Invalid Galaxy band selection. Available bands: %s" % band_mapping.keys())

    doSwope = True
    if s_exp_time <= 0.0:
        doSwope = False
        print("SKIPPING Swope!")

    doThacher = True
    if t_exp_time <= 0.0:
        doThacher = False
        print("SKIPPING Thacher!")

    doGalaxies = True
    if g_exp_time <= 0.0:
        doGalaxies = False
        print("SKIPPING Galaxies!")

    if not doSwope and not doThacher and not doGalaxies:
        is_error = True
        print("Skipping all tile extractions!")

    if extinct <= 0.0:
        is_error = True
        print("Extinction must be a valid float > 0.0")

    if not (prob_type == _4D or prob_type == _2D):
        is_error = True
        print("Prob type must either be `4D` or `2D`")

    if s_cum_prob > 0.95 or s_cum_prob < 0.20:
        is_error = True
        print("Swope cumulative prob must be between 0.2 and 0.95")

    if t_cum_prob > 0.95 or t_cum_prob < 0.20:
        is_error = True
        print("Thacher cumulative prob must be between 0.2 and 0.95")

    if galaxies_detector not in valid_detector_names:
        is_error = True
        print("Invalid galaxies detector. Must be one of: %s" % valid_detector_names)

    is_Swope_box_query = False
    if s_min_ra != -1 and s_max_ra != -1 and s_min_dec != -1 and s_max_dec != -1:
        # if any of these are non-default values, then we will validate them all
        is_Swope_box_query = True

        if s_min_ra < 0.0 or s_min_ra >= s_max_ra:
            is_error = True
            print("Swope Min RA must be >= 0.0 and be < Swope Max RA")

        if s_max_ra > 360.0 or s_max_ra <= s_min_ra:
            is_error = True
            print("Swope Max RA must be <= 360.0 and be > Swope Min RA")

        if s_min_dec < -90.0 or s_min_dec >= s_max_dec:
            is_error = True
            print("Swope Min Dec must be >= -90.0 and be < Swope Max Dec")

        if s_max_dec > 90.0 or s_max_dec <= s_min_dec:
            is_error = True
            print("Swope Max Dec must be <= 90.0 and be > Swope Min Dec")

    is_Thacher_box_query = False
    if t_min_ra != -1 and t_max_ra != -1 and t_min_dec != -1 and t_max_dec != -1:
        # if any of these are non-default values, then we will validate them all
        is_Thacher_box_query = True

        if t_min_ra < 0.0 or t_min_ra >= t_max_ra:
            is_error = True
            print("Thacher Min RA must be >= 0.0 and be < Thacher Max RA")

        if t_max_ra > 360.0 or t_max_ra <= t_min_ra:
            is_error = True
            print("Thacher Max RA must be <= 360.0 and be > Thacher Min RA")

        if t_min_dec < -90.0 or t_min_dec >= t_max_dec:
            is_error = True
            print("Thacher Min Dec must be >= -90.0 and be < Thacher Max Dec")

        if t_max_dec > 90.0 or t_max_dec <= t_min_dec:
            is_error = True
            print("Thacher Max Dec must be <= 90.0 and be > Thacher Min Dec")

    is_Galaxy_box_query = False
    if g_min_ra != -1 and g_max_ra != -1 and g_min_dec != -1 and g_max_dec != -1:
        # if any of these are non-default values, then we will validate them all
        is_Galaxy_box_query = True

        if g_min_ra < 0.0 or g_min_ra >= g_max_ra:
            is_error = True
            print("Galaxy Min RA must be >= 0.0 and be < Galaxy Max RA")

        if g_max_ra > 360.0 or g_max_ra <= g_min_ra:
            is_error = True
            print("Galaxy Max RA must be <= 360.0 and be > Galaxy Min RA")

        if g_min_dec < -90.0 or g_min_dec >= g_max_dec:
            is_error = True
            print("Galaxy Min Dec must be >= -90.0 and be < Galaxy Max Dec")

        if g_max_dec > 90.0 or g_max_dec <= g_min_dec:
            is_error = True
            print("Galaxy Max Dec must be <= 90.0 and be > Galaxy Min Dec")
    # endregion

    if is_error:
        print("Exiting...")
        return 1

    # region DB QUERIES
    # Get Map ID
    healpix_map_select = "SELECT id FROM HealpixMap WHERE GWID = '%s' and Filename = '%s'"

    band_select = "SELECT id, Name, F99_Coefficient FROM Band WHERE `Name`='%s'"
    detector_select_by_name = "SELECT id, Name, Deg_width, Deg_height, Deg_radius, Area, MinDec, MaxDec, ST_AsText(Poly) FROM Detector WHERE Name='%s'"

    # 4D NON-BOX QUERY
    # REPLACEMENT PARAMETERS: (band_F99, detector_id, healpix_map_id, self.options.cum_prob, band_F99,
    # self.options.extinct)
    tile_select_4D = '''
                    SELECT 
                        running_tile_prob.id as static_tile_id,  
                        st.FieldName, 
                        st.RA, 
                        st._Dec, 
                        running_tile_prob.net_prob, 
                        running_tile_prob.mean_pixel_dist, 
                        running_tile_prob.EBV*%s as A_lambda, 
                        running_tile_prob.cum_prob 
                    FROM 
                        (SELECT 
                            aggregate_tile.id, 
                            aggregate_tile.Detector_id, 
                            aggregate_tile.net_prob, 
                            aggregate_tile.mean_pixel_dist, 
                            aggregate_tile.EBV, 
                            SUM(aggregate_tile.net_prob) OVER(ORDER BY SUM(aggregate_tile.net_prob) DESC) AS cum_prob 
                         FROM 
                            (SELECT 
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
                            WHERE d.id = %s and hp.HealpixMap_id = %s 
                            GROUP BY 
                                st.id, 
                                st.Detector_id, 
                                sp_ebv.EBV) aggregate_tile 
                        GROUP BY 
                            aggregate_tile.id, 
                            aggregate_tile.Detector_id, 
                            aggregate_tile.net_prob, 
                            aggregate_tile.mean_pixel_dist, 
                            aggregate_tile.EBV 
                        ORDER BY 
                            aggregate_tile.net_prob DESC) as running_tile_prob 
                    JOIN StaticTile st on st.id = running_tile_prob.id 
                    JOIN Detector d on d.id = running_tile_prob.Detector_id 
                    WHERE 
                        running_tile_prob.cum_prob <= %s AND 
                        running_tile_prob.EBV*%s <= %s AND 
                        st._Dec BETWEEN d.MinDec AND d.MaxDec 
                    ORDER BY 
                        running_tile_prob.net_prob DESC 
                    '''

    # 2D NON-BOX QUERY
    # REPLACEMENT PARAMETERS: (band_F99, detector_id, healpix_map_id, self.options.cum_prob, band_F99,
    # self.options.extinct)
    tile_select_2D = '''
                    SELECT 
                        running_tile_prob.id as static_tile_id,  
                        st.FieldName, 
                        st.RA, 
                        st._Dec, 
                        running_tile_prob.net_prob, 
                        running_tile_prob.mean_pixel_dist, 
                        running_tile_prob.EBV*%s as A_lambda, 
                        running_tile_prob.cum_prob 
                    FROM 
                        (SELECT 
                            aggregate_tile.id, 
                            aggregate_tile.Detector_id, 
                            aggregate_tile.net_prob, 
                            aggregate_tile.mean_pixel_dist, 
                            aggregate_tile.EBV, 
                            SUM(aggregate_tile.net_prob) OVER(ORDER BY SUM(aggregate_tile.net_prob) DESC) AS cum_prob 
                         FROM 
                            (SELECT 
                                st.id, 
                                st.Detector_id, 
                                SUM(hp.prob) as net_prob, 
                                AVG(hp.Mean) as mean_pixel_dist, 
                                sp_ebv.EBV 
                            FROM 
                                StaticTile st 
                            JOIN StaticTile_HealpixPixel st_hp on st_hp.StaticTile_id = st.id  
                            JOIN HealpixPixel hp on hp.id = st_hp.HealpixPixel_id 
                            JOIN SkyPixel_EBV sp_ebv on sp_ebv.N128_SkyPixel_id = st.N128_SkyPixel_id 
                            JOIN Detector d on d.id = st.Detector_id 
                            WHERE d.id = %s and hp.HealpixMap_id = %s 
                            GROUP BY 
                                st.id, 
                                st.Detector_id, 
                                sp_ebv.EBV) aggregate_tile 
                        GROUP BY 
                            aggregate_tile.id, 
                            aggregate_tile.Detector_id, 
                            aggregate_tile.net_prob, 
                            aggregate_tile.mean_pixel_dist, 
                            aggregate_tile.EBV 
                        ORDER BY 
                            aggregate_tile.net_prob DESC) as running_tile_prob 
                    JOIN StaticTile st on st.id = running_tile_prob.id 
                    JOIN Detector d on d.id = running_tile_prob.Detector_id 
                    WHERE 
                        running_tile_prob.cum_prob <= %s AND 
                        running_tile_prob.EBV*%s <= %s AND 
                        st._Dec BETWEEN d.MinDec AND d.MaxDec 
                    ORDER BY 
                        running_tile_prob.net_prob DESC 
                '''

    # 4D BOX QUERY
    # REPLACEMENT PARAMETERS: (band_F99, detector_id, healpix_map_id, self.options.cum_prob, band_F99,
    # self.options.extinct, s/t_min_ra, s/t_max_ra, s/t_min_dec, s/t_max_dec)
    box_tile_select_4D = '''
                    SELECT 
                        running_tile_prob.id as static_tile_id,  
                        st.FieldName, 
                        st.RA, 
                        st._Dec, 
                        running_tile_prob.net_prob, 
                        running_tile_prob.mean_pixel_dist, 
                        running_tile_prob.EBV*%s as A_lambda, 
                        running_tile_prob.cum_prob 
                    FROM 
                        (SELECT 
                            aggregate_tile.id, 
                            aggregate_tile.Detector_id, 
                            aggregate_tile.net_prob, 
                            aggregate_tile.mean_pixel_dist, 
                            aggregate_tile.EBV, 
                            SUM(aggregate_tile.net_prob) OVER(ORDER BY SUM(aggregate_tile.net_prob) DESC) AS cum_prob 
                         FROM 
                            (SELECT 
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
                            WHERE d.id = %s and hp.HealpixMap_id = %s 
                            GROUP BY 
                                st.id, 
                                st.Detector_id, 
                                sp_ebv.EBV) aggregate_tile 
                        GROUP BY 
                            aggregate_tile.id, 
                            aggregate_tile.Detector_id, 
                            aggregate_tile.net_prob, 
                            aggregate_tile.mean_pixel_dist, 
                            aggregate_tile.EBV 
                        ORDER BY 
                            aggregate_tile.net_prob DESC) as running_tile_prob 
                    JOIN StaticTile st on st.id = running_tile_prob.id 
                    JOIN Detector d on d.id = running_tile_prob.Detector_id 
                    WHERE 
                        running_tile_prob.cum_prob <= %s AND 
                        running_tile_prob.EBV*%s <= %s AND  
                        st.RA BETWEEN %s AND %s AND 
                        st._Dec BETWEEN %s AND %s 
                    ORDER BY 
                        running_tile_prob.net_prob DESC 
                    '''

    # 2D BOX QUERY
    # REPLACEMENT PARAMETERS: (band_F99, detector_id, healpix_map_id, self.options.cum_prob, band_F99,
    # self.options.extinct, s/t_min_ra, s/t_max_ra, s/t_min_dec, s/t_max_dec)
    box_tile_select_2D = '''
                    SELECT 
                        running_tile_prob.id as static_tile_id,  
                        st.FieldName, 
                        st.RA, 
                        st._Dec, 
                        running_tile_prob.net_prob, 
                        running_tile_prob.mean_pixel_dist, 
                        running_tile_prob.EBV*%s as A_lambda, 
                        running_tile_prob.cum_prob 
                    FROM 
                        (SELECT 
                            aggregate_tile.id, 
                            aggregate_tile.Detector_id, 
                            aggregate_tile.net_prob, 
                            aggregate_tile.mean_pixel_dist, 
                            aggregate_tile.EBV, 
                            SUM(aggregate_tile.net_prob) OVER(ORDER BY SUM(aggregate_tile.net_prob) DESC) AS cum_prob 
                         FROM 
                            (SELECT 
                                st.id, 
                                st.Detector_id, 
                                SUM(hp.Prob) as net_prob, 
                                AVG(hp.Mean) as mean_pixel_dist, 
                                sp_ebv.EBV 
                            FROM 
                                StaticTile st 
                            JOIN StaticTile_HealpixPixel st_hp on st_hp.StaticTile_id = st.id  
                            JOIN HealpixPixel hp on hp.id = st_hp.HealpixPixel_id 
                            JOIN SkyPixel_EBV sp_ebv on sp_ebv.N128_SkyPixel_id = st.N128_SkyPixel_id 
                            JOIN Detector d on d.id = st.Detector_id 
                            WHERE d.id = %s and hp.HealpixMap_id = %s 
                            GROUP BY 
                                st.id, 
                                st.Detector_id, 
                                sp_ebv.EBV) aggregate_tile 
                        GROUP BY 
                            aggregate_tile.id, 
                            aggregate_tile.Detector_id, 
                            aggregate_tile.net_prob, 
                            aggregate_tile.mean_pixel_dist, 
                            aggregate_tile.EBV 
                        ORDER BY 
                            aggregate_tile.net_prob DESC) as running_tile_prob 
                    JOIN StaticTile st on st.id = running_tile_prob.id 
                    JOIN Detector d on d.id = running_tile_prob.Detector_id 
                    WHERE 
                        running_tile_prob.cum_prob <= %s AND 
                        running_tile_prob.EBV*%s <= %s AND  
                        st.RA BETWEEN %s AND %s AND 
                        st._Dec BETWEEN %s AND %s 
                    ORDER BY 
                        running_tile_prob.net_prob DESC 
                    '''

    # GALAXIES SELECT
    # REPLACEMENT PARAMETERS: (band_F99, healpix_map_id, band_F99, self.options.extinct, g_detector_min_dec,
    # g_detector_max_dec)
    galaxies_select = '''
                SELECT 
                    g.id, 
                    g.Name_GWGC, 
                    g.Name_HyperLEDA, 
                    g.Name_2MASS, 
                    g.RA, 
                    g._Dec, 
                    g.z_dist, 
                    g.B, 
                    g.K, 
                    hp_g_w.GalaxyProb, 
                    sp_ebv.EBV*%s AS A_lambda 
                FROM 
                    Galaxy g 
                JOIN HealpixPixel_Galaxy_Weight hp_g_w on hp_g_w.Galaxy_id = g.id 
                JOIN HealpixPixel hp on hp.id = hp_g_w.HealpixPixel_id 
                JOIN SkyPixel_EBV sp_ebv on sp_ebv.N128_SkyPixel_id = hp.N128_SkyPixel_id 
                WHERE 
                    hp.HealpixMap_id = %s AND 
                    sp_ebv.EBV*%s <= %s AND 
                    g._Dec BETWEEN %s AND %s; 
            '''

    # GALAXIES SELECT
    # REPLACEMENT PARAMETERS: (band_F99, healpix_map_id, band_F99, self.options.extinct, g_min_ra, g_max_ra,
    # g_min_dec, g_max_dec)
    box_galaxies_select = '''
                SELECT 
                    g.id, 
                    g.Name_GWGC, 
                    g.Name_HyperLEDA, 
                    g.Name_2MASS, 
                    g.RA, 
                    g._Dec, 
                    g.z_dist, 
                    g.B, 
                    g.K, 
                    hp_g_w.GalaxyProb, 
                    sp_ebv.EBV*%s AS A_lambda 
                FROM 
                    Galaxy g 
                JOIN HealpixPixel_Galaxy_Weight hp_g_w on hp_g_w.Galaxy_id = g.id 
                JOIN HealpixPixel hp on hp.id = hp_g_w.HealpixPixel_id 
                JOIN SkyPixel_EBV sp_ebv on sp_ebv.N128_SkyPixel_id = hp.N128_SkyPixel_id 
                WHERE 
                    hp.HealpixMap_id = %s AND 
                    sp_ebv.EBV*%s <= %s AND 
                    g.RA BETWEEN %s AND %s AND 
                    g._Dec BETWEEN %s AND %s; 
            '''
    # endregion

    SWOPE = "SWOPE"
    THACHER = "THACHER"
    GALAXY = galaxies_detector

    hpx_path = "%s/%s" % (formatted_healpix_dir, healpix_file)
    s_band_name = band_mapping[s_band]
    t_band_name = band_mapping[t_band]
    g_band_name = band_mapping[g_band]

    healpix_map_id = int(query_db([healpix_map_select % (gw_id, healpix_file)])[0][0][0])

    # E.g. Non box:
    #   <directory>/SWOPE_4D_0.9_LALInference_0.fits.gz.csv
    # Box:
    #   <directory>/SWOPE_2D_0.8_bayestar.fits.gz_box.csv
    box_file_formatter = "%s/%s_%s_%s_%s_box.csv"
    non_box_file_formatter = "%s/%s_%s_%s_%s.csv"
    box_galaxies_formatter = "%s/%s_%s_box.csv"
    non_box_galaxies_formatter = "%s/%s_%s.csv"

    s_formatted_output_path = ""
    t_formatted_output_path = ""
    g_formatted_output_path = ""
    r_formatted_output_path = "%s/%s_%s.reg" % (formatted_healpix_dir, "REGION",
                                                healpix_file.replace(",", "_"))

    if doSwope:
        print("Extracting Swope...")
        s_detector_result = query_db([detector_select_by_name % SWOPE])[0][0]
        s_detector_id = int(s_detector_result[0])

        s_band_result = query_db([band_select % s_band_name])[0][0]
        s_band_F99 = float(s_band_result[2])

        s_select_to_execute = ""
        if is_Swope_box_query:
            s_formatted_output_path = box_file_formatter % (formatted_healpix_dir, SWOPE, prob_type,
                                                            s_cum_prob, healpix_file.replace(",", "_"))
            if prob_type == _4D:
                s_select_to_execute = box_tile_select_4D % (s_band_F99, s_detector_id, healpix_map_id,
                                                            s_cum_prob, s_band_F99, extinct,
                                                            s_min_ra, s_max_ra, s_min_dec, s_max_dec)
            else:
                s_select_to_execute = box_tile_select_2D % (s_band_F99, s_detector_id, healpix_map_id,
                                                            s_cum_prob, s_band_F99, extinct,
                                                            s_min_ra, s_max_ra, s_min_dec, s_max_dec)
        else:
            s_formatted_output_path = non_box_file_formatter % (
                formatted_healpix_dir, SWOPE, prob_type, s_cum_prob, healpix_file.replace(",", "_"))
            if prob_type == _4D:
                s_select_to_execute = tile_select_4D % (s_band_F99, s_detector_id, healpix_map_id,
                                                        s_cum_prob, s_band_F99, extinct)
            else:
                s_select_to_execute = tile_select_2D % (s_band_F99, s_detector_id, healpix_map_id,
                                                        s_cum_prob, s_band_F99, extinct)

        s_tile_result = query_db([s_select_to_execute])[0]



        with open(s_formatted_output_path, 'w') as csvfile:
            csvwriter = csv.writer(csvfile)

            cols = []
            cols.append('# FieldName')
            cols.append('FieldRA')
            cols.append('FieldDec')
            cols.append('Telscope')
            cols.append('Filter')
            cols.append('ExpTime')
            cols.append('Priority')
            cols.append('Status')
            cols.append('A_lambda')
            csvwriter.writerow(cols)

            for i, row in enumerate(s_tile_result):
                c = coord.SkyCoord(row[2], row[3], unit=(u.deg, u.deg))
                coord_str = get_sexigesimal_string(c)

                cols = []

                cols.append(row[1])
                cols.append(coord_str[0])
                cols.append(coord_str[1])
                cols.append(SWOPE)
                cols.append(s_band)
                cols.append(str(s_exp_time))
                cols.append(row[4])
                cols.append('False')
                cols.append("%0.3f" % float(row[6]))
                csvwriter.writerow(cols)

            print("Done writing out %s" % SWOPE)

    if doThacher:
        print("Extracting Thacher...")
        t_detector_result = query_db([detector_select_by_name % THACHER])[0][0]
        t_detector_id = int(t_detector_result[0])

        t_band_result = query_db([band_select % t_band_name])[0][0]
        t_band_F99 = float(t_band_result[2])

        t_select_to_execute = ""
        if is_Thacher_box_query:
            t_formatted_output_path = box_file_formatter % (formatted_healpix_dir, THACHER, prob_type,
                                                            t_cum_prob, healpix_file.replace(",", "_"))
            if prob_type == _4D:
                t_select_to_execute = box_tile_select_4D % (t_band_F99, t_detector_id, healpix_map_id,
                                                            t_cum_prob, t_band_F99, extinct,
                                                            t_min_ra, t_max_ra, t_min_dec, t_max_dec)
            else:
                t_select_to_execute = box_tile_select_2D % (t_band_F99, t_detector_id, healpix_map_id,
                                                            t_cum_prob, t_band_F99, extinct,
                                                            t_min_ra, t_max_ra, t_min_dec, t_max_dec)
        else:
            t_formatted_output_path = non_box_file_formatter % (
                formatted_healpix_dir, THACHER, prob_type, t_cum_prob, healpix_file.replace(",", "_"))
            if prob_type == _4D:
                t_select_to_execute = tile_select_4D % (t_band_F99, t_detector_id, healpix_map_id,
                                                        t_cum_prob, t_band_F99, extinct)
            else:
                t_select_to_execute = tile_select_2D % (t_band_F99, t_detector_id, healpix_map_id,
                                                        t_cum_prob, t_band_F99, extinct)

        t_tile_result = query_db([t_select_to_execute])[0]

        with open(t_formatted_output_path, 'w') as csvfile:
            csvwriter = csv.writer(csvfile)

            cols = []
            cols.append('# FieldName')
            cols.append('FieldRA')
            cols.append('FieldDec')
            cols.append('Telscope')
            cols.append('Filter')
            cols.append('ExpTime')
            cols.append('Priority')
            cols.append('Status')
            cols.append('A_lambda')
            csvwriter.writerow(cols)

            for i, row in enumerate(t_tile_result):
                c = coord.SkyCoord(row[2], row[3], unit=(u.deg, u.deg))
                coord_str = get_sexigesimal_string(c)

                cols = []

                cols.append(row[1])
                cols.append(coord_str[0])
                cols.append(coord_str[1])
                cols.append(THACHER)
                cols.append(t_band)
                cols.append(str(t_exp_time))
                cols.append(row[4])
                cols.append('False')
                cols.append("%0.3f" % float(row[6]))
                csvwriter.writerow(cols)

            print("Done writing out %s" % THACHER)

    if doGalaxies:
        print("Extracting Galaxies for %s..." % GALAXY)
        g_detector_result = query_db([detector_select_by_name % GALAXY])[0][0]
        g_detector_min_dec = float(g_detector_result[6])
        g_detector_max_dec = float(g_detector_result[7])

        g_band_result = query_db([band_select % g_band_name])[0][0]
        g_band_F99 = float(g_band_result[2])

        g_select_to_execute = ""
        if is_Galaxy_box_query:
            g_formatted_output_path = box_galaxies_formatter % (formatted_healpix_dir, "GALAXIES",
                                                                healpix_file.replace(",", "_"))
            # (f99, healpix_map_id, f99, extinct, ra1, ra2, dec1, dec2)
            g_select_to_execute = box_galaxies_select % (g_band_F99, healpix_map_id, g_band_F99, extinct,
                                                         g_min_ra, g_max_ra, g_min_dec, g_max_dec)
        else:
            g_formatted_output_path = non_box_galaxies_formatter % (formatted_healpix_dir, "GALAXIES",
                                                                    healpix_file.replace(",", "_"))
            # (f99, healpix_mpa_id, f99, extinction, min_dec, max_dec)
            g_select_to_execute = galaxies_select % (g_band_F99, healpix_map_id, g_band_F99, extinct,
                                                     g_detector_min_dec, g_detector_max_dec)

        galaxies_result = query_db([g_select_to_execute])[0]

        # Sort/limit galaxies_result
        sorted_galaxies_result = sorted(galaxies_result, key=lambda x: float(x[9]), reverse=True)
        galaxies_to_write = []
        if num_gal > 0:
            galaxies_to_write = sorted_galaxies_result[0:num_gal]

        with open(g_formatted_output_path, 'w') as csvfile:

            csvwriter = csv.writer(csvfile)

            cols = []
            cols.append('# FieldName')
            cols.append('FieldRA')
            cols.append('FieldDec')
            cols.append('Telscope')
            cols.append('Filter')
            cols.append('ExpTime')
            cols.append('Priority')
            cols.append('Status')
            cols.append('A_lambda')
            csvwriter.writerow(cols)

            for i, row in enumerate(galaxies_to_write):

                c = coord.SkyCoord(row[4], row[5], unit=(u.deg, u.deg))
                coord_str = get_sexigesimal_string(c)

                cols = []

                Name_GWGC = row[1]
                Name_HyperLEDA = row[2]
                Name_2MASS = row[3]

                field_name = ""
                if Name_GWGC is not None:
                    field_name = Name_GWGC
                elif Name_HyperLEDA is not None:
                    field_name = "LEDA" + Name_HyperLEDA
                elif Name_2MASS is not None:
                    field_name = Name_2MASS

                if field_name == "":
                    raise ("No field name!")

                cols.append(field_name)
                cols.append(coord_str[0])
                cols.append(coord_str[1])
                cols.append(GALAXY)
                cols.append(g_band)
                cols.append(str(g_exp_time))
                cols.append(row[9])
                cols.append('False')
                cols.append("%0.3f" % float(row[10]))
                csvwriter.writerow(cols)

            print("Done w/ Galaxies")

        with open(r_formatted_output_path, 'w') as csvfile:

            csvfile.write("# Region file format: DS9 version 4.0 global\n\n")
            csvfile.write("global color=lightgreen\n")
            csvfile.write("ICRS\n")

            for i, row in enumerate(galaxies_to_write):

                Name_GWGC = row[1]
                Name_HyperLEDA = row[2]
                Name_2MASS = row[3]

                field_name = ""
                if Name_GWGC is not None:
                    field_name = Name_GWGC
                elif Name_HyperLEDA is not None:
                    field_name = "LEDA" + Name_HyperLEDA
                elif Name_2MASS is not None:
                    field_name = Name_2MASS

                if field_name == "":
                    raise ("No field name!")

                csvfile.write('circle(%s,%s,120") # width=4 text="%s"\n' % (row[4], row[5], field_name))

            print("Done w/ Region File")

    t2 = time.time()
    print("\n********* start DEBUG ***********")
    print("Extract Tiles execution time: %s" % (t2 - t1))
    print("********* end DEBUG ***********\n")

    return s_formatted_output_path
