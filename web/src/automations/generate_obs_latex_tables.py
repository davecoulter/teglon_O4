from collections import OrderedDict
import argparse
import time

from astropy import units as u
import astropy.coordinates as coord
from astropy.table import Table
from astropy.time import Time
from astropy.io import ascii
from dustmaps.sfd import SFDQuery
import healpy as hp

import json

from web.src.utilities.Database_Helpers import *

from web.src.utilities.filesystem_utilties import *

utilities_base_dir = "/app/web/src/utilities"
pickle_output_dir = "%s/pickles/" % utilities_base_dir

def build_latex_table(gw_id, healpix_dir, healpix_file, telescope_ids_to_include):

    formatted_healpix_dir = healpix_dir
    formatted_healpix_dir = formatted_healpix_dir.replace("{GWID}", gw_id)

    healpix_map_select = "SELECT id, RescaledNSIDE, t_0 FROM HealpixMap WHERE GWID = '%s' and Filename = '%s'"
    healpix_map_result = query_db([healpix_map_select % (gw_id, healpix_file)])[0][0]
    healpix_map_id = int(healpix_map_result[0])

    select_tiles = '''
        SELECT
            D.id,
            D.Name,
            OT.RA,
            OT._Dec,
            B.`Name`,
            OT.MJD,
            OT.Mag_Lim
        FROM
            ObservedTile OT
        JOIN Detector D on OT.Detector_id = D.id
        JOIN Band B on OT.Band_id = B.id
        WHERE
            OT.HealpixMap_id = %s AND
            D.id IN (%s)
        ORDER BY
            D.Area
        DESC;
    '''
    tile_result = query_db([select_tiles % (healpix_map_id, telescope_ids_to_include)])[0]

    abbreviation_lookup = {
        "SWOPE":"S",
        "ANDICAM-CCD":"A-CCD",
        "THACHER": "T",
        "NICKEL": "N",
        "ANDICAM-IR": "A-IR",
        # "MOSFIRE": "M",
        "SDSS u": "u",
        "SDSS g": "g",
        "SDSS r": "r",
        "SDSS i": "i",
        "SDSS z": "z",
        "Landolt B": "B",
        "Landolt V": "V",
        "Landolt I": "I",
        "UKIRT J": "J",
        "UKIRT H": "H",
        "UKIRT K": "K"
    }

    output_file = "%s/%s_fields.tex" % (formatted_healpix_dir, gw_id)
    with open(output_file, 'w') as fout:
        print('''\\startlongtable
        \\begin{deluxetable}
        {lccccc}
        \\tabletypesize{\\scriptsize}
        \\tablecaption{1M2H Optical/IR Imaging of the GW190425 Localization Region\\label{tab:observations}}
        \\tablewidth{0pt}
        \\tablehead{
        \\colhead{Source\\tablenotemark{\\footnotesize a}} & 
        \\colhead{$\\alpha$} & 
        \\colhead{$\\delta$} & 
        \\colhead{Date} &
        \\colhead{Filter} &
        \\colhead{Magnitude Limit\\tablenotemark{\\footnotesize b}} \\\\
        &
        (J2000) &
        (J2000) &
        (MJD) &
        &
        (3$\\sigma$)
        }
        \\startdata''', file=fout)

        num_rows = len(tile_result)
        for i, tile in enumerate(tile_result):
            row_num = i + 1
            detector_name = abbreviation_lookup[tile[1]]
            ra_deg = float(tile[2])
            dec_deg = float(tile[3])
            ra_str, dec_str = coord.SkyCoord(ra_deg, dec_deg, unit=(u.deg, u.deg)).to_string(style='hmsdms',
                                                                                             sep=":",
                                                                                             precision=2).split(" ")
            filter_name = abbreviation_lookup[tile[4]]
            mjd = tile[5]
            lim_mag = float(tile[6])

            line_terminator = ""
            if row_num < num_rows:
                line_terminator = "\\\\"

            # S & 01:27:55.560 & -34:41:50.280 & 58715.2170 & $r$ & 20.64 \\
            print("%s & %s & %s & %0.4f & $%s$ & %0.2f %s" % (detector_name, ra_str, dec_str, mjd, filter_name,
                                                                lim_mag, line_terminator), file=fout)
        print('''\\enddata
        
    \\tablecomments{We only include data that 1M2H acquired and reduced. For LCO data referred to in \\autoref{sec:data}, these data will be published in \\lcogtpub. For all other data, see the curated pointings on the Treasure Map \\citep{Wyatt20}.}
    \\tablenotetext{a}{Surveys correspond to Swope (S), Thacher (T), Nickel (N), ANDICAM-CCD (A-CCD), and ANDICAM-IR (A-IR). Imaging as described in \\autoref{sec:data}.}
    \\tablenotetext{b}{In-band 3$\\sigma$ limit for the reported image as described in \\autoref{sec:data} and \\autoref{sec:model_comparison}. All magnitudes are on the AB system \\citep{Oke83}.}
    \\end{deluxetable}''', file=fout)

def build_teglon_perf_table(gw_id, healpix_dir, healpix_file, telescope_ids_to_include):
    formatted_healpix_dir = healpix_dir
    formatted_healpix_dir = formatted_healpix_dir.replace("{GWID}", gw_id)

    healpix_map_select = "SELECT id, RescaledNSIDE, t_0 FROM HealpixMap WHERE GWID = '%s' and Filename = '%s'"
    healpix_map_result = query_db([healpix_map_select % (gw_id, healpix_file)])[0][0]
    healpix_map_id = int(healpix_map_result[0])
    healpix_map_nside = int(healpix_map_result[1])
    area_per_pix = hp.nside2pixarea(healpix_map_nside, degrees=True) # Should be 0.052455852825697924
    print("Area_per_pix: %s" % area_per_pix)

    select_pix_stats_by_detector = '''
        WITH UniquePixels As (
            SELECT
                DISTINCT HP.id            as _2d_Pix_id,
                HPC.id           as _4d_Pix_id,
                HP.Pixel_Index,
                HP.Prob          as _2d_Prob,
                HPC.NetPixelProb as _4d_Prob,
                HP.Mean          as Mean_Dist,
                HP.Stddev        as Dist_Err,
                HPC.PixelCompleteness,
                D.Name as detect_name
            FROM HealpixPixel HP
            JOIN HealpixPixel_Completeness HPC on HPC.HealpixPixel_id = HP.id
            JOIN ObservedTile_HealpixPixel OTHP on HP.id = OTHP.HealpixPixel_id
            JOIN ObservedTile OT on OT.id = OTHP.ObservedTile_id
            JOIN Detector D on OT.Detector_id = D.id
            WHERE HP.HealpixMap_id = %s)
            SELECT
                ut.detect_name,
                SUM(ut._2d_Prob) * 100 as Total_2D_Percent,
                SUM(ut._4d_Prob) * 100 as Total_4D_Percent
            FROM 
                UniquePixels ut
            GROUP BY
                ut.detect_name;
    '''

    select_tile_stats_by_detector = '''
        WITH All_Tiles AS (
            SELECT
                D.`Name` as detector_name,
                D.id as detector_id,
                OT.RA,
                OT._Dec,
                OT.MJD,
                B.`Name` as filter_name,
                OT.Mag_Lim
            FROM ObservedTile OT
            JOIN Detector D on D.id = OT.Detector_id
            JOIN Band B on B.id = OT.Band_id
            WHERE OT.HealpixMap_id = %s)
        SELECT 
            _at.detector_name,
            D.Area as FOV,
            COUNT(_at.MJD) AS TotalTiles
        FROM All_Tiles _at
        JOIN Detector D on D.id = _at.detector_id
        GROUP BY
            _at.detector_name, 
            D.Area
        ORDER BY
            D.Area;
    '''

    select_all_pix_synopsis = '''
        WITH All_Unique AS (
            SELECT
                DISTINCT HP.id            as _2d_Pix_id,
                HPC.id           as _4d_Pix_id,
                HP.Pixel_Index,
                HP.Prob          as _2d_Prob,
                HPC.NetPixelProb as _4d_Prob
            FROM HealpixPixel HP
            JOIN HealpixPixel_Completeness HPC on HPC.HealpixPixel_id = HP.id
            JOIN ObservedTile_HealpixPixel OTHP on HP.id = OTHP.HealpixPixel_id
            WHERE HP.HealpixMap_id = %s)
            SELECT
                COUNT(_2d_Pix_id) * %s as TotalArea,# Area of 
                SUM(_2d_Prob)*100,
                SUM(_4d_Prob)*100
            FROM All_Unique;
    '''

    pix_by_detector_result = query_db([select_pix_stats_by_detector % (healpix_map_id)])[0]
    tiles_by_detector_result = query_db([select_tile_stats_by_detector % (healpix_map_id)])[0]
    all_pix_result = query_db([select_all_pix_synopsis % (healpix_map_id, area_per_pix)])[0]

    by_detector = OrderedDict()

    total_tiles = 0
    for row in tiles_by_detector_result:
        detect_name = row[0]
        fov = float(row[1])
        num_tiles = int(row[2])
        if detect_name not in by_detector:
            by_detector[detect_name] = []

        by_detector[detect_name] += [detect_name, fov, num_tiles]
        total_tiles += num_tiles

    for row in pix_by_detector_result:
        detect_name = row[0]
        detect_2D = float(row[1])
        detect_4D = float(row[2])
        boost = detect_4D/detect_2D
        by_detector[detect_name] += [detect_2D, detect_4D, boost]

    total_2D = 0.0
    total_4D = 0.0
    total_area = 0.0
    for row in all_pix_result:
        total_area = float(row[0])
        total_2D = float(row[1])
        total_4D = float(row[2])

    total_boost = total_4D/total_2D

    composed_rows = []
    for d_name, data_row in by_detector.items():
        composed_rows.append((data_row[0], data_row[1], data_row[2], data_row[3], data_row[4], data_row[5]))

    composed_rows.append(('All Tiles', total_area, total_tiles, total_2D, total_4D, total_boost))

    output_file = "%s/%s_search_statistics.tex" % (formatted_healpix_dir, gw_id)
    with open(output_file, 'w') as fout:
        print('''\\begin{deluxetable*}
            {lccccc}
            \\tabletypesize{\\scriptsize}
            \\tablecaption{GW190425 Search Synopsis\\label{tab:search_synopsis}}
            \\tablewidth{0pt}
            \\tablehead{
            \\colhead{Search Instrument} & 
            \\colhead{FOV} & 
            \\colhead{\# of Images} & 
            \\colhead{Total 2D Probability} & 
            \\colhead{Total Resampled Probability} & 
            \\colhead{Efficiency Increase} & \\\\
            &
            (deg$^{2}$) &
            &
            (\%) &
            (\%) & 
            }
            \\startdata''', file=fout)

        num_rows = len(composed_rows)
        for i, row in enumerate(composed_rows):
            row_num = i + 1
            detect_name = row[0]
            fov = row[1]
            num = row[2]
            _2d = row[3]
            _4d = row[4]
            boost = row[5]

            if row_num < num_rows:
                # Swope & 0.25 & 204 & 0.58 & 1.87 \\
                print("%s & %0.4f & %s & %0.2f & %0.2f & %0.2f \\\\" % (detect_name, fov, num, _2d, _4d, boost), file=fout)
            else:
                print("{\\bf %s} & {\\bf %0.0f}\\tablenotemark{\\footnotesize a} & {\\bf %s} & {\\bf %0.2f} & {\\bf %0.2f} & {\\bf %0.2f} \\\\" %
                      (detect_name, fov, num, _2d, _4d, boost), file=fout)

        print('''\\enddata
        \\tablecomments{}
        \\tablenotetext{a}{Total unique area covered by all observations.}
        \end{deluxetable*}''', file=fout)


if __name__ == "__main__":
    start = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('--gw_id', default="", type=str, help='LIGO superevent name, e.g. `S190425z`')
    parser.add_argument('--healpix_dir', default="./web/events/{GWID}", type=str,
                        help='Directory for where to look for the tile file.')
    parser.add_argument('--healpix_file', default="bayestar.fits.gz", type=str, help='healpix filename')
    parser.add_argument('--telescope_ids_to_include', default="", type=str, help='''Comma-delimited list of telescope 
    ids to include in the Latex table output''')


    args = parser.parse_args()
    build_latex_table(**vars(args))
    build_teglon_perf_table(**vars(args))

    end = time.time()
    duration = (end - start)
    print("\n********* start DEBUG ***********")
    print("Teglon `generate_obs_latex_tables` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")