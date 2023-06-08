from collections import OrderedDict
import argparse
import time

from astropy import units as u
import astropy.coordinates as coord
from astropy.table import Table
from astropy.time import Time
from astropy.io import ascii
from dustmaps.sfd import SFDQuery

import json

from web.src.objects.Detector import *
from web.src.objects.Tile import *
# from web.src.utilities.Database_Helpers import *

from web.src.utilities.filesystem_utilties import *

table_keywords = ['input_file_name', 'detector_name']
target_table_names = ('source', 'field_name', 'ra', 'dec', 'filter', 'mjd', 'exp_time', 'mag_lim', 'position_angle',
                      'x0', 'x0_err', 'a', 'a_err', 'n', 'n_err')
target_table_row = [['X' * 40], ['X' * 40], [0.], [0.], ['X' * 40], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.]]

utilities_base_dir = "/app/web/src/utilities"
pickle_output_dir = "%s/pickles/" % utilities_base_dir

def blank_target_table():
    tab = Table(target_table_row, names=target_table_names)
    return tab[:0].copy()


def is_number(val):
    try:
        val = float(val)
        return True
    except ValueError:
        return False


def parse_coord(ra, dec):
    if not (is_number(ra) and is_number(dec)) and (':' not in str(ra) and ':' not in str(dec)):
        return None

    if ':' in str(ra) and ':' in str(dec):
        # Input RA/DEC are sexagesimal
        unit = (u.hourangle, u.deg)
    else:
        unit = (u.deg, u.deg)

    try:
        c = coord.SkyCoord(ra, dec, frame='icrs', unit=unit)
        return(c)
    except ValueError:
        return(None)


def convert_obs_file_to_escv(gw_id, healpix_dir, tile_file, tele):

    obs_tiles_dir = healpix_dir + "/observed_tiles/{tile_file}"
    formatted_obs_tiles_dir = obs_tiles_dir.replace("{GWID}", gw_id)
    file_path = formatted_obs_tiles_dir.replace("{tile_file}", tile_file)

    # print("\tLoading existing E(B-V)...")
    # ebv = None
    # if os.path.exists(pickle_output_dir + "ebv.pkl"):
    #     with open(pickle_output_dir + "ebv.pkl", 'rb') as handle:
    #         ebv = pickle.load(handle)
    # else:
    #     select_ebv = '''
    #     SELECT
    #         ebv.EBV
    #     FROM SkyPixel_EBV ebv
    #     JOIN SkyPixel sp on sp.id = ebv.N128_SkyPixel_id
    #     ORDER BY sp.Pixel_Index;
    #     '''
    #     ebv = np.asarray(query_db([select_ebv])[0])

    # Check detector name
    # detector_select_by_name = "SELECT id, Name FROM Detector WHERE Name='{detect_name}';"
    # detector_result = query_db([detector_select_by_name.format(detect_name=tele)])[0]
    #
    # detect_name = tele
    # if len(detector_result) != 1:
    #     print("\nUnknown/ambiguous detector: %s" % tele)
    #     print("Exiting...")
    #     return 1
    # else:
    #     detect_name = detector_result[0][1]

    detect_name = tele

    band_mapping = {
        "u": "SDSS u",
        "g": "SDSS g",
        "r": "SDSS r",
        "i": "SDSS i",
        "z": "SDSS z",
        "B": "Landolt B",
        "V": "Landolt V",
        "R": "Landolt R",
        "I": "Landolt I",
        "J": "UKIRT J",
        "H": "UKIRT H",
        "K": "UKIRT K",
        "Clear": "Clear",
        "open": "Clear"
    }

    table_to_convert = ascii.read(file_path, delimiter=' ')

    # Sanitize columns
    for key in table_to_convert.keys():
        newkey = key.lower().replace(' ', '_')
        table_to_convert.rename_column(key, newkey)

    position_angle_in_keys = False
    for key in table_to_convert.keys():
        if key in ['fieldra', 'r.a.', 'right_ascension', 'RA']:
            table_to_convert.rename_column(key, 'ra')
        if key in ['fielddec', 'declination', 'dec.', 'Dec']:
            table_to_convert.rename_column(key, 'dec')
        if key in ['fieldname', 'object', 'name', 'field', 'Object']:
            table_to_convert.rename_column(key, 'field_name')
        if key in ['Exp_Time']:
            table_to_convert.rename_column(key, 'exp_time')
        if key in ['p.a.', 'angle', 'P.A.']:
            position_angle_in_keys = True
            table_to_convert.rename_column(key, 'position_angle')
        if key in ['band_pass', 'band', 'Filter']:
            table_to_convert.rename_column(key, 'filter')
        if key in ['obs_date', 'date', 'MJD']:
            table_to_convert.rename_column(key, 'mjd')
        if key in ['lim_mag', 'mag']:
            table_to_convert.rename_column(key, 'mag_lim')

    output_table = blank_target_table()

    for row in table_to_convert:
        new_row = {}

        new_row['field_name'] = str(row['field_name'])

        c = parse_coord(row['ra'], row['dec'])
        new_row['ra'] = c.ra.degree
        new_row['dec'] = c.dec.degree

        # # get EBV for this sight line
        # sfd = SFDQuery()
        # ebv = sfd(c)
        # new_row['EBV'] = ebv

        # Try parsing date:
        tm = None
        try:
            tm = Time(row['mjd'], format='mjd')
        except:
            print("\nError! Date field not in MJD format!")
            print("Exiting...")
            return 1
        new_row['mjd'] = tm.to_value(format="mjd")

        new_row['filter'] = band_mapping[row['filter']]
        new_row['exp_time'] = row['exp_time']

        if position_angle_in_keys:
            new_row['position_angle'] = row['position_angle']
        else:
            new_row['position_angle'] = None

        new_row['x0'] = None
        new_row['x0_err'] = None
        new_row['a'] = None
        new_row['a_err'] = None
        new_row['n'] = None
        new_row['n_err'] = None

        if detect_name in ["SWOPE", "THACHER", "NICKEL", "ANDICAM-CCD", "ANDICAM-IR"]:
            new_row['source'] = row['file']
            new_row['mag_lim'] = row['x0']
            new_row['x0'] = row['x0']
            new_row['x0_err'] = row['x0_err']
            new_row['a'] = row['a']
            new_row['a_err'] = row['a_err']
            new_row['n'] = row['n']
            new_row['n_err'] = row['n_err']
        else:
            new_row['source'] = row['source']
            new_row['mag_lim'] = row['mag_lim']

        output_table.add_row(new_row)

    comments = OrderedDict()
    comments["input_file_name"] = tile_file
    comments["detector_name"] = detect_name
    output_table.meta['comments'] = comments

    file_output = file_path + ".ecsv"
    output_table.write(file_output, overwrite=True, format='ascii.ecsv')
    chmod_outputfile(file_output)

if __name__ == "__main__":
    start = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('--gw_id', default="", type=str, help='LIGO superevent name, e.g. `S190425z`')
    parser.add_argument('--healpix_dir', default="./web/events/{GWID}", type=str,
                        help='Directory for where to look for the tile file.')
    parser.add_argument('--tile_file', default="", type=str, help='Input tile filename.')
    parser.add_argument('--tele', default="", type=str,
                        help='''Telescope name for the telescope to extract files for. Must be a name known 
                        in the Teglon db.''')

    args = parser.parse_args()
    convert_obs_file_to_escv(**vars(args))

    end = time.time()
    duration = (end - start)
    print("\n********* start DEBUG ***********")
    print("Teglon `extract_tiles` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")


