# region imports
import matplotlib as mpl
print(mpl.__version__)


import sys
sys.path.append('../')

import os
import sys
import pickle
# from scipy.interpolate import interp1d, interp2d, spline
from scipy.interpolate import interp1d, interp2d
from web.src.objects.Detector import *
from web.src.objects.Tile import *
from web.src.objects.Pixel_Element import *
from web.src.utilities.Database_Helpers import *

from collections import OrderedDict
from astropy.table import Table
from astropy.time import Time, TimeDelta
import pytz
import csv
from astroplan import moon_illumination, moon_phase_angle
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.patheffects as path_effects
from configparser import RawConfigParser
import json

# endregion

from astropy.cosmology import LambdaCDM
from astropy.coordinates import Distance

from web.src.utilities.filesystem_utilties import *


class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        return (parser)

    def main(self):

        plot_models = False
        plot_data = True

        # print("Hello World!")
        fig = plt.figure(figsize=(15, 15), dpi=600)
        ax = fig.add_subplot(111)

        if plot_models:
            ### Villar Models ###

            # _0817_blue = "./web/models/kne/21/villar_0.0275_0.2493_1.0000.dat" # from Kilpatrick
            _0817_blue = "./web/models/kne/8/villar_0.0182_0.2650_1.0000.dat"  # from Villar
            blue_model_table = Table.read(_0817_blue, format='ascii.ecsv')
            blue_mask = blue_model_table['sdss_r'] != np.nan
            blue_model_time = np.asarray(blue_model_table['time'][blue_mask])

            blue_r = np.asarray(blue_model_table['sdss_r'][blue_mask])
            blue_f_r = interp1d(blue_model_time, blue_r, fill_value="extrapolate")

            # _0817_red = "./web/models/kne/7/villar_0.0338_0.1553_10.0000.dat" # from Kilpatrick
            _0817_red = "./web/models/kne/32/villar_0.0120_0.1397_10.0000.dat" # From Villar
            red_model_table = Table.read(_0817_red, format='ascii.ecsv')
            red_mask = red_model_table['sdss_r'] != np.nan
            red_model_time = np.asarray(red_model_table['time'][red_mask])

            red_r = np.asarray(red_model_table['sdss_r'][red_mask])
            red_f_r = interp1d(red_model_time, red_r, fill_value="extrapolate")

            ### Metzger Models ###

            # Get the Kilpatrick model #s because we want to compare efficiency plots
            _0817_blue_metzger = "./web/models/kne/metzger_models/8/Metzger_0.024890_0.256897_0.450000.dat"  # from Kilpatrick
            blue_model_table_metzger = Table.read(_0817_blue_metzger, format='ascii.ecsv')
            blue_mask_metzger = blue_model_table_metzger['sdss_r'] != np.nan
            blue_model_time_metzger = np.asarray(blue_model_table_metzger['time'][blue_mask_metzger])

            blue_r_metzger = np.asarray(blue_model_table_metzger['sdss_r'][blue_mask_metzger])
            blue_f_r_metzger = interp1d(blue_model_time_metzger, blue_r_metzger, fill_value="extrapolate")

            _0817_red_metzger = "./web/models/kne/metzger_models/4/Metzger_0.038208_0.143448_0.100000.dat"  # from Kilpatrick
            red_model_table_metzger = Table.read(_0817_red_metzger, format='ascii.ecsv')
            red_mask_metzger = red_model_table_metzger['sdss_r'] != np.nan
            red_model_time_metzger = np.asarray(red_model_table_metzger['time'][red_mask_metzger])

            red_r_metzger = np.asarray(red_model_table_metzger['sdss_r'][red_mask_metzger])
            red_f_r_metzger = interp1d(red_model_time_metzger, red_r_metzger, fill_value="extrapolate")

            # _0817_blue = "./web/models/kne/8/villar_0.0182_0.2650_1.0000.dat"  # from Villar
            # _0817_red = "./web/models/kne/32/villar_0.0120_0.1397_10.0000.dat" # From Villar
            ax.plot(red_model_time, red_f_r(red_model_time), color='red', linewidth=2.0,
                    label='Villar r-band Red: Mej = 0.0120, vej = 0.1397, kappa = 10.0')
            ax.plot(blue_model_time, blue_f_r(blue_model_time), color='blue', linewidth=2.0,
                    label='Villar r-band Blue: Mej = 0.0182, vej = 0.2650, kappa = 1.0')
            # Model - Add luminosities together:
            m_tot = -2.5 * np.log10(10 ** (-0.4 * red_f_r(red_model_time)) + 10 ** (-0.4 * blue_f_r(red_model_time)))
            ax.plot(red_model_time, m_tot, color='black', linestyle="-", linewidth=2.0,
                    label='Villar r-band Red + Blue')

            # _0817_blue_metzger = "./web/models/kne/metzger_models/8/Metzger_0.024890_0.256897_0.450000.dat"  # from Kilpatrick
            # _0817_red_metzger = "./web/models/kne/metzger_models/4/Metzger_0.038208_0.143448_0.100000.dat"  # from Kilpatrick
            ax.plot(red_model_time_metzger, red_f_r_metzger(red_model_time_metzger), color='red', linewidth=2.0,
                    linestyle="--",
                    label='Metzger r-band Red: Mej = 0.038208, vej = 0.143448, Ye = 0.10')
            ax.plot(blue_model_time_metzger, blue_f_r_metzger(blue_model_time_metzger), color='blue', linewidth=2.0,
                    linestyle="--",
                    label='Metzger r-band Blue: Mej = 0.024890, vej = 0.256897, Ye = 0.45')
            # Model - Add luminosities together:
            m_tot_metzger = -2.5 * np.log10(10 ** (-0.4 * red_f_r_metzger(red_model_time_metzger)) + 10 ** (
                        -0.4 * blue_f_r_metzger(red_model_time_metzger)))
            ax.plot(red_model_time_metzger, m_tot_metzger, color='black', linestyle="--", linewidth=2.0,
                    label='Metzger r-band Red + Blue')

        if plot_data:

            ### Get Kilonova.space data
            gw170817_filters = ["g", "r", "i"]
            gw170817_json_file = open("./web/models/kne/GW170817.json")
            gw170817_json_data = json.load(gw170817_json_file)
            json_KN_merger_mjd = float(gw170817_json_data["GW170817"]["timeofmerger"][0]["value"])
            sss17a_json_phot = gw170817_json_data["GW170817"]["photometry"]

            gw170817_r_band_phot = []
            gw170817_r_band_phot_err = []
            gw170817_r_band_phot_delta_mjd = []
            for p in sss17a_json_phot:
                if ('band' in p) and ('u_time' in p) and ('upperlimit' not in p) and ('system' in p) and ('e_magnitude' in p) and (p['u_time'] == "MJD") and ('model' not in p) and (p['system'] == 'AB'):
                    if (p['band'] == 'r'):
                        gw170817_r_band_phot_delta_mjd.append(float(p['time']) - json_KN_merger_mjd)
                        gw170817_r_band_phot.append(float(p['magnitude']))
                        gw170817_r_band_phot_err.append(float(p['e_magnitude']))

            sorted_indices = np.argsort(gw170817_r_band_phot_delta_mjd)
            gw170817_mjd_sorted = np.asarray(gw170817_r_band_phot_delta_mjd)[sorted_indices]
            gw170817_r_sorted = np.asarray(gw170817_r_band_phot)[sorted_indices]
            gw170817_err_sorted = np.asarray(gw170817_r_band_phot_err)[sorted_indices]

            # convert to app map
            mu_0817 = 5*np.log10(43.2 * 1e6)-5
            gw170817_r_abs_mag = gw170817_r_sorted - mu_0817

            ### Grab all r-band limits we have in the db
            r_band_select = '''
                    SELECT
                        ot.MJD, 
                        ot.Mag_Lim,
                        d.Name
                    FROM ObservedTile ot
                    JOIN Band b on b.id = ot.Band_id
                    JOIN Detector d on d.id = ot.Detector_id
                    WHERE Band_id = 19
                    ORDER BY MJD;
                    '''

            GW190425_merger_mjd = 58598.346134259256
            mu_0425 = 5*np.log10(156 * 1e6)-5
            print(mu_0425)
            # r_band_limits_abs_mag = []
            # r_band_limits_mjd = []
            r_band_detector_abs_mag = {}
            r_band_detector_mjd = {}
            rb_results = query_db([r_band_select])[0]
            for rbr in rb_results:

                if rbr[2] not in r_band_detector_abs_mag:
                    r_band_detector_abs_mag[rbr[2]] = []
                    r_band_detector_mjd[rbr[2]] = []
                else:
                    r_band_detector_mjd[rbr[2]].append(float(rbr[0]) - GW190425_merger_mjd)
                    r_band_detector_abs_mag[rbr[2]].append(float(rbr[1]) - mu_0425)





            # # Name, A_lambda
            # SDSS u, 0.5126587253999999
            # SDSS g, 0.3994601958
            # SDSS r, 0.276344701
            # SDSS i, 0.2053537428
            # SDSS z, 0.1527454518
            # Landolt B, 0.43852336359999994
            # Landolt V, 0.3316136412
            # Landolt R, 0.2623158234
            # Landolt I, 0.18201259299999997
            # UKIRT J, 0.0857454674
            # UKIRT H, 0.054301431399999996
            # UKIRT K, 0.0365234572
            # Clear, 0.110054126
            # ATLAS cyan, 0.35151045967199995
            # ATLAS orange, 0.25633734464759994
            # PS1 w, 0.28311726260000003




            # GW170817 r-band Data
            ax.errorbar(gw170817_mjd_sorted, gw170817_r_abs_mag, yerr=gw170817_err_sorted, color='black', fmt='*',
                        linestyle="None", label='AT 2017gfo r-band')


            # plot limits
            for key, value in r_band_detector_abs_mag.items():
                detect_name = key
                detect_abs_mag = value
                detect_mdj = r_band_detector_mjd[key]

                ax.plot(detect_mdj, detect_abs_mag, linestyle="None", marker="v",
                        label='%s r-band limits' % key, alpha=0.5)
                # ax.plot(r_band_limits_mjd, r_band_limits_abs_mag, color='green', linestyle="None", marker="v",
                #         label='All r-band limits', alpha=0.5)

        ax.set_ylabel("Absolute Mag. (AB)", fontsize=32, labelpad=9.0)
        ax.set_xlabel("Days From Merger", fontsize=32)

        # ax.set_ylim([16, 28])
        ax.set_ylim([-18, -10])
        ax.set_xlim([0, 14.0])
        ax.invert_yaxis()
        ax.legend()

        ax.tick_params(axis='both', which='major', labelsize=24, length=12.0, width=2)
        ax.grid(linestyle=":")



        # Get LC for GW190425 red model
        # Get LC for Swope red component
        # Get LC for Swope blue component
        # Get real LC for Swope
        # Get GRB170817A LC
        # Get all data upper limits


        output_file = "S190425_LCs.png"
        output_path = "./web/src/analysis/%s" % output_file

        fig.savefig(output_path, bbox_inches='tight', format="png")
        chmod_outputfile(output_path)
        plt.close('all')
        print("... Done.")





if __name__ == "__main__":
    useagestring = """python plot_light_curves.py

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
    print("Teglon `plot_light_curves` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")
