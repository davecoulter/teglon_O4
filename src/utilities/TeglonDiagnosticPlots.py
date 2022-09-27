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


from src.utilities.HEALPix_Helpers import *
from src.utilities.Database_Helpers import *

from src.objects.Detector import *
from src.objects.Tile import *
from src.objects.SQL_Polygon import *
from src.objects.Pixel_Element import *
from src.objects.Completeness_Objects import *


# from mpl_toolkits.basemap import Basemap

import multiprocessing as mp

class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--plot_pix_relations', action="store_true", default=False,
                          help='Plot parent-child pixel relations (Default = False)')

        parser.add_option('--plot_MWE', action="store_true", default=False,
                          help='Plot EBV (Default = False)')

        parser.add_option('--plot_completeness', action="store_true", default=False,
                          help='Expensive! Plot all completeness slices (Default = False)')

        return (parser)

    def main(self):

        # Sanity
        if not self.options.plot_pix_relations and not self.options.plot_MWE and not self.options.plot_completeness:
            print("Nothing to do!")
            return 0

        def plot_completeness(d1, d2, pixels, pixel_outline, file_name):
            lon0 = 180.0
            fig = plt.figure(figsize=(8, 8), dpi=180)
            ax = fig.add_subplot(111)
            m = Basemap(projection='moll', lon_0=lon0)

            norm = colors.Normalize(0.0, 1.0)

            values = []
            for p in pixels:
                pprob = p.prob
                if p.prob <= 0.0:
                    pprob = 0.0
                elif p.prob > 1.0:
                    pprob = 1.0

                values.append(pprob)

            for i, p in enumerate(pixels):
                p.plot(m, ax, facecolor=plt.cm.viridis(values[i]), edgecolor='None', linewidth=0.5, alpha=0.8,
                       zorder=9900)

            for p in pixel_outline:
                p.plot(m, ax, facecolor='None', edgecolor='k', linewidth=0.5)

            meridians = np.arange(0., 360., 60.)
            m.drawparallels(np.arange(-90., 91., 30.), fontsize=14, labels=[True, True, False, False], dashes=[2, 2],
                            linewidth=0.5, xoffset=2500000)
            m.drawmeridians(meridians, labels=[False, False, False, False], dashes=[2, 2], linewidth=0.5)

            for mer in meridians[1:]:
                plt.annotate("%0.0f" % mer, xy=m(mer, 0), xycoords='data', fontsize=14, zorder=9999)

            sm = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.viridis)
            sm.set_array([])  # can be an empty list

            tks = np.linspace(0.0, 1.0, 11)
            tks_strings = []

            for t in tks:
                tks_strings.append('%0.0f' % (t * 100))

            top_left = 1.05
            delta_y = 0.04
            ax.annotate('[%0.0f - %0.0f] Mpc' % (d1, d2), xy=(0.5, top_left), xycoords='axes fraction', fontsize=16,
                        ha='center')

            cb = fig.colorbar(sm, ax=ax, ticks=tks, orientation='horizontal', fraction=0.08951, pad=0.02, alpha=0.80)
            cb.ax.set_xticklabels(tks_strings, fontsize=14)
            cb.set_label("% Complete per Pixel", fontsize=14, labelpad=10.0)
            cb.outline.set_linewidth(1.0)

            for axis in ['top', 'bottom', 'left', 'right']:
                ax.spines[axis].set_linewidth(2.0)

            ax.invert_xaxis()

            fig.savefig(file_name, bbox_inches='tight')  # ,dpi=840
            plt.close('all')

            print("... %s Done." % file_name)

        diag_output_dir = "./diagnostic_plots/"
        if not os.path.exists(diag_output_dir):
            os.makedirs(diag_output_dir)

        nside2 = 2
        nside128 = 128
        sky_pixels = None
        print("\tLoading sky pixels...")
        with open('./pickles/sky_pixels.pkl', 'rb') as handle:
            sky_pixels = pickle.load(handle)

        if self.options.plot_pix_relations:
            print("\n\nTest pixel parent-child relationships...")
            # Pick a random nside32 pixel and plot all of its parents
            test_index = 110435
            print("Plotting n128 index[%s]" % test_index)
            n128 = sky_pixels[nside128][test_index]
            n64 = n128.parent_pixel
            n32 = n64.parent_pixel
            n16 = n32.parent_pixel
            n8 = n16.parent_pixel
            n4 = n8.parent_pixel
            n2 = n4.parent_pixel

            lon0 = 180.0
            fig = plt.figure(figsize=(10, 10), dpi=300)
            ax = fig.add_subplot(111)
            m = Basemap(projection='moll', lon_0=lon0)

            n2.plot(m, ax, facecolor='blue', edgecolor='blue', linewidth=0.5)
            n4.plot(m, ax, facecolor='red', edgecolor='red', linewidth=0.5)
            n8.plot(m, ax, facecolor='green', edgecolor='green', linewidth=0.5)
            n16.plot(m, ax, facecolor='magenta', edgecolor='magenta', linewidth=0.5)
            n32.plot(m, ax, facecolor='orange', edgecolor='orange', linewidth=0.5)
            n64.plot(m, ax, facecolor='pink', edgecolor='pink', linewidth=0.5)
            n128.plot(m, ax, facecolor='cyan', edgecolor='cyan', linewidth=0.5)

            meridians = np.arange(0., 360., 60.)
            m.drawparallels(np.arange(-90., 91., 30.), fontsize=14, labels=[True, True, False, False], dashes=[2, 2],
                            linewidth=0.5, xoffset=2500000)
            m.drawmeridians(meridians, labels=[False, False, False, False], dashes=[2, 2], linewidth=0.5)

            for mer in meridians[1:]:  # np.str(mer)
                plt.annotate("%0.0f" % mer, xy=m(mer, 0), xycoords='data', fontsize=14, zorder=9999)

            for axis in ['top', 'bottom', 'left', 'right']:
                ax.spines[axis].set_linewidth(2.0)

            ax.invert_xaxis()

            fig.savefig('%s/Pixel_Relation_Check.png' % diag_output_dir, bbox_inches='tight')  # ,dpi=840
            plt.close('all')

            print("... Done.")

        if self.options.plot_MWE:
            print("Plotting MWE E(B-V) extinction...")

            ebv = None
            print("\tLoading EBV...")
            with open('./pickles/ebv.pkl', 'rb') as handle:
                ebv = pickle.load(handle)

            print("E(B-V) data length: %s" % len(ebv))
            min_ebv = np.min(ebv)
            max_ebv = np.max(ebv)

            lon0 = 180.0
            fig = plt.figure(figsize=(10, 10), dpi=300)
            ax = fig.add_subplot(111)
            m = Basemap(projection='moll', lon_0=lon0)

            norm = colors.LogNorm(min_ebv, max_ebv)
            for i, (pi, pe) in enumerate(sky_pixels[nside128].items()):
                pe.plot(m, ax, facecolor=plt.cm.inferno(norm(ebv[i])), edgecolor=plt.cm.inferno(norm(ebv[i])),
                        linewidth=0.5, alpha=0.8, zorder=9900)

            meridians = np.arange(0., 360., 60.)
            m.drawparallels(np.arange(-90., 91., 30.), fontsize=14, labels=[True, True, False, False], dashes=[2, 2],
                            linewidth=0.5, xoffset=2500000)
            m.drawmeridians(meridians, labels=[False, False, False, False], dashes=[2, 2], linewidth=0.5)

            for mer in meridians[1:]:
                plt.annotate("%0.0f" % mer, xy=m(mer, 0), xycoords='data', fontsize=14, zorder=9999)

            sm = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.inferno)
            sm.set_array([])  # can be an empty list

            tks = np.logspace(np.log10(min_ebv), np.log10(max_ebv), 11)
            tks_strings = []
            for t in tks:
                tks_strings.append('%0.2f' % t)

            cb = fig.colorbar(sm, ax=ax, ticks=tks, orientation='horizontal', fraction=0.08951, pad=0.02, alpha=0.80)
            cb.ax.set_xticklabels(tks_strings, fontsize=10)
            cb.set_label("E(B-V) per pixel", fontsize=14, labelpad=10.0)
            cb.outline.set_linewidth(1.0)

            for axis in ['top', 'bottom', 'left', 'right']:
                ax.spines[axis].set_linewidth(2.0)
            ax.invert_xaxis()

            fig.savefig("%s/MWE_Debug_Plot.png" % diag_output_dir, bbox_inches='tight')
            plt.close('all')

        if self.options.plot_completeness:
            sky_distance_select = '''
                SELECT 
                    id,
                    D1,
                    D2,
                    NSIDE
                FROM SkyDistance
                ORDER BY D1;
            '''

            sky_completeness_select = '''
                SELECT
                    sp.Pixel_Index,
                    sc.SmoothedCompleteness 
                FROM SkyCompleteness sc
                JOIN SkyPixel sp on sp.id = sc.SkyPixel_id
                WHERE SkyDistance_id = %s
                ORDER By sp.Pixel_Index;
            '''

            sky_distances = batch_query([sky_distance_select])[0]

            j = 0
            this_nside = nside2
            for sd_index, d in enumerate(sky_distances):

                curr_iter = sd_index + 1

                print("\n**** Creating Completeness Slice %s/%s ****\n" % (curr_iter, len(sky_distances)))
                completeness_start = time.time()

                sd_id = d[0]
                sd_d1 = float(d[1])
                sd_d2 = float(d[2])
                sd_nside = int(d[3])
                completeness_result = query_db([sky_completeness_select % sd_id])[0]

                pixels = []
                for cr in completeness_result:
                    pix_ind = int(cr[0])
                    pix_val = float(cr[1])
                    pixels.append(Pixel_Element(pix_ind, sd_nside, pix_val))

                # check for output directory
                completeness_output_dir = "./%s/smoothed_completeness_plots/" % diag_output_dir
                if not os.path.exists(completeness_output_dir):
                    os.makedirs(completeness_output_dir)

                file_base_name = completeness_output_dir + '/sky_completeness_%s.png'
                fname = file_base_name % str(j).zfill(3)

                print("\tPlot [%0.2f to %0.2f] ..." % (sd_d1, sd_d2))
                plot_completeness(sd_d1, sd_d2, pixels, pixels, fname)

                next_nside = int(sky_distances[sd_index + 1][3])
                if (sd_index+1) < len(sky_distances) and next_nside > this_nside:
                    j += 1
                    fname = file_base_name % str(j).zfill(3)
                    this_nside = next_nside

                    outline_pixels = []
                    for sp_key, sp_val in sky_pixels[this_nside].items():
                        outline_pixels.append(sp_val)

                    plot_completeness(sd_d1, sd_d2, pixels, outline_pixels, fname)

                j += 1

                completeness_end = time.time()
                completeness_duration = (completeness_end - completeness_start)
                print("\n**** Done with Completeness Slice %s/%s - time: %s [sec] ****\n" %
                      (curr_iter, len(sky_distances), completeness_duration))

if __name__ == "__main__":
    useagestring = """python TeglonDiagnosticPlots.py [options]

Example non-box command line:
python TeglonDiagnosticPlots.py --plot_pix_relations --plot_MWE --plot_completeness

All args default to FALSE. You must pass a flag to do the operation. Note that --plot_completeness takes ~ 60 minutes to 
complete,
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
    print("Teglon `TeglonDiagnosticPlots` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")