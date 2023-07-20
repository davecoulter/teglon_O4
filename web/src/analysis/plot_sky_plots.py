import os.path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patheffects as path_effects
import copy
from matplotlib.ticker import ScalarFormatter
mpl.use("Agg")

from collections import OrderedDict

from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy.coordinates import get_sun
from astroplan import Observer
from dustmaps.sfd import SFDQuery
from ligo.skymap.postprocess.util import find_greedy_credible_levels
import ligo.skymap.plot
import pickle

from web.src.objects.Detector import *
from web.src.objects.Tile import *
from web.src.utilities.Database_Helpers import *
from web.src.utilities.filesystem_utilties import *

start = time.time()



plot_2D = False
plot_4D = True
plot_sun = True
plot_mwe = True

plot_tiles = False
plot_east = False
plot_west = False

plot_big = True
plot_small = True

plot_pixels_by_epoch = False
plot_pixels_by_depth = False
plot_pixels_by_delta_t = False
plot_pixels_by_filters = False


nside128 = 128
# Set up external resource directories
utilities_base_dir = "/app/web/src/utilities"
pickle_output_dir = "%s/pickles/" % utilities_base_dir

gw_id = "S190425z"
healpix_file = "GW190425_PublicationSamples_flattened.fits.gz"

healpix_map_select = "SELECT id, RescaledNSIDE, t_0 FROM HealpixMap WHERE GWID = '%s' and Filename = '%s'"
healpix_map_result = query_db([healpix_map_select % (gw_id, healpix_file)])[0][0]
healpix_map_id = int(healpix_map_result[0])
healpix_map_nside = int(healpix_map_result[1])
healpix_map_event_time = Time(healpix_map_result[2], format="gps")
time_of_trigger = Time(healpix_map_event_time.to_datetime(), scale="utc")


pix_2D_pkl = "./web/src/analysis/2d_pix.pkl"
pix_4D_pkl = "./web/src/analysis/4d_pix.pkl"
_2D_pix_dict = None
_4D_pix_dict = None

if not os.path.exists(pix_2D_pkl):
    print("Building pixel pickles...")

    _2D_pix_dict = OrderedDict()
    _4D_pix_dict = OrderedDict()

    select_pix = '''
        SELECT
            hp.id,
            hp.HealpixMap_id,
            hp.Pixel_Index,
            hp.Prob,
            hpc.NetPixelProb,
            hp.Distmu,
            hp.Distsigma,
            hp.Mean,
            hp.Stddev,
            hp.Norm,
            hp.N128_SkyPixel_id
        FROM HealpixPixel hp
        JOIN HealpixPixel_Completeness hpc on hpc.HealpixPixel_id = hp.id
        WHERE hp.HealpixMap_id = %s
        ORDER BY hp.Pixel_Index;
    '''

    pixels_to_select = select_pix % healpix_map_id
    map_pix_result = query_db([pixels_to_select])[0]

    for mpr in map_pix_result:

        pix_id = int(mpr[0])
        pix_index = int(mpr[2])
        pix_2d = float(mpr[3])
        pix_4d = float(mpr[4])
        mean_dist = float(mpr[7])
        stddev_dist = float(mpr[8])

        if pix_index not in _2D_pix_dict:
            _2D_pix_dict[pix_index] = Pixel_Element(pix_index, healpix_map_nside, pix_2d, pix_id, mean_dist, stddev_dist)
            _4D_pix_dict[pix_index] = Pixel_Element(pix_index, healpix_map_nside, pix_4d, pix_id, mean_dist, stddev_dist)

    with open(pix_2D_pkl, 'wb') as handle:
        pickle.dump(_2D_pix_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(pix_4D_pkl, 'wb') as handle:
        pickle.dump(_4D_pix_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
else:
    print("Loading pixel pickles...")

    with open(pix_2D_pkl, 'rb') as handle:
        _2D_pix_dict = pickle.load(handle)

    with open(pix_4D_pkl, 'rb') as handle:
        _4D_pix_dict = pickle.load(handle)


map_2d_pix = [float(pix.prob) for pix_i, pix in _2D_pix_dict.items()]
map_4d_pix = [float(pix.prob) for pix_i, pix in _4D_pix_dict.items()]



_90_50_levels_2d = find_greedy_credible_levels(np.asarray(map_2d_pix))
_90_50_levels_4d = find_greedy_credible_levels(np.asarray(map_4d_pix))


if not (plot_4D or plot_2D):
    # plt.rc('axes', edgecolor='white')
    if not plot_tiles:
        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'astro degrees globe',
                                                 'center': '247.5d +30d'}, dpi=600)
    elif plot_tiles:
        if plot_east:
            fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'astro degrees globe',
                                                     # 'center': '247.5d +30d'}, figsize=(4,4), dpi=600)
                                                     'center': '232.0d +24.0d'}, figsize=(4, 4), dpi=600)
        if plot_west:
            fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'astro degrees globe',
                                                     'center': '75.0d -28d'}, figsize=(4, 4), dpi=600)

    # ax.tick_params(axis='both', colors='white')

else:
    plt.rcParams["axes.axisbelow"] = False
    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'astro hours mollweide'}, figsize=(8,8), dpi=600)

# ax.grid(color='gray', linestyle=':', linewidth=0.15)
ax.grid(color='silver', linestyle=':', linewidth=0.15)


global_alpha = 1.0

if plot_2D or plot_4D:

    class ScalarFormatterClass(ScalarFormatter):
        def _set_format(self):
            self.format = "%1.1f"
    cbformat = ScalarFormatterClass()  # create the formatter

    ax.contour_hpx(_90_50_levels_2d, colors=['None', 'white', 'white'], levels=[0.0, 0.5, 0.9],
                   linewidths=0.25, alpha=1.0,
                   path_effects=[path_effects.withStroke(linewidth=0.75, foreground='black')])
    norm = None
    tks = None

    if plot_2D:
        subplot_type = "all_sky_2D"
        ax.contourf_hpx(_90_50_levels_2d, cmap=plt.cm.inferno_r, levels=np.linspace(0.0, 1.0, 50), linewidths=0.0,
                        alpha=0.9)

        norm = colors.Normalize(min(map_2d_pix) * 100, max(map_2d_pix) * 100)
        tks = np.linspace(min(map_2d_pix) * 100, max(map_2d_pix) * 100, 5)
        cbformat.set_powerlimits((-3, 0))  # set the limits for sci. not.

    if plot_4D:

        subplot_type = "all_sky_4D"
        ax.contourf_hpx(_90_50_levels_4d, cmap=plt.cm.inferno_r, levels=np.linspace(0.0, 1.0, 50), linewidths=0.0,
                        alpha=0.9)

        norm = colors.Normalize(min(map_4d_pix) * 100, max(map_4d_pix) * 100)
        tks = np.linspace(min(map_4d_pix) * 100, max(map_4d_pix) * 100, 5)
        cbformat.set_powerlimits((-1, 0))  # set the limits for sci. not.

    lon = ax.coords[0]
    lat = ax.coords[1]
    lon.set_ticklabel(color='silver', size=12,
                      path_effects=[path_effects.withStroke(linewidth=0.75, foreground='black')], zorder=9999)
    lat.set_ticklabel(color='black', size=12, zorder=9999)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.inferno)
    sm.set_array([])  # can be an empty list

    cb = fig.colorbar(sm, ax=ax, ticks=tks, orientation='horizontal',
                      fraction=0.08951, pad=0.08, alpha=0.80,
                      aspect=40, format=cbformat)

    cb.set_label("% per Pixel", fontsize=14, labelpad=4.0)
    cb.outline.set_linewidth(1.0)

    plt.ylabel(r'$\mathrm{Declination}$', fontsize=16, labelpad=0.5)
    plt.xlabel(r'$\mathrm{Right\;Ascension}$', fontsize=16, labelpad=6.5)

if plot_mwe:
    # Dust
    print("\tLoading existing mwe...")
    ebv = None
    if os.path.exists(pickle_output_dir + "ebv.pkl"):
        with open(pickle_output_dir + "ebv.pkl", 'rb') as handle:
            ebv = pickle.load(handle)
    else:
        print('ebv.pkl does not exist! Creating...')
        theta, phi = hp.pix2ang(nside=nside128, ipix=np.arange(hp.nside2npix(nside128)))
        pix_coord = coord.SkyCoord(ra=np.rad2deg(phi), dec=np.rad2deg(0.5 * np.pi - theta), unit=(u.deg, u.deg))

        print("Retrieving dust info")
        sfd = SFDQuery()
        ebv = sfd(pix_coord)

        with open(pickle_output_dir + 'ebv.pkl', 'wb') as handle:
            pickle.dump(ebv, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # mwe_norm = colors.Normalize(np.min(ebv), 2.0)

    # ax.contourf_hpx(ebv, cmap='Blues', levels=np.linspace(0.5, 2.0, 50), linewidths=0.2, alpha=0.5)
    # ax.contourf_hpx(ebv, cmap='Blues', levels=np.logspace(np.log10(0.5), np.log10(np.max(ebv)), 50), linewidths=0.2, alpha=0.5)
    # ax.contourf_hpx(ebv, cmap='Blues', levels=np.linspace(np.min(ebv), np.max(ebv), 200), linewidths=0.2, alpha=1.0)
    ax.contour_hpx(ebv, colors='dodgerblue', levels=[0.5, np.max(ebv)], linewidths=0.35, alpha=global_alpha)
    ax.contourf_hpx(ebv, cmap='Blues', levels=[0.5, np.max(ebv)], alpha=0.25)

    # fmt = {}
    # strs = ['', r'$A_{r} = 0.5 mag$']
    # for l, s in zip(mw_contour.levels, strs):
    #     fmt[l] = s
    # c_labels = ax.clabel(mw_contour, inline=True, fontsize=8, fmt=fmt,
    #                      levels=mw_contour.levels, colors='dodgerblue', use_clabeltext=True,
    #                      inline_spacing=60)
    # [txt.set_path_effects([path_effects.withStroke(linewidth=0.75,
    #                                                foreground='k')]) for txt in c_labels]

if plot_tiles:

    subplot_type = "all_tiles_{hemisphere}_{size}"

    if plot_west:
        subplot_type = subplot_type.replace("{hemisphere}", "west")
    else:
        subplot_type = subplot_type.replace("{hemisphere}", "east")

    if plot_small:
        subplot_type = subplot_type.replace("{size}", "small")
    else:
        subplot_type = subplot_type.replace("{size}", "big")

    # ax.contour_hpx(_90_50_levels_2d, colors=['None', 'white', 'white'], levels=[0.0, 0.5, 0.9],
    #                linewidths=0.25, alpha=1.0,
    #                path_effects=[path_effects.withStroke(linewidth=0.75, foreground='black')])

    # ax.contourf_hpx(_90_50_levels_4d, cmap=plt.cm.Greys_r, levels=np.linspace(0.0, 1.0, 50), linewidths=0.0,
    #                 alpha=0.9)

    # ax.contourf_hpx(_90_50_levels_2d, cmap=plt.cm.Greys_r, levels=np.linspace(0.0, 1.0, 500), linewidths=0.0,
    #                 alpha=0.25)

    _2d_contours = ax.contour_hpx(_90_50_levels_2d, colors=['None', 'black', 'black'], levels=[0.0, 0.5, 0.9],
                   linewidths=0.25, alpha=1.0,
                   path_effects=[path_effects.withStroke(linewidth=0.75, foreground='black')])

    plt.setp(_2d_contours.collections, path_effects=[path_effects.withStroke(linewidth=1.0, foreground='black',
                                                                             alpha=0.1)])


    detector_select = '''
        SELECT id, Name, ST_AsText(Poly) FROM Detector; 
    '''
    detector_result = query_db([detector_select])[0]
    detectors = {int(dr[0]): Detector(dr[1], Detector.get_detector_vertices_from_teglon_db(dr[2]),
                                      detector_id=int(dr[0])) for dr in detector_result}
    observed_tile_select = '''
        SELECT 
            OT.id, 
            OT.RA, 
            OT._Dec, 
            OT.Detector_id,
            D.Area 
        FROM ObservedTile OT
        JOIN Detector D on OT.Detector_id = D.id
        ORDER BY D.Area DESC;
    '''
    ot_result = query_db([observed_tile_select])[0]

    tile_color_dict = {
        1: "black", # Swope
        2: "magenta", # ANDICAM-CCD
        3: "red", # THACHER
        4: "limegreen", # NICKEL
        7: "gold", #LCOGT
        11: "turquoise", # CSS
        12: "violet", # Swift
        13: "brown", # MMT
        46: "royalblue", # ZTF
        45: "orange", # GOTO-4
        51: "slateblue" # ANDICAM-IR
    }

    ZTF_id = 20
    KAIT_id = 6
    ZTF_approx_id = 46

    big_tile_list = []
    little_tile_list = []

    for ot in ot_result:
        detector_id = int(ot[3])
        area = float(ot[4])

        if detector_id == ZTF_id:
            detector_id = ZTF_approx_id

        if detector_id == KAIT_id:
            continue

        t = Tile(float(ot[1]), float(ot[2]), nside=healpix_map_nside,
             detector=detectors[detector_id], plot_color=tile_color_dict[detector_id])

        if area > 1.0:
            big_tile_list.append(t)
        else:
            little_tile_list.append(t)

    if plot_big:
        for t in big_tile_list:
            t.plot2(ax, facecolor=t.plot_color, edgecolor='None', alpha=0.1)
            t.plot2(ax, facecolor='None', edgecolor=t.plot_color, linewidth=0.35, alpha=1.0)
    if plot_small:
        for t in little_tile_list:
            t.plot2(ax, facecolor=t.plot_color, edgecolor='None', alpha=0.1)
            t.plot2(ax, facecolor='None', edgecolor=t.plot_color, linewidth=0.35, alpha=1.0)


if plot_pixels_by_epoch:

    subplot_type = "visits"

    pix_stats_pkl = "./web/src/analysis/pix_stats.pkl"
    pix_stat_result = None
    if not os.path.exists(pix_stats_pkl):
        print("Selecting pixel stats...")

        select_pixel_stats = '''
            SELECT
                hp.Pixel_Index,
                ot.Mag_Lim,
                B.Name,
                (ot.MJD - 58598.346134259256) as delta_t,
                ot.Mag_Lim - (5*LOG10(hp.Mean*1e6)-5) + SPE.EBV*B.F99_Coefficient as Extincted_Abs_Mag_Lim,
                hp.Mean as pix_distance,
                SPE.EBV*B.F99_Coefficient as A_lambda,
                hp.Prob as 2D_prob,
                HPC.NetPixelProb as 4D_prob,
                D.Name
            FROM ObservedTile_HealpixPixel ot_hp
            JOIN HealpixPixel hp ON hp.id = ot_hp.HealpixPixel_id
            JOIN HealpixPixel_Completeness HPC on ot_hp.HealpixPixel_id = HPC.HealpixPixel_id
            JOIN ObservedTile ot on ot_hp.ObservedTile_id = ot.id
            JOIN SkyPixel SP on hp.N128_SkyPixel_id = SP.id
            JOIN SkyPixel_EBV SPE on SP.id = SPE.N128_SkyPixel_id
            JOIN Band B on ot.Band_id = B.id
            JOIN Detector D on ot.Detector_id = D.id
            ORDER BY Pixel_Index, delta_t ASC, Extincted_Abs_Mag_Lim ASC, ot.Mag_Lim DESC, B.Name;
        '''
        pix_stat_result = query_db([select_pixel_stats])[0]

        with open(pix_stats_pkl, 'wb') as handle:
            pickle.dump(pix_stat_result, handle, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        print("Loading pixel stats...")

        with open(pix_stats_pkl, 'rb') as handle:
            pix_stat_result = pickle.load(handle)

    # calculate pixels by number of visits
    visits = {}
    for psr in pix_stat_result:
        p_index = int(psr[0])

        if p_index not in visits:
            visits[p_index] = 0.0

        visits[p_index] += 1.0

    pix_visits = {}
    max_visit = 0
    for pix_ind, visits in visits.items():
        if visits > max_visit:
            max_visit = visits

        pix_visits[pix_ind] = Pixel_Element(pix_ind, healpix_map_nside, prob=visits)
    print("Max visits: %s" % max_visit)

    # norm = colors.Normalize(1.0, max_visit)
    norm = colors.LogNorm(1.0, max_visit)
    for pix_ind, pix in pix_visits.items():
        clr = plt.cm.inferno(norm(pix.prob))
        pix.plot2(ax, facecolor=clr, linewidth=0.05, alpha=1.0)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.inferno)
    sm.set_array([])  # can be an empty list
    tks = np.logspace(np.log10(1.0), np.log10(max_visit), 5)
    # tks = np.logspace(np.log10(-worst_depth), np.log10(-best_depth), 5)
    tks_strings = []
    for t in tks:
        # tks_strings.append('%0.1f' % (-t))
        tks_strings.append('%0.0f' % (t))

    cb = fig.colorbar(sm, ax=ax, ticks=tks, orientation='horizontal',
                      fraction=0.025, pad=0.05, alpha=0.80,
                      aspect=40)
    cb.ax.set_xticklabels(tks_strings, fontsize=14, color="white")
    cb.set_label("# Epochs", fontsize=14, labelpad=4.0, color="white")
    cb.outline.set_linewidth(1.0)
    cb.ax.xaxis.set_tick_params(which='both', color='white')

    # ax.spines['bottom'].set_color('white')
    # ax.spines['top'].set_color('white')
    # ax.spines['right'].set_color('white')
    # ax.spines['left'].set_color('white')
    for _ax in [ax]:
        ax.set_facecolor('white')
        for key in ['ra', 'dec']:
            ax.coords[key].set_auto_axislabel(False)

    # [v.set_color('white') for k, v in ax.spines.items()]
    #
    # ax.spines['bottom'].set_color('red')
    # ax.spines['top'].set_color('red')
    # ax.xaxis.label.set_color('white')
    # ax.tick_params(axis='both', colors='white')

if plot_pixels_by_depth:

    subplot_type = "depth"

    pix_stats_pkl = "./web/src/analysis/pix_stats.pkl"
    pix_stat_result = None
    if not os.path.exists(pix_stats_pkl):
        print("Selecting pixel stats...")

        select_pixel_stats = '''
            SELECT
                hp.Pixel_Index,
                ot.Mag_Lim,
                B.Name,
                (ot.MJD - 58598.346134259256) as delta_t,
                # ot.Mag_Lim - (5*LOG10(hp.Mean*1e6)-5) + SPE.EBV*B.F99_Coefficient as Extincted_Abs_Mag_Lim,
                ot.Mag_Lim - SPE.EBV*B.F99_Coefficient as Extinction_Corrected_Mag_Lim,
                hp.Mean as pix_distance,
                SPE.EBV*B.F99_Coefficient as A_lambda,
                hp.Prob as 2D_prob,
                HPC.NetPixelProb as 4D_prob,
                D.Name
            FROM ObservedTile_HealpixPixel ot_hp
            JOIN HealpixPixel hp ON hp.id = ot_hp.HealpixPixel_id
            JOIN HealpixPixel_Completeness HPC on ot_hp.HealpixPixel_id = HPC.HealpixPixel_id
            JOIN ObservedTile ot on ot_hp.ObservedTile_id = ot.id
            JOIN SkyPixel SP on hp.N128_SkyPixel_id = SP.id
            JOIN SkyPixel_EBV SPE on SP.id = SPE.N128_SkyPixel_id
            JOIN Band B on ot.Band_id = B.id
            JOIN Detector D on ot.Detector_id = D.id
            WHERE SPE.EBV*B.F99_Coefficient < 0.5
            ORDER BY Pixel_Index, delta_t ASC, Extinction_Corrected_Mag_Lim ASC, ot.Mag_Lim DESC, B.Name;
        '''
        pix_stat_result = query_db([select_pixel_stats])[0]

        with open(pix_stats_pkl, 'wb') as handle:
            pickle.dump(pix_stat_result, handle, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        print("Loading pixel stats...")

        with open(pix_stats_pkl, 'rb') as handle:
            pix_stat_result = pickle.load(handle)

    # calculate pixels by number of visits
    depths = {}
    for psr in pix_stat_result:
        p_index = int(psr[0])
        ext_corr_depth = float(psr[4])

        if p_index not in depths:
            depths[p_index] = []

        depths[p_index].append(ext_corr_depth)

    pix_depths = {}
    best_depth = -99999
    worst_depth = 99999
    for pix_ind, depth_list in depths.items():
        _best_depth = max(depth_list)
        _worst_depth = min(depth_list)
        if _best_depth > best_depth:
            best_depth = _best_depth
        if _worst_depth < worst_depth:
            worst_depth = _worst_depth

        # only keep the deepest pixel
        pix_depths[pix_ind] = Pixel_Element(pix_ind, healpix_map_nside, prob=_best_depth)
        # if not pix_ind in pix_depths:
        #     pix_depths[pix_ind] = Pixel_Element(pix_ind, healpix_map_nside, prob=_best_depth)
        # else:
        #     current_depth = pix_depths[pix_ind].prob
        #     if _best_depth > current_depth:
        #         pix_depths[pix_ind] = Pixel_Element(pix_ind, healpix_map_nside, prob=_best_depth)

    print("Best depth: %s mag" % best_depth)
    print("Worst depth: %s mag" % worst_depth)

    norm = colors.Normalize(worst_depth, best_depth)
    for pix_ind, pix in pix_depths.items():
        try:
            clr = plt.cm.inferno(norm(pix.prob))
        except Exception as e:
            print(pix.prob)
            print(e)
            raise Exception("Stop!")
        pix.plot2(ax, facecolor=clr, linewidth=0.05, alpha=1.0)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.inferno)
    sm.set_array([])  # can be an empty list
    tks = np.linspace(worst_depth, best_depth, 5)
    tks_strings = ["%0.1f" % t for t in tks]
    cb = fig.colorbar(sm, ax=ax, ticks=tks, orientation='horizontal',
                      fraction=0.025, pad=0.05, alpha=0.80,
                      aspect=40)
    cb.ax.set_xticklabels(tks_strings, fontsize=14, color="white")
    cb.set_label("Depth (mag)", fontsize=14, labelpad=4.0, color="white")
    cb.outline.set_linewidth(1.0)
    cb.ax.xaxis.set_tick_params(which='both', color='white')

if plot_pixels_by_delta_t:
    subplot_type = "delta_t"

    pix_stats_pkl = "./web/src/analysis/pix_stats.pkl"
    pix_stat_result = None
    if not os.path.exists(pix_stats_pkl):
        print("Selecting pixel stats...")

        select_pixel_stats = '''
                SELECT
                    hp.Pixel_Index,
                    ot.Mag_Lim,
                    B.Name,
                    (ot.MJD - 58598.346134259256) as delta_t,
                    # ot.Mag_Lim - (5*LOG10(hp.Mean*1e6)-5) + SPE.EBV*B.F99_Coefficient as Extincted_Abs_Mag_Lim,
                    ot.Mag_Lim - SPE.EBV*B.F99_Coefficient as Extinction_Corrected_Mag_Lim,
                    hp.Mean as pix_distance,
                    SPE.EBV*B.F99_Coefficient as A_lambda,
                    hp.Prob as 2D_prob,
                    HPC.NetPixelProb as 4D_prob,
                    D.Name
                FROM ObservedTile_HealpixPixel ot_hp
                JOIN HealpixPixel hp ON hp.id = ot_hp.HealpixPixel_id
                JOIN HealpixPixel_Completeness HPC on ot_hp.HealpixPixel_id = HPC.HealpixPixel_id
                JOIN ObservedTile ot on ot_hp.ObservedTile_id = ot.id
                JOIN SkyPixel SP on hp.N128_SkyPixel_id = SP.id
                JOIN SkyPixel_EBV SPE on SP.id = SPE.N128_SkyPixel_id
                JOIN Band B on ot.Band_id = B.id
                JOIN Detector D on ot.Detector_id = D.id
                WHERE SPE.EBV*B.F99_Coefficient < 0.5
                ORDER BY Pixel_Index, delta_t ASC, Extinction_Corrected_Mag_Lim ASC, ot.Mag_Lim DESC, B.Name;
            '''
        pix_stat_result = query_db([select_pixel_stats])[0]

        with open(pix_stats_pkl, 'wb') as handle:
            pickle.dump(pix_stat_result, handle, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        print("Loading pixel stats...")

        with open(pix_stats_pkl, 'rb') as handle:
            pix_stat_result = pickle.load(handle)

    times = {}
    for psr in pix_stat_result:
        p_index = int(psr[0])
        delta_t = float(psr[3])

        if p_index not in times:
            times[p_index] = []

        times[p_index].append(delta_t)

    pix_delta_t = {}
    best_delta_t = 99999
    worst_delta_t = -99999
    for pix_ind, delta_t_list in times.items():
        _best_delta_t = min(delta_t_list)
        _worst_delta_t = max(delta_t_list)
        if _best_delta_t < best_delta_t:
            best_delta_t = _best_delta_t
        if _worst_delta_t > worst_delta_t:
            worst_delta_t = _worst_delta_t

        # only keep the deepest pixel
        pix_delta_t[pix_ind] = Pixel_Element(pix_ind, healpix_map_nside, prob=_best_delta_t)

    print("Best delta t: %s mag" % best_delta_t)
    print("Worst delta t: %s mag" % worst_delta_t)

    norm = colors.LogNorm(best_delta_t, worst_delta_t)
    for pix_ind, pix in pix_delta_t.items():
        try:
            clr = plt.cm.inferno_r(norm(pix.prob))
        except Exception as e:
            print(pix.prob)
            print(e)
            raise Exception("Stop!")
        pix.plot2(ax, facecolor=clr, linewidth=0.05, alpha=1.0)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.inferno_r)
    sm.set_array([])  # can be an empty list
    tks = np.logspace(np.log10(best_delta_t), np.log10(worst_delta_t), 5)
    tks_strings = ["%0.1f" % t for t in tks]
    cb = fig.colorbar(sm, ax=ax, ticks=tks, orientation='horizontal',
                      fraction=0.025, pad=0.05, alpha=0.80,
                      aspect=40)
    cb.ax.set_xticklabels(tks_strings, fontsize=14, color="white")
    cb.set_label(r"t-t$_{0}$ (days)", fontsize=14, labelpad=4.0, color="white")
    cb.outline.set_linewidth(1.0)
    cb.ax.xaxis.set_tick_params(which='both', color='white')

if plot_pixels_by_filters:
    subplot_type = "filters"

    pix_stats_pkl = "./web/src/analysis/pix_stats.pkl"
    pix_stat_result = None
    if not os.path.exists(pix_stats_pkl):
        print("Selecting pixel stats...")

        select_pixel_stats = '''
                    SELECT
                        hp.Pixel_Index,
                        ot.Mag_Lim,
                        B.Name,
                        (ot.MJD - 58598.346134259256) as delta_t,
                        # ot.Mag_Lim - (5*LOG10(hp.Mean*1e6)-5) + SPE.EBV*B.F99_Coefficient as Extincted_Abs_Mag_Lim,
                        ot.Mag_Lim - SPE.EBV*B.F99_Coefficient as Extinction_Corrected_Mag_Lim,
                        hp.Mean as pix_distance,
                        SPE.EBV*B.F99_Coefficient as A_lambda,
                        hp.Prob as 2D_prob,
                        HPC.NetPixelProb as 4D_prob,
                        D.Name
                    FROM ObservedTile_HealpixPixel ot_hp
                    JOIN HealpixPixel hp ON hp.id = ot_hp.HealpixPixel_id
                    JOIN HealpixPixel_Completeness HPC on ot_hp.HealpixPixel_id = HPC.HealpixPixel_id
                    JOIN ObservedTile ot on ot_hp.ObservedTile_id = ot.id
                    JOIN SkyPixel SP on hp.N128_SkyPixel_id = SP.id
                    JOIN SkyPixel_EBV SPE on SP.id = SPE.N128_SkyPixel_id
                    JOIN Band B on ot.Band_id = B.id
                    JOIN Detector D on ot.Detector_id = D.id
                    WHERE SPE.EBV*B.F99_Coefficient < 0.5
                    ORDER BY Pixel_Index, delta_t ASC, Extinction_Corrected_Mag_Lim ASC, ot.Mag_Lim DESC, B.Name;
                '''
        pix_stat_result = query_db([select_pixel_stats])[0]

        with open(pix_stats_pkl, 'wb') as handle:
            pickle.dump(pix_stat_result, handle, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        print("Loading pixel stats...")

        with open(pix_stats_pkl, 'rb') as handle:
            pix_stat_result = pickle.load(handle)

    bands = {}
    for psr in pix_stat_result:
        p_index = int(psr[0])
        band = psr[2]

        if p_index not in bands:
            bands[p_index] = []

        if band not in bands[p_index]:
            bands[p_index].append(band)

    pix_band = {}
    max_bands = -99999
    for pix_ind, band_list in bands.items():

        if len(band_list) > max_bands:
            max_bands = len(band_list)

        pix_band[pix_ind] = Pixel_Element(pix_ind, healpix_map_nside, prob=len(band_list))

    print("Most bands: %s" % max_bands)

    norm = colors.Normalize(1.0, max_bands)
    for pix_ind, pix in pix_band.items():
        try:
            clr = plt.cm.inferno(norm(pix.prob))
        except Exception as e:
            print(pix.prob)
            print(e)
            raise Exception("Stop!")
        pix.plot2(ax, facecolor=clr, linewidth=0.05, alpha=1.0)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.inferno)
    sm.set_array([])  # can be an empty list
    tks = np.linspace(1, max_bands, 7)
    tks_strings = ["%0.0f" % t for t in tks]
    cb = fig.colorbar(sm, ax=ax, ticks=tks, orientation='horizontal',
                      fraction=0.025, pad=0.05, alpha=0.80,
                      aspect=40)
    cb.ax.set_xticklabels(tks_strings, fontsize=14, color="white")
    cb.set_label("# of filters", fontsize=14, labelpad=4.0, color="white")
    cb.outline.set_linewidth(1.0)
    cb.ax.xaxis.set_tick_params(which='both', color='white')

    # ax = plt.subplots(1, 1, subplot_kw={'projection': 'astro degrees globe', 'center': '247.5d +30d'})
    ax_zoom_rect = plt.axes(
        # [0.85, 0.2, 0.6, 0.6],
        [0.8, 0.11, 0.8, 0.8],
        projection='astro degrees zoom',
        center='245d +20d',
        radius='9 deg')

    # for pix_ind, pix in pix_band.items():
    #     try:
    #         clr = plt.cm.inferno(norm(pix.prob))
    #     except Exception as e:
    #         print(pix.prob)
    #         print(e)
    #         raise Exception("Stop!")
    #     pix.plot2(ax_zoom_rect, facecolor=clr, linewidth=0.05, alpha=1.0)
    # pix_4d = {}
    # max_prob = -999999
    # min_prob = 999999
    # for psr in pix_stat_result:
    #     p_index = int(psr[0])
    #     prob = float(psr[8])
    #     if prob > max_prob:
    #         max_prob = prob
    #     if prob < min_prob:
    #         min_prob = prob
    #     if p_index not in pix_4d:
    #         pix_4d[p_index] = Pixel_Element(p_index, healpix_map_nside, prob=prob)

    # print("Len pix: %s" % len(pix_4d))
    # print("Min: %s Max: %s" % (min_prob, max_prob))

    # norm = colors.Normalize(min_prob, max_prob)
    # for pix_ind, pix in pix_4d.items():
    #     clr = plt.cm.Greys(norm(pix.prob))
    #     pix.plot2(ax_zoom_rect, facecolor=clr, linewidth=0.05, alpha=1.0)
    max_4d = np.max(map_4d_pix)
    min_4d = np.min(map_4d_pix)
    ax_zoom_rect.contourf_hpx(_90_50_levels_4d, cmap=plt.cm.Greys_r, levels=np.linspace(0.0, 1.0, 100),
                    alpha=0.3)  # linewidths=0.2,

    tks2 = np.linspace(min_4d, max_4d, 5)
    tks_strings2 = ["%0.2f%%" % (t*100) for t in tks2]
    print(tks_strings2)
    norm2 = colors.Normalize(min_4d, max_4d)
    sm2 = plt.cm.ScalarMappable(norm=norm2, cmap=plt.cm.Greys)
    cb2 = fig.colorbar(sm2, ax=ax_zoom_rect, ticks=tks2, orientation='horizontal',
                      fraction=0.025, pad=0.15, alpha=0.80,
                      aspect=40)
    cb2.ax.set_xticklabels(tks_strings2, fontsize=12, color="white")
    cb2.set_label("Teglon Prob", fontsize=14, labelpad=4.0, color="white")
    cb2.outline.set_linewidth(1.0)
    cb2.ax.xaxis.set_tick_params(which='both', color='white')


    detector_select = '''
            SELECT id, Name, ST_AsText(Poly) FROM Detector;
        '''
    detector_result = query_db([detector_select])[0]
    detectors = {int(dr[0]): Detector(dr[1], Detector.get_detector_vertices_from_teglon_db(dr[2]),
                                      detector_id=int(dr[0])) for dr in detector_result}
    observed_tile_select = '''
            SELECT id, RA, _Dec, Detector_id FROM ObservedTile WHERE Detector_id IN (1,2,3,4,51);
        '''
    _colors = {
        1: "red",
        2: "orange",
        3: "green",
        4: "blue",
        51: "brown"

    }
    ot_result = query_db([observed_tile_select])[0]
    tile_list = [Tile(float(ot[1]), float(ot[2]), nside=healpix_map_nside,
                      detector=detectors[int(ot[3])], tile_id=int(ot[3])) for ot in ot_result]
    for t in tile_list:
        t.plot2(ax_zoom_rect, facecolor="None", edgecolor=_colors[t.id], linewidth=0.5, alpha=1.0, zorder=9999)

    ax.mark_inset_axes(ax_zoom_rect)
    ax.connect_inset_axes(ax_zoom_rect, 'upper left')
    ax.connect_inset_axes(ax_zoom_rect, 'lower left')

    ax_zoom_rect.grid(color='white', linestyle=':', linewidth=0.15)
    ax_zoom_rect.tick_params(axis='both', colors='white')

    for ax in [ax_zoom_rect]:
        ax.set_facecolor('none')
        for key in ['ra', 'dec']:
            ax.coords[key].set_auto_axislabel(False)

if plot_sun:
    # Sun Contours
    time_of_trigger = Time(healpix_map_event_time.to_datetime(), scale="utc") # Sun Position
    theta, phi = hp.pix2ang(nside=64, ipix=np.arange(hp.nside2npix(64)))
    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)
    all_sky_coords = coord.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))

    detector_geography = {
        "SWOPE": EarthLocation(lat=-29.0182 * u.deg, lon=-70.6926 * u.deg, height=2402 * u.m),
        "THACHER": EarthLocation(lat=34.46479 * u.deg, lon=-121.6431 * u.deg, height=630.0 * u.m),
        "NICKEL": EarthLocation(lat=37.3414 * u.deg, lon=-121.6429 * u.deg, height=1283 * u.m),
        "T80S_T80S-Cam": EarthLocation(lat=-70.8035 * u.deg, lon=-70.8035 * u.deg, height=2207.0 * u.m),
    }

    observatory = Observer(location=detector_geography["SWOPE"])
    sun_set = observatory.sun_set_time(time_of_trigger, which='nearest', horizon=-18 * u.deg)
    sun_rise = observatory.sun_rise_time(sun_set, which='next', horizon=-18 * u.deg)
    sun_coord = get_sun(sun_set)
    hours_of_the_night = int((sun_rise - sun_set).to_value('hr'))
    time_next = sun_set

    sun_separations = sun_coord.separation(all_sky_coords)
    deg_sep = np.asarray([sp.deg for sp in sun_separations])

    # twi18 = ax.contour_hpx(deg_sep, colors='gold', levels=[0, 20.0, 40.0, 60.0], alpha=global_alpha, linewidths=0.35)
    twi18 = ax.contour_hpx(deg_sep, colors='gold', levels=[0, 45.0], alpha=global_alpha, linewidths=0.35)
    # ax.contourf_hpx(deg_sep, colors='cornsilk', levels=[0, 20.0, 40.0, 60.0], linewidths=0.2, alpha=0.25)

    if plot_tiles:
        ax.contourf_hpx(deg_sep, colors='cornsilk', levels=[0, 45.0], linewidths=0.2, alpha=1.0)
    else:
        ax.contourf_hpx(deg_sep, colors='cornsilk', levels=[0, 45.0], linewidths=0.2, alpha=0.25) # alpha=0.25

    ax.plot(sun_coord.ra.degree, sun_coord.dec.degree, transform=ax.get_transform('world'), marker="o", markersize=8,
            markeredgecolor="darkorange", markerfacecolor="gold", linewidth=0.05)

    fmt = {}
    strs = ['', r'$45\degree$']
    for l, s in zip(twi18.levels, strs):
        fmt[l] = s
    c_labels = ax.clabel(twi18, inline=True, fontsize=8, fmt=fmt,
                         levels=twi18.levels, colors=('gold', 'gold', 'gold'), use_clabeltext=True,
                         inline_spacing=60) # , zorder=9999
    [txt.set_path_effects([path_effects.withStroke(linewidth=0.75,
                                                   foreground='k')]) for txt in c_labels]
    plt.setp(twi18.collections, path_effects=[path_effects.withStroke(linewidth=1.0, foreground='black',
                                                                             alpha=0.1)])

# output_file = "S190425_all_tiles.svg"
# output_path = "./web/src/analysis/%s" % output_file
# fig.savefig(output_path, bbox_inches='tight', format="svg")
# chmod_outputfile(output_path)
# plt.close('all')
# print("... Done.")



output_file = "S190425_by_%s.png" % subplot_type
output_path = "./web/src/analysis/%s" % output_file
# fig.savefig(output_path, bbox_inches='tight', format="png", transparent=True)
fig.savefig(output_path, bbox_inches='tight', format="png")
chmod_outputfile(output_path)
plt.close('all')
print("... Done.")


end = time.time()
duration = (end - start)
print("\n********* start DEBUG ***********")
print("Teglon `plot_mollweide` execution time: %s" % duration)
print("********* end DEBUG ***********\n")