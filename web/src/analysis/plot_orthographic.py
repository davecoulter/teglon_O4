import os.path

import matplotlib as mpl
print(mpl.__version__)
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
from matplotlib import colors

from web.src.objects.Detector import *
from web.src.objects.Tile import *
from web.src.utilities.Database_Helpers import *
from web.src.utilities.filesystem_utilties import *

start = time.time()


plot_dust = True
plot_contours = True
plot_tiles = True
plot_cutout = True
plot_pixels = True
plot_sun = True


nside128 = 128
# Set up external resource directories
utilities_base_dir = "/app/web/src/utilities"
pickle_output_dir = "%s/pickles/" % utilities_base_dir

gw_id = "S190425z"
healpix_file = "GW190425_PublicationSamples_flattened.fits.gz"


fig = plt.figure(figsize=(8, 8), dpi=600)
axes = {
    'east': fig.add_subplot(221, projection='astro degrees globe', center='232.0d +24.0d'),
    # 'west': fig.add_subplot(224, projection='astro degrees globe', center='75.0d -28d')
    'west': plt.axes(
        [0.125, 0.13, 0.35, 0.35], # for position 3
        projection='astro degrees globe', center='75.0d -28d')
}



for ax_name, ax in axes.items():
    ax.grid(color='black', linestyle=':', linewidth=0.4)

    # modified_labels = []
    # for item in ax.get_yticklabels():
    #     lbl = item.get_text()
    #     if "60" in lbl:
    #         lbl = ""
    #     modified_labels.append(lbl)
    #
    # ax.set_xticklabels(modified_labels)



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


# Plot dust
if plot_dust:
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

    for ax_name, ax in axes.items():
        ax.contour_hpx(ebv, colors='dodgerblue', levels=[0.5, np.max(ebv)], linewidths=0.35, alpha=1.0)
        ax.contourf_hpx(ebv, cmap='Blues', levels=[0.5, np.max(ebv)], alpha=0.25)

if plot_contours:
    for ax_name, ax in axes.items():
        _2d_contours = ax.contour_hpx(_90_50_levels_2d, colors=['None', 'black', 'black'], levels=[0.0, 0.5, 0.9],
                                      linewidths=0.60, alpha=1.0, zorder=9999)

        plt.setp(_2d_contours.collections, path_effects=[path_effects.withStroke(linewidth=1.0, foreground='black',
                                                                                 alpha=0.1)], zorder=9999)

# region SELECT Tiles
# detector_select = '''
#             SELECT
#                 d.id,
#                 d.Name,
#                 d.ST_AsText(Poly)
#             FROM Detector d
#             JOIN DetectorPlot dp on dp.Detector_id = d.id
#             JOIN ObservedTile ot on ot.Detector_id = d.id
#             WHERE ot.HealpixMap_id = 1;
#         '''
detector_select = '''
    SELECT 
        DISTINCT d.id, 
        d.Name, 
        ST_AsText(d.Poly),
        dp.Display_Name,
        dp.Color,
        dp.Size,
        dp.Linewidth,
        dp.Alpha,
        dp.zOrder
    FROM Detector d
    JOIN DetectorPlot dp on dp.Detector_id = d.id;
'''
detector_result = query_db([detector_select])[0]

detect_display_name_dict = { int(dr[0]): dr[3] for dr in detector_result }
tile_color_dict = { int(dr[0]): dr[4] for dr in detector_result }
tile_marker_size_dict = { int(dr[0]): int(dr[5]) for dr in detector_result }
tile_linewidth_dict = { int(dr[0]): float(dr[6]) for dr in detector_result }
tile_alpha_dict = { int(dr[0]): float(dr[7]) for dr in detector_result }
tile_zorder_dict = { int(dr[0]): int(dr[8]) for dr in detector_result }
patheffect_coefficient = 1.08

detectors = {}
for dr in detector_result:
    id = int(dr[0])
    if id in detect_display_name_dict:
        detectors[id] = Detector(detect_display_name_dict[id],
                                 Detector.get_detector_vertices_from_teglon_db(dr[2]),
                                  detector_id=id)

# Hack, don't get PS1 sky cells...
observed_tile_select1 = '''
        SELECT 
            OT.id, 
            OT.RA, 
            OT._Dec, 
            OT.Detector_id,
            D.Area
        FROM ObservedTile OT
        JOIN Detector D on OT.Detector_id = D.id
        # WHERE D.id = 111
        ORDER BY D.Area DESC
    '''
ot_result1 = query_db([observed_tile_select1])[0]

# ZTF_id = 20
# KAIT_id = 6
# ZTF_approx_id = 46
ZTF_id = 74
KAIT_id = 60
ZTF_approx_id = 100

big_tile_list = []
little_tile_list = []
for ot in ot_result1:
    detector_id = int(ot[3])
    area = float(ot[4])

    if detector_id == ZTF_id:
        detector_id = ZTF_approx_id

    # if detector_id == KAIT_id:
    #     continue

    t = Tile(float(ot[1]), float(ot[2]), nside=healpix_map_nside,
             detector=detectors[detector_id], plot_color=tile_color_dict[detector_id])
    if area > 1.0:
        big_tile_list.append(t)
    else:
        little_tile_list.append(t)

# observed_tile_select2 = '''
#         SELECT
#             OT.id,
#             OT.RA,
#             OT._Dec,
#             OT.Detector_id,
#             D.Area
#         FROM ObservedTile OT
#         JOIN Detector D on OT.Detector_id = D.id
#         WHERE D.id = 111 AND FieldName LIKE '%.044%';
#     '''
# ot_result2 = query_db([observed_tile_select2])[0]
# for ot in ot_result2:
#     detector_id = 112 # Hack for Projection Cells
#
#     t = Tile(float(ot[1]), float(ot[2]), nside=healpix_map_nside,
#              detector=detectors[detector_id], plot_color=tile_color_dict[detector_id])
#
#     big_tile_list.append(t)

# endregion

x1, y1 = (0, 0)
legend_labels = []
if plot_tiles:

    ax_west = axes["west"]
    all_tile_list = big_tile_list + little_tile_list

    for t in all_tile_list:
        r, g, b = colors.to_rgb(t.plot_color)
        edg_clr = (r, g, b, 1.0)
        face_clr = (r, g, b, 0.1)
        # face_clr = "None"

        if t.detector.name not in legend_labels:
            legend_labels.append(t.detector.name)
            mkr = 's'
            ax.plot(x1, y1, marker=mkr, alpha=1.0, markeredgecolor=edg_clr,
                    markerfacecolor="None", markersize=tile_marker_size_dict[t.detector.id],
                    label=t.detector.name, linestyle='None', markeredgewidth=2.0)



    for ax_name, ax in axes.items():

        for t in big_tile_list:
            # r, g, b = colors.to_rgb(t.plot_color)
            r, g, b = colors.to_rgb(t.plot_color)
            edg_clr = (r, g, b, 1.0)
            face_clr = (r, g, b, 0.1)

            # if t.detector.id == 11:
            #     for tt in all_tile_list:
            #         if tt.detector.id == 54:  # "Pan-STARRS1_skycell":
            #             _r, _g, _b = colors.to_rgb(tt.plot_color)
            #             _edg_clr = (_r, _g, _b, 1.0)
            #             _face_clr = (_r, _g, _b, 0.1)
            #
            #             tt.plot2(ax, facecolor=_face_clr, edgecolor=_edg_clr, linewidth=0.05, alpha=0.25)

            # t.plot2(ax, facecolor=face_clr, edgecolor=edg_clr, linewidth=0.25)
            pe_stroke = tile_linewidth_dict[t.detector.id] * patheffect_coefficient
            t.plot2(ax, facecolor="None", edgecolor=edg_clr, linewidth=tile_linewidth_dict[t.detector.id],
                    alpha=tile_alpha_dict[t.detector.id], zorder=tile_zorder_dict[t.detector.id],
                    path_effects=[path_effects.withStroke(linewidth=pe_stroke, foreground='black',
                                                                             alpha=tile_alpha_dict[t.detector.id])])



        for t in little_tile_list:
            r, g, b = colors.to_rgb(t.plot_color)
            edg_clr = (r, g, b, 1.0)
            face_clr = (r, g, b, 0.1)

            # if t.detector.id == 54:  # "Pan-STARRS1_skycell":
            #     continue

            # t.plot2(ax, facecolor=face_clr, edgecolor=edg_clr, linewidth=0.25)
            pe_stroke = tile_linewidth_dict[t.detector.id] * patheffect_coefficient
            t.plot2(ax, facecolor="None", edgecolor=edg_clr, linewidth=tile_linewidth_dict[t.detector.id],
                    alpha=tile_alpha_dict[t.detector.id], zorder=tile_zorder_dict[t.detector.id],
                    path_effects=[path_effects.withStroke(linewidth=pe_stroke, foreground='black',
                                                                             alpha=tile_alpha_dict[t.detector.id])])


if plot_cutout:
    axes["inset"] = plt.axes(
        # [0.15, 0.14, 0.35, 0.35], # for position 3
        # [0.60, 0.27, 0.35, 0.35], # for mid position
        # [0.55, 0.135, 0.35, 0.35],  # for position 4
        [0.545, 0.5, 0.35, 0.35],  # for position 2
        projection='astro degrees zoom',
        center='245d +20d',
        radius='12 deg', zorder=9999)
        # radius='10 deg')

    contour_fmt_dict = {
        0.5: "50%",
        0.9: "90%"
    }

    _2d_contours = axes["inset"].contour_hpx(_90_50_levels_2d, colors=['None', 'black', 'black'],
                                             levels=[0.0, 0.5, 0.9], linewidths=0.6, alpha=1.0, zorder=9999)
    # path_effects = [path_effects.withStroke(linewidth=1.0, foreground='black')]

    plt.setp(_2d_contours.collections, path_effects=[path_effects.withStroke(linewidth=2.5, foreground='black',
                                                                             alpha=0.1)])

    clbls = axes["inset"].clabel(_2d_contours, inline=True, fontsize=9, fmt=contour_fmt_dict, inline_spacing=10.0, zorder=9999)

    if plot_pixels:

        pix_to_plot = []
        central_coord = coord.SkyCoord(245.0, 20.0, unit=(u.deg, u.deg))
        max_prob = -9999
        min_prob = 0
        # min_prob = 1e-8

        all_sky_pix_indices = np.asarray(np.arange(hp.nside2npix(healpix_map_nside)))
        # all_sky_pix_indices = np.asarray([0,1,2])
        theta, phi = hp.pix2ang(nside=healpix_map_nside, ipix=all_sky_pix_indices)
        ra = np.rad2deg(phi)
        dec = np.rad2deg(0.5 * np.pi - theta)
        all_coords = coord.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))
        coord_lookup = { ind: c for ind, c in zip(np.arange(hp.nside2npix(healpix_map_nside)), all_coords)}
        seps = central_coord.separation(all_coords).degree
        in_window_pix_index = np.where(seps < 16.5)
        allowed_pix = all_sky_pix_indices[in_window_pix_index]

        print("# Resolved coordinates: %s" % len(coord_lookup))
        print("# allowed_pix: %s" % len(allowed_pix))

        # for pix_index, pix_element in _4D_pix_dict.items():
        for pix_index, pix_element in _2D_pix_dict.items():

            if pix_index in allowed_pix:
                pix_to_plot.append(pix_element)
                # if pix_element.prob < min_prob:
                #     min_prob = pix_element.prob
                if pix_element.prob > max_prob:
                    max_prob = pix_element.prob

        print("Total Pixels to Plot: %s" % len(pix_to_plot))
        print("Min probs: %s" % min_prob)
        print("Max probs: %s" % max_prob)

        norm = colors.Normalize(min_prob, max_prob)
        cb_norm = colors.Normalize(min_prob * 100, max_prob * 100)
        # norm = colors.LogNorm(min_prob, max_prob)

        class ScalarFormatterClass(ScalarFormatter):
            def _set_format(self):
                self.format = "%1.1f"
        cbformat = ScalarFormatterClass()
        cbformat.set_powerlimits((-3, 0))

        tks = np.linspace(min_prob*100, max_prob*100, 5)
        # tks_strings = ["%0.2f%%" % (t * 100) for t in tks]
        sm = plt.cm.ScalarMappable(norm=cb_norm, cmap=plt.cm.Greys)
        cb = fig.colorbar(sm, ax=axes["inset"], ticks=tks,
                          orientation='horizontal',
                          location='bottom',
                          fraction=0.025, pad=0.15, alpha=0.5,
                          aspect=38,
                          format=cbformat)
        # cb.ax.set_xticklabels(tks_strings, fontsize=12)
        cb.set_label("2D % per Pixel", labelpad=4.0) # fontsize=14,
        cb.outline.set_linewidth(1.0)

        for i, p in enumerate(pix_to_plot):
            if i % 1000 == 0:
                print("Plotted %s..." % i)

            clr = plt.cm.Greys(norm(p.prob))
            p.plot2(axes["inset"], edgecolor='None', facecolor=clr, linewidth=0.25, alpha=0.5)

    # for t in all_tile_list:
    #     if t.detector.id == 54:  # "Pan-STARRS1_skycell":
    #         t.plot2(ax, facecolor=face_clr, edgecolor=edg_clr, linewidth=0.25, alpha=1.0)

    for t in big_tile_list:
        # t.plot2(axes["inset"], facecolor=t.plot_color, edgecolor='None', alpha=0.1)
        pe_stroke = tile_linewidth_dict[t.detector.id] * patheffect_coefficient
        t.plot2(axes["inset"], facecolor='None', edgecolor=t.plot_color, linewidth=tile_linewidth_dict[t.detector.id],
                    alpha=tile_alpha_dict[t.detector.id], zorder=tile_zorder_dict[t.detector.id],
                    path_effects=[path_effects.withStroke(linewidth=pe_stroke, foreground='black',
                                                                             alpha=tile_alpha_dict[t.detector.id])])

    for t in little_tile_list:

        alpha1 = 0.1
        alpha2 = 1.0
        # if t.detector.id == 54:  # "Pan-STARRS1_skycell":
        #     alpha1 = 0.05
        #     alpha2 = 0.1
        # if t.detector.id == 54:  # "Pan-STARRS1_skycell":
        #     continue

        # t.plot2(axes["inset"], facecolor=t.plot_color, edgecolor='None', alpha=alpha1)
        pe_stroke = tile_linewidth_dict[t.detector.id] * patheffect_coefficient
        t.plot2(axes["inset"], facecolor='None', edgecolor=t.plot_color, linewidth=tile_linewidth_dict[t.detector.id],
                    alpha=tile_alpha_dict[t.detector.id], zorder=tile_zorder_dict[t.detector.id],
                    path_effects=[path_effects.withStroke(linewidth=pe_stroke, foreground='black',
                                                                             alpha=tile_alpha_dict[t.detector.id])])

    axes["east"].mark_inset_axes(axes["inset"], zorder=9999)
    # axes["east"].mark_inset_circle(axes["inset"], '245d +20d', '20 deg')
    # axes["east"].connect_inset_circle(axes["inset"], '245d +20d', '12 deg')

    # from matplotlib.patches import ConnectionPatch
    # xy = (x[i], y[i])
    # con = ConnectionPatch(xyA=xy, xyB=xy, coordsA="data", coordsB="data",
    #                       axesA=ax2, axesB=ax1, color="red")

    # Why doesn't this code work?
    # axes["east"].connect_inset_axes(axes["inset"], 'upper left')
    # axes["east"].connect_inset_axes(axes["inset"], 'lower left')

    axes["inset"].grid(color='black', linestyle=':', linewidth=0.25)
    axes["inset"].tick_params(axis='both', colors='black')



    for ax in [axes["inset"]]:
        ax.set_facecolor('none')
        for key in ['ra', 'dec']:
            ax.coords[key].set_auto_axislabel(False)

if plot_sun:

    ax = axes["west"]
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
    twi18 = ax.contour_hpx(deg_sep, colors='gold', levels=[0, 45.0], alpha=1.0, linewidths=0.35)
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
    c_labels = ax.clabel(twi18, inline=True, fontsize=9, fmt=fmt,
                         levels=twi18.levels, colors=('gold', 'gold', 'gold'), use_clabeltext=True,
                         inline_spacing=60) # , zorder=9999
    [txt.set_path_effects([path_effects.withStroke(linewidth=0.75,
                                                   foreground='k')]) for txt in c_labels]
    plt.setp(twi18.collections, path_effects=[path_effects.withStroke(linewidth=1.0, foreground='black',
                                                                             alpha=0.1)])

axes["east"].set_ylabel('Eastern Spur', labelpad=0.5, fontsize=16)
axes["west"].set_ylabel('Western Spur', labelpad=0.5, fontsize=16)
axes["west"].legend(bbox_to_anchor=(1.25, 0.80), ncols=2, frameon=False, columnspacing=0.0, labelspacing=1.2,
                    prop={'size': 10}) # , bbox_transform=fig.transFigure



output_file = "S190425_figure_2.png"
output_path = "./web/src/analysis/%s" % output_file

fig.savefig(output_path, bbox_inches='tight', format="png")
chmod_outputfile(output_path)
plt.close('all')
print("... Done.")


end = time.time()
duration = (end - start)
print("\n********* start DEBUG ***********")
print("Teglon `plot_orthographic` execution time: %s" % duration)
print("********* end DEBUG ***********\n")