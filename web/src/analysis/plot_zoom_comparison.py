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


plot_contour = True
plot_mwe = True
plot_sun = True
plot_pix = True
all_sky_colorbar = True
is_debug = False



class ScalarFormatterClass(ScalarFormatter):
    def _set_format(self):
        self.format = "%1.1f"

nside128 = 128
# Set up external resource directories
utilities_base_dir = "/app/web/src/utilities"
pickle_output_dir = "%s/pickles/" % utilities_base_dir

gw_id = "S190425z"
healpix_file = "GW190425_PublicationSamples_flattened.fits.gz"

central_coord = coord.SkyCoord(245.0, 20.0, unit=(u.deg, u.deg))

fig = plt.figure(figsize=(8, 8), dpi=600)
axes = {
        '2D': fig.add_subplot(223, projection='astro degrees zoom',
                              center='%sd %sd' % (central_coord.ra.deg, central_coord.dec.deg), radius='12 deg'),
        '4D': fig.add_subplot(224, projection='astro degrees zoom',
                              center='%sd %sd' % (central_coord.ra.deg, central_coord.dec.deg), radius='12 deg')
}
if all_sky_colorbar:
    axes["all_sky"] = plt.axes([0.1855, 0.28, 0.65, 0.65], projection='astro hours mollweide')
else:
    axes["all_sky"] = plt.axes([0.1855, 0.31, 0.65, 0.65], projection='astro hours mollweide')


plt.rcParams["axes.axisbelow"] = False

for ax_name, ax in axes.items():
    ax.grid(color='silver', linestyle=':', linewidth=0.6, path_effects=[path_effects.withStroke(linewidth=0.8, foreground='black',
                                                                                 alpha=1.0)])

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

if plot_contour:
    contour_fmt_dict = {
        0.5: "50%",
        0.9: "90%"
    }

    axes["all_sky"].contour_hpx(_90_50_levels_2d, colors=['None', 'white', 'white'], levels=[0.0, 0.5, 0.9],
                       linewidths=0.25, alpha=1.0,
                       path_effects=[path_effects.withStroke(linewidth=0.75, foreground='black')])
    axes["all_sky"].contourf_hpx(_90_50_levels_4d, cmap=plt.cm.inferno_r, levels=np.linspace(0.0, 1.0, 50), linewidths=0.0,
                            alpha=0.9)


    axes["2D"].contourf_hpx(_90_50_levels_4d, cmap=plt.cm.inferno_r, levels=np.linspace(0.0, 1.0, 50),
                                 linewidths=0.0,
                                 alpha=0.9)
    _2d_contours = axes["2D"].contour_hpx(_90_50_levels_2d, colors=['None', 'white', 'white'], levels=[0.0, 0.5, 0.9],
                                linewidths=0.25, alpha=1.0)
    plt.setp(_2d_contours.collections, path_effects=[path_effects.withStroke(linewidth=2.5, foreground='white',
                                                                             alpha=0.1)])
    clbls = axes["2D"].clabel(_2d_contours, inline=True, fontsize=9, fmt=contour_fmt_dict, inline_spacing=10.0)



    axes["4D"].contourf_hpx(_90_50_levels_4d, cmap=plt.cm.inferno_r, levels=np.linspace(0.0, 1.0, 50),
                                 linewidths=0.0,
                                 alpha=0.9)
    _4d_contours = axes["4D"].contour_hpx(_90_50_levels_2d, colors=['None', 'white', 'white'], levels=[0.0, 0.5, 0.9],
                           linewidths=0.25, alpha=1.0)
    plt.setp(_4d_contours.collections, path_effects=[path_effects.withStroke(linewidth=2.5, foreground='white',
                                                                             alpha=0.1)])
    clbls = axes["4D"].clabel(_4d_contours, inline=True, fontsize=9, fmt=contour_fmt_dict, inline_spacing=10.0)

cbformat = ScalarFormatterClass()
norm = colors.Normalize(min(map_4d_pix) * 100, max(map_4d_pix) * 100)
tks = np.linspace(min(map_4d_pix) * 100, max(map_4d_pix) * 100, 5)
cbformat.set_powerlimits((-1, 0))  # set the limits for sci. not.
lon = axes["all_sky"].coords[0]
lat = axes["all_sky"].coords[1]
lon.set_ticklabel(color='silver', size=12,
                  path_effects=[path_effects.withStroke(linewidth=0.75, foreground='black')], zorder=9999)
lat.set_ticklabel(color='black', size=12, zorder=9999)
if all_sky_colorbar:
    sm = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.inferno)
    sm.set_array([])  # can be an empty list
    cb = fig.colorbar(sm, ax=axes["all_sky"], ticks=tks, orientation='horizontal',
                      fraction=0.08951, pad=0.08, alpha=0.80,
                      aspect=40, format=cbformat,
                      location='top')

    cb.set_label("% per Pixel", fontsize=14, labelpad=4.0)
    cb.outline.set_linewidth(1.0)
axes["all_sky"].set_ylabel(r'$\mathrm{Declination}$', fontsize=16, labelpad=0.5)
axes["all_sky"].set_xlabel(r'$\mathrm{Right\;Ascension}$', fontsize=16, labelpad=-8.5)

axes["all_sky"].mark_inset_axes(axes["2D"], path_effects=[path_effects.withStroke(linewidth=1.5, foreground='white',
                                                                                 alpha=1.0)])

axes["all_sky"].connect_inset_axes(axes["2D"], 'upper left',
                           path_effects=[path_effects.withStroke(linewidth=1.5, foreground='white',
                                                                                 alpha=1.0)]) # color="black",
axes["all_sky"].connect_inset_axes(axes["4D"], 'upper right',
                           path_effects=[path_effects.withStroke(linewidth=1.5, foreground='white',
                                                                                 alpha=1.0)])



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

    axes["all_sky"].contour_hpx(ebv, colors='dodgerblue', levels=[0.5, np.max(ebv)], linewidths=0.35, alpha=1.0)
    axes["all_sky"].contourf_hpx(ebv, cmap='Blues', levels=[0.5, np.max(ebv)], alpha=0.25)

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

    twi18 = axes["all_sky"].contour_hpx(deg_sep, colors='gold', levels=[0, 45.0], alpha=1.0, linewidths=0.35)
    axes["all_sky"].contourf_hpx(deg_sep, colors='cornsilk', levels=[0, 45.0], linewidths=0.2, alpha=0.25)  # alpha=0.25

    axes["all_sky"].plot(sun_coord.ra.degree, sun_coord.dec.degree, transform=axes["all_sky"].get_transform('world'),
                         marker="o", markersize=8,
            markeredgecolor="darkorange", markerfacecolor="gold", linewidth=0.05)

    fmt = {}
    strs = ['', r'$45\degree$']
    for l, s in zip(twi18.levels, strs):
        fmt[l] = s
    c_labels = axes["all_sky"].clabel(twi18, inline=True, fontsize=8, fmt=fmt,
                         levels=twi18.levels, colors=('gold', 'gold', 'gold'), use_clabeltext=True,
                         inline_spacing=60) # , zorder=9999
    [txt.set_path_effects([path_effects.withStroke(linewidth=0.75,
                                                   foreground='k')]) for txt in c_labels]
    plt.setp(twi18.collections, path_effects=[path_effects.withStroke(linewidth=1.0, foreground='black',
                                                                             alpha=0.1)])



def plot_pixels_within_radius(ax, allowed_pix_indices, pixel_collection, is_log, label_text):

    pix_to_plot = []

    max_prob = -9999
    min_prob = 9999
    if is_log:
        min_prob = 1e-8

    print("# Resolved coordinates: %s" % len(coord_lookup))
    print("# allowed_pix: %s" % len(allowed_pix))

    for pix_index, pix_element in pixel_collection.items():

        if pix_index in allowed_pix_indices:
            pix_to_plot.append(pix_element)

            if not is_log:
                if pix_element.prob < min_prob:
                    min_prob = pix_element.prob
            if pix_element.prob > max_prob:
                max_prob = pix_element.prob

    print("Total Pixels to Plot: %s" % len(pix_to_plot))
    print("Min probs: %s" % min_prob)
    print("Max probs: %s" % max_prob)

    norm = colors.Normalize(min_prob, max_prob)
    alpha_norm = colors.LogNorm(1e-8, max_prob)
    cb_norm = colors.Normalize(min_prob * 100, max_prob * 100)
    if is_log:
        norm = colors.LogNorm(min_prob, max_prob)
        cb_norm = colors.LogNorm(min_prob * 100, max_prob * 100)

    cbformat = ScalarFormatterClass()
    cbformat.set_powerlimits((-1, 0))

    tks = np.linspace(min_prob*100, max_prob*100, 5)
    if is_log:
        tks = np.linspace(np.log10(min_prob * 100), np.log10(max_prob * 100), 5)
    sm = plt.cm.ScalarMappable(norm=cb_norm, cmap=plt.cm.inferno)
    cb = fig.colorbar(sm, ax=ax, ticks=tks,
                      orientation='horizontal',
                      location='bottom',
                      fraction=0.025, pad=0.15, alpha=1.0,
                      aspect=38,
                      format=cbformat)
    cb.set_label(label_text, labelpad=4.0) # fontsize=14,
    cb.outline.set_linewidth(1.0)

    for i, p in enumerate(pix_to_plot):
        if i % 1000 == 0:
            print("Plotted %s..." % i)

        clr = plt.cm.inferno(norm(p.prob))

        thresh = 0.15e-4
        alph = alpha_norm(thresh)
        if p.prob > thresh:
            alph = alpha_norm(p.prob)
        p.plot2(ax, edgecolor='None', facecolor=clr, linewidth=0.1, alpha=alph)

    for _ax in [ax]:
        _ax.set_facecolor('none')
        for key in ['ra', 'dec']:
            _ax.coords[key].set_auto_axislabel(False)

radius = 16.5
all_sky_pix_indices = np.asarray(np.arange(hp.nside2npix(healpix_map_nside)))

if is_debug:
    all_sky_pix_indices = all_sky_pix_indices[0:10] # DEBUG

theta, phi = hp.pix2ang(nside=healpix_map_nside, ipix=all_sky_pix_indices)
ra = np.rad2deg(phi)
dec = np.rad2deg(0.5 * np.pi - theta)
all_coords = coord.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))
coord_lookup = { ind: c for ind, c in zip(np.arange(hp.nside2npix(healpix_map_nside)), all_coords)}
seps = central_coord.separation(all_coords).degree
in_window_pix_index = np.where(seps < radius)
allowed_pix = all_sky_pix_indices[in_window_pix_index]

if plot_pix:
    plot_pixels_within_radius(axes['2D'], allowed_pix, _2D_pix_dict, is_log=False, label_text='Original % Per Pixel')
    plot_pixels_within_radius(axes['4D'], allowed_pix, _4D_pix_dict, is_log=False, label_text='Resampled % Per Pixel')

output_file = "S190425_figure_4.png"
output_path = "./web/src/analysis/%s" % output_file

fig.savefig(output_path, bbox_inches='tight', format="png")
chmod_outputfile(output_path)
plt.close('all')
print("... Done.")

end = time.time()
duration = (end - start)
print("\n********* start DEBUG ***********")
print("Teglon `plot_zoom_comparison` execution time: %s" % duration)
print("********* end DEBUG ***********\n")