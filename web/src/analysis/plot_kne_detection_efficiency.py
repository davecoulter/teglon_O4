import matplotlib

matplotlib.use("Agg")
# matplotlib.use("TkAgg")
from matplotlib.ticker import ScalarFormatter

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.pyplot import cm
from matplotlib.patches import CirclePolygon
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
from scipy.special import erf
from scipy.optimize import minimize, minimize_scalar
import scipy.stats as st
from scipy.integrate import simps, quad
from scipy import interpolate
# from scipy.interpolate import interp1d, interp2d, spline
import scipy as sp

print(sp.__version__)

from astropy.table import Table
import csv
from collections import OrderedDict
import matplotlib.patheffects as path_effects

import scipy.ndimage as ndimage

from scipy.interpolate import griddata, interp2d
from matplotlib import ticker

import astropy.constants as c


def find_nearest_index(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


C_CGS = c.c.cgs.value
M_SUN_CGS = c.M_sun.cgs.value
G_CGS = c.G.cgs.value
R_NS_CGS = 20e5  # NS radius cm == 20 km

is_log = False
blue_kn = False  # else, red_kn
is_poster_plot = False

is_Ye_0_05 = False
is_Ye_0_10 = False
is_Ye_0_20 = False
is_Ye_0_30 = False
is_Ye_0_40 = False
is_Ye_0_45 = False

if blue_kn:
    is_Ye_0_45 = True
else:
    is_Ye_0_10 = True

ye_thresh = 0.05  # default
start_file = 0
end_file = 112

if is_Ye_0_10:
    ye_thresh = 0.10
    start_file = 112
    end_file = 224
elif is_Ye_0_20:
    ye_thresh = 0.20
    start_file = 224
    end_file = 336
elif is_Ye_0_30:
    ye_thresh = 0.30
    start_file = 336
    end_file = 448
elif is_Ye_0_40:
    ye_thresh = 0.40
    start_file = 448
    end_file = 560
elif is_Ye_0_45:
    ye_thresh = 0.45
    start_file = 560
    end_file = 672

model_tups = []

# for i in range(start_file, end_file):
#     file_num = i+1
#     results_table = Table.read("../Events/S190814bv/ModelDetection/Detection_KNe_%i.prob" % file_num,
#                                format='ascii.ecsv')
#
#     vej = list(results_table['vej'])
#     mej = list(results_table['mej'])
#     ye = list(results_table['ye'])
#     prob = list(results_table['Prob'])
#
#     for j, y in enumerate(ye):
#         if y == ye_thresh:
#             model_tups.append((vej[j], mej[j], prob[j]))

# 0817 + 0425 props
# red_kn_props = (0.15, 0.035, 0.1)
# blue_kn_props = (0.25, 0.025, 0.45)  # vej, mej, ye



# # Debug - low velocity, high mass
# blue_kn_props = (0.1, 0.3, 0.45)  # vej, mej, ye
# # Debug - low velocity, low mass
# blue_kn_props = (0.1, 0.1, 0.45)  # vej, mej, ye
# # Debug - high velocity, high mass
# blue_kn_props = (0.3, 0.3, 0.45)  # vej, mej, ye
# # Debug - high velocity, low mass
# blue_kn_props = (0.3, 0.1, 0.45)  # vej, mej, ye
#
# # Debug - low velocity, high mass
# red_kn_props = (0.1, 0.3, 0.1)  # vej, mej, ye
# # Debug - low velocity, low mass
# red_kn_props = (0.1, 0.1, 0.1)  # vej, mej, ye
# # Debug - high velocity, high mass
# red_kn_props = (0.3, 0.3, 0.1)  # vej, mej, ye
# Debug - high velocity, low mass
red_kn_props = (0.3, 0.1, 0.1)  # vej, mej, ye



closest_blue_sep = 9999
closest_red_sep = 9999
closest_blue = None
closest_red = None
closest_sub_dir = ""

# for i in range(18):

v_critical = np.sqrt(1.0 - (1.0 / ((3 * G_CGS * 0.5 * M_SUN_CGS) / (5.0 * R_NS_CGS * C_CGS ** 2.0) + 1.0)) ** 2.0)

utilized_mej = []
utilized_vej = []
test = []

utilized_prob = []
for i in range(36):
    sub_dir = i + 1
    file_num = i + 1
    f_path = "./web/events/S190425z/model_detection/kne/Detection_KNe_%i.prob" % file_num
    # f_path = "../Events/S190814bv/ModelDetection/Detection_KNe_%i.prob" % file_num
    # f_path = "../Events/S190425z/ModelDetection/Detection_KNe_%i.prob" % file_num
    # f_path = "../Events/S190425z/ModelDetection.old/Detection_KNe_%i.prob" % file_num
    results_table = Table.read(f_path, format='ascii.ecsv')

    vej = list(results_table['vej'])
    mej = list(results_table['mej'])
    ye = list(results_table['ye'])
    prob = list(results_table['Prob'])

    for j, y in enumerate(ye):

        if y == ye_thresh:
            model_tups.append((vej[j], mej[j], prob[j]))

            if mej[j] >= 0.5 and vej[j] >= v_critical:
                test.append((vej[j], mej[j], prob[j]))
                utilized_prob.append(prob[j])
                utilized_mej.append(mej[j])
                utilized_vej.append(vej[j])

            if blue_kn:
                blue_sep = np.sqrt((blue_kn_props[0] - vej[j]) ** 2 + (blue_kn_props[1] - mej[j]) ** 2 + (
                        blue_kn_props[1] - ye[j]) ** 2)
                if blue_sep < closest_blue_sep:
                    closest_blue_sep = blue_sep
                    closest_blue = [vej[j], mej[j], ye[j]]
                    closest_sub_dir = sub_dir
            else:
                red_sep = np.sqrt(
                    (red_kn_props[0] - vej[j]) ** 2 + (red_kn_props[1] - mej[j]) ** 2 + (red_kn_props[1] - ye[j]) ** 2)
                if red_sep < closest_red_sep:
                    closest_red_sep = red_sep
                    closest_red = [vej[j], mej[j], ye[j]]
                    closest_sub_dir = sub_dir

print("Max Raw Prob %s" % np.max(utilized_prob))

if blue_kn:
    print("Blue")
    print(closest_blue_sep)
    print(closest_blue)
    print(closest_sub_dir)
else:
    print("Red")
    print(closest_red_sep)
    print(closest_red)
    print(closest_sub_dir)

print("Done reading in data... %s rows" % len((model_tups)))

# Debug exit
raise Exception("Stop!")


all_vej = np.asarray([mt[0] for mt in model_tups])
all_mej = np.asarray([mt[1] for mt in model_tups])
points = [(mt[0], mt[1]) for mt in model_tups]
values = [mt[2] for mt in model_tups]

model_vej = np.logspace(np.log10(np.min(all_vej)), np.log10(np.max(all_vej)), 1000)
model_mej = np.logspace(np.log10(np.min(all_mej)), np.log10(np.max(all_mej)), 1000)
grid_vej, grid_mej = np.meshgrid(model_vej, model_mej)
# grid_vej, grid_mej = np.mgrid[model_vej[0]:model_vej[-1]:]
grid_prob = griddata(points, values, (grid_vej, grid_mej), method='linear', rescale=True)
min_prob = np.min(values)
max_prob = np.max(values)

print("Max Grid Prob %s" % max_prob)

gamma = 1.0 / np.sqrt(1.0 - model_vej ** 2)
grid_Mass_KE = grid_mej * (gamma - 1.0) * M_SUN_CGS * C_CGS ** 2  # ergs

# default to linear

# Convert to percent for colorbar...
blue_norm = colors.Normalize(0.0, max_prob * 100)
red_norm = colors.Normalize(0.0, max_prob * 100)

blue_tks = np.linspace(0.0, max_prob * 100, 6)
red_tks = np.linspace(0.0, max_prob * 100, 6)

# tks = [min_prob, 1.0e-9, 2.0e-9, 3.0e-9, 4.0e-9, 5.0e-9, max_prob]
if is_log:
    blue_norm = colors.LogNorm(min_prob, max_prob)
    red_norm = colors.LogNorm(min_prob, max_prob)

    blue_tks = np.logspace(np.log10(min_prob), np.log10(max_prob), 5)
    red_tks = np.logspace(np.log10(min_prob), np.log10(max_prob), 5)

fig = plt.figure(figsize=(10, 10), dpi=600)
ax = fig.add_subplot(111)
ax.set_xscale('log')
ax.set_yscale('log')


# sorted_tups = sorted(model_tups, key=lambda x: x[2], reverse=False)
# for mt in sorted_tups:
#     clr = plt.cm.viridis(norm(mt[2]))
#     ax.plot(mt[0], mt[1], marker='s', color=clr, markersize=5.0, alpha=1.0) # alpha=0.5


def bound_mass(beta):
    mass_Msol = ((5.0 * R_NS_CGS * C_CGS ** 2) / (3.0 * G_CGS)) * (1.0 / np.sqrt(1.0 - beta ** 2) - 1.0) / M_SUN_CGS
    return mass_Msol


b_mass = bound_mass(model_vej)
ceiling = np.full(len(b_mass), 5e-1)

# ax.fill_between(model_vej, ceiling, b_mass, zorder=9900, hatch="\\\\", edgecolor="gray", facecolor="none", linewidth=0.0)
ax.fill_between(model_vej, ceiling, b_mass, zorder=9990, color="gray")
ax.plot(model_vej, b_mass, linestyle="--", color="black", label="Binding Energy", zorder=9990)

if not is_poster_plot:
    ke_mass_lvls = [1e49, 1e50, 1e51, 1e52]  # ergs
    manual_locations = [
        (3.5e-2, 3.0e-3),  # 1e49
        (6.3e-2, 1e-2),  # 1e50
        # (1.355e-1, 8e-2), # 1e51
        (4.1e-1, 1e-2),  # 1e51
        (2.2e-1, 1.6e-1)]  # 1e52

    ke_lines = ax.contour(grid_vej, grid_mej, grid_Mass_KE, levels=ke_mass_lvls, colors='white',
                          locator=ticker.LogLocator(), zorder=8888)
    plt.setp(ke_lines.collections, path_effects=[path_effects.withStroke(linewidth=3.0, foreground='black')])
    # plt.setp(ke_lines.collections, path_effects=[path_effects.withStroke(linewidth=3.0, foreground='white')])

    ke_fmt_dict = {
        1e49: r"$10^{49}$ ergs",
        1e50: r"$10^{50}$ ergs",
        1e51: r"$10^{51}$ ergs",
        1e52: r"$10^{52}$ ergs",
        # 1e53:r"$10^{53}$ ergs"
    }

    # KNe label = 18 font; path effect linewidth = 1.0
    ke_clbls = ax.clabel(ke_lines, inline=True, fontsize=24, fmt=ke_fmt_dict, inline_spacing=10.0,
                         manual=manual_locations)
    plt.setp(ke_clbls, path_effects=[path_effects.withStroke(linewidth=3.0, foreground='black')], zorder=8888)
    # plt.setp(ke_clbls, path_effects=[path_effects.withStroke(linewidth=0.5, foreground='white')], zorder=9990)

# cnt = ax.contourf(grid_vej, grid_mej, grid_prob, cmap=plt.cm.inferno,
#                   levels=np.logspace(np.log10(min_prob), np.log10(max_prob), 1000), zorder=8800)  # , rasterized=True

cnt = ax.contourf(grid_vej, grid_mej, grid_prob, cmap=plt.cm.inferno,
                  levels=np.linspace(min_prob, max_prob, 1000), zorder=8800)  # , rasterized=True

# This corrects the aliasing for eps plots
for c in cnt.collections:
    c.set_edgecolor("face")

if blue_kn:

    # test = ndimage.gaussian_filter(grid_prob, sigma=3.0, order=0)
    test = grid_prob

    CS_lines = ax.contour(grid_vej, grid_mej, test, colors="mediumturquoise", linewidths=2.0,
                          levels=[0.0, 0.01, 0.1, 0.25], zorder=9900)  # , 0.35
    # CS_lines = ax.contour(grid_vej, grid_mej, test, colors="red", linewidths=2.0,
    #                       levels=[0.0, 0.01, 0.1, 0.25, 0.3, 0.4, 0.5, max_prob], zorder=9900)  # , 0.35
    plt.setp(CS_lines.collections, path_effects=[path_effects.withStroke(linewidth=3.0, foreground='black')])

    fmt_dict = {
        0.0: "",
        0.01: "1%",
        0.10: "10%",
        0.25: "25%",
        # 0.35:"35%"
    }
    # 0814
    # manual_locations = [
    #     (1.2e-1, 3e-2), # 1%
    #     (2.2e-1, 7e-2), # 10%
    #     (4e-1, 2e-1) # 25%
    # ]
    # 0425
    manual_locations = [
        (2.0e-1, 7e-2),  # 25%
        (1.5e-1, 1e-2),  # 10%
        (4.0e-2, 1.5e-3)  # 1%
    ]
    clbls = ax.clabel(CS_lines, inline=True, fontsize=24, fmt=fmt_dict, inline_spacing=100.0, manual=manual_locations)
    plt.setp(clbls, path_effects=[path_effects.withStroke(linewidth=3.0, foreground='black')], zorder=9900)
else:
    print("here")

    test = ndimage.gaussian_filter(grid_prob, sigma=2.8, order=0)

    # grid_prob
    CS_lines = ax.contour(grid_vej, grid_mej, test, colors="mediumturquoise", linewidths=2.0,  # zorder=9900)
                          # levels=[0.0, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6], zorder=9900)

                          # For 0814
                          # levels=[0.0, 0.4e-7, 2e-7], zorder=9900)

                          # For 0425
                          levels=[0.0, 1.0e-2], zorder=9900)

    CS_lines.set_clim(min_prob, max_prob)
    plt.setp(CS_lines.collections, path_effects=[path_effects.withStroke(linewidth=3.0, foreground='black')])

    # For 0814
    # fmt_dict = {
    #     0.0: "",
    #     # 1e-6: r"$1 \times 10^{-4}$ %",
    #     # 2e-6: r"$3 \times 10^{-4}$ %",
    #     # 5e-6: r"$6 \times 10^{-4}$ %",
    #     # 1e-12: r"$1 \times 10^{-10}$ %",
    #     # 1e-11: r"$1 \times 10^{-9}$ %",
    #     # 1e-10: r"$1 \times 10^{-8}$ %",
    #     4e-8: r"$4 \times 10^{-6}$ %",
    #     2e-7: r"$2 \times 10^{-5}$ %",
    # }

    # For 0425
    fmt_dict = {
        0.0: "",
        1e-2: r"1 %",
    }

    clbls = ax.clabel(CS_lines, inline=True, fontsize=24, fmt=fmt_dict, inline_spacing=120.0)
    plt.setp(clbls, path_effects=[path_effects.withStroke(linewidth=3.0, foreground='black')], zorder=9900)

# ax.text(8.8e-2, 2.3e-1, r"$U_{R_{20}} \geq KE$", rotation=31, fontsize=20, zorder=9999)
if not is_poster_plot:
    ax.text(6e-2, 1.2e-1, r"$U_{R_{20}} \geq KE$", rotation=45, fontsize=24, zorder=9999)
else:
    ax.text(4.0e-2, 2.0e-1, "Unphysical Values", rotation=45, fontsize=24, zorder=9999)

sm = None
sm_norm = None
if blue_kn:
    sm_norm = blue_norm
else:
    sm_norm = red_norm

sm = plt.cm.ScalarMappable(norm=sm_norm, cmap=plt.cm.inferno)
sm.set_array([])  # can be an empty list


class ScalarFormatterClass(ScalarFormatter):
    def _set_format(self):
        self.format = "%1.0f%%"


cbformat = ScalarFormatterClass()  # create the

tks = []
tks_strings = []
if blue_kn:

    # tks_strings = [
    #     "0%",
    #     "6%",
    #     "12%",
    #     "19%",
    #     "25%",
    #     "31%"
    # ]
    tks = blue_tks
else:
    # tks_strings = [
    #     "0%",
    #     "1.1%",
    #     "0.02%",
    #     "0.03%",
    #     "0.04%",
    #     "2.8%"
    # ]

    # tks_strings = ["%0.0f" % t for t in tks]

    tks = red_tks
    # cbformat.set_powerlimits((-2, 0))  # set the limits for sci. not.

# Turn this back on
# tks_strings = ["0%", "11%", "22%", "33%", "44%"]


# cb.ax.set_yticklabels(tks_strings, fontsize=24)


cb = fig.colorbar(sm, ax=ax, orientation='vertical', fraction=0.05, pad=0.02, alpha=0.80, format=cbformat)
cb.set_ticks(tks)
cb.ax.tick_params(labelsize=24)
# cb.ax.yaxis.get_offset_text().set_fontsize(24)


if blue_kn:

    if not is_poster_plot:
        cb.set_label("", fontsize=24, labelpad=-5.0)
    else:
        cb.set_label("Probability to Detect", fontsize=32, labelpad=11.0)
else:
    # For 0814
    # cb.set_label(r"$\times 10^{-4}$ %", fontsize=24, labelpad=2.0)

    # NO LABEL FOR 0425
    pass

    # cb.set_label(r"$\times 10^{-6}$ %", fontsize=24, labelpad=5.0)
cb.ax.tick_params(length=8.0, width=2.0)
# cb.ax.locator_params(nbins=5)


# ax.text(1.5e-1, 2e-4, r"$Y_e$ = " + "%0.2f" % ye_thresh, fontsize=24, color="white",
if not is_poster_plot:
    ax.text(8e-2, 2e-3, r"$Y_e$ = " + "%0.2f" % ye_thresh, fontsize=32, color="white",
            path_effects=[path_effects.withStroke(linewidth=3.5, foreground='black')], zorder=9999)

if blue_kn:
    # SSS17a - Blue KN
    #deepskyblue
    e = ax.errorbar(0.25, 0.025, fmt="*", mfc="black", mec="black", ecolor="black",
                    elinewidth=1.0, ms=24.0, zorder=9999, mew=1.5)
    print("Prob of blue KN: %s" % griddata(points, values, (0.25, 0.025)))




    # Test from Charlie
    # ax.errorbar(0.15, 0.06, fmt="*", mfc="deepskyblue", mec="black", ms=24.0, zorder=9999, mew=1.5)

    if not is_poster_plot:

        ax.text(0.23, 0.0125, "AT 2017gfo-like\nBlue Component", fontsize=20, color="white",
                # ax.text(0.23, 0.0135, "SSS17a-like\nBlue Component", fontsize=20, color="white",
                ha="center", zorder=9999,
                path_effects=[path_effects.withStroke(linewidth=10.0, foreground='black')])
    else:
        ax.text(0.23, 0.0135, "Neutron Star\nMerger", fontsize=20, color="white",
                ha="center", zorder=9999,
                path_effects=[path_effects.withStroke(linewidth=2.0, foreground='black')])

else:
    # An r-Process Kilonova Associated with the Short-Hard GRB 130603B
    # https://arxiv.org/abs/1306.3960
    # M
    # ej ≈ 0.03 - 0.08
    # M ☉ for v ej ≈ 0.1-0.3c.

    e = ax.errorbar(0.2, 0.055, xerr=0.1, yerr=0.025, fmt="^", mfc="black", mec="white", ecolor="white", elinewidth=1.0,
                    capsize=5.0, ms=18.0, zorder=9999, mew=1.5)  #

    e[1][0].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])
    e[1][1].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])
    e[1][2].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])
    e[1][3].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])

    e[2][0].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])
    e[2][1].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])

    ax.text(0.31, 0.023, "GRB 130603B\nAfterglow Excess", fontsize=20, color="white",
            ha="center", zorder=9999, path_effects=[path_effects.withStroke(linewidth=10.0, foreground='black')])

    # SSS17a - Red KN
    # From Charlie's SSS17a paper:
    # A good description to the panchromatic data is obtained by summing a red kilonova component with
    # Mej = 0.035 ± 0.15 M⊙, vk = 0.15 ± 0.03 c, log(Xlan) = −2.0 ± 0.5,
    e = ax.errorbar(0.15, 0.035, xerr=0.03, yerr=0.015, fmt="*", mfc="black", mec="white", ms=24.0, zorder=9999, mew=1.5,
                    ecolor="white", elinewidth=1.0, capsize=5.0)
    print("Prob of red KN: %s" % griddata(points, values, (0.15, 0.035)))

    e[1][0].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])
    e[1][1].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])
    e[1][2].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])
    e[1][3].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])

    e[2][0].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])
    e[2][1].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])

    ax.text(0.17, 0.012, "AT 2017gfo-like\nRed Component", fontsize=20, ha="center", zorder=9999, color="white",
            path_effects=[path_effects.withStroke(linewidth=10.0, foreground='black')])
    # ax.text(0.17, 0.012, "SSS17a-like\nRed Component", fontsize=20, ha="center", zorder=9999, color="white",
    #         path_effects = [path_effects.withStroke(linewidth=2.0, foreground='black')])
    # path_effects = [path_effects.withStroke(linewidth=1.25, foreground='white')]

# ax.grid(which='both', axis='both', linestyle=':', color="gray", zorder=1)

for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(2.0)
    ax.spines[axis].set_zorder(9999)

cb.outline.set_linewidth(2.0)

# ax.set_zorder(9999)
ax.tick_params(axis='both', which='major', labelsize=24, length=12.0, width=2)
ax.tick_params(axis='both', which='minor', labelsize=24, length=8.0, width=2)

# ax.set_ylim([1e-4, 5e-1])
ax.set_ylim([1e-3, 5e-1])

if not is_poster_plot:
    ax.set_ylabel(r'$\mathrm{M_{ej}}$ ($\mathrm{M_{\odot}}$)', fontsize=32, labelpad=9.0)
else:
    ax.set_ylabel('Ejecta Mass (Solar Masses)', fontsize=32, labelpad=9.0)

if not is_poster_plot:
    ax.set_xlabel(r'$\mathrm{\beta_{ej}}$', fontsize=32, labelpad=-10.0)
else:
    ax.set_xlabel('Ejecta Velocity (% c)', fontsize=32, labelpad=4.0)

# model_vej = np.logspace(np.log10(np.min(all_vej)), np.log10(np.max(all_vej)), 1000)
# model_mej = np.logspace(np.log10(np.min(all_mej)), np.log10(np.max(all_mej)), 1000)
# grid_vej, grid_mej = np.meshgrid(model_vej, model_mej)
#
# points = [(mt[0], mt[1]) for mt in model_tups]
# values = [mt[2] for mt in model_tups]
# grid_prob = griddata(points, values, (grid_vej, grid_mej), method='linear', rescale=True)
#
#
#
vej_25 = 0.25
vej_max = 0.5

mej_30 = 0.3
mej_max = 0.5
#
# fixed_vej = np.full_like(log_xx, GRB170817A_E_k_iso_FOE) # fixed E_iso in the shape of n
# fixed_mej = np.full_like(log_yy, GRB170817A_n) # fixed n in the shape of E_iso
#
# prob_fixed_E = _2d_func(log_xx, fixed_E)
# prob_fixed_n = _2d_func(fixed_n, log_yy)

## Get numbers for Charlie:
prob1 = griddata(points, values, (vej_25, mej_max), method='linear', rescale=True)
prob2 = griddata(points, values, (vej_25, mej_30), method='linear', rescale=True)
prob3 = griddata(points, values, (vej_max, 0.025), method='linear', rescale=True)

print(prob1)
print(prob2)
print(prob3)

test = 1

if not is_poster_plot:
    fig.savefig('./web/events/S190425z/model_detection/kne/0425_KNe_Prob2Detect_ye_%0.3f.png' % ye_thresh,
                bbox_inches='tight')
# else:
#     fig.savefig('poster_KNE_Sensitivity.png', bbox_inches='tight', transparent=True)
# fig.savefig('new_KNE_Sensitivity_Test_ye_%0.3f.eps' % ye_thresh, bbox_inches='tight')
plt.close('all')
print("... Done.")
