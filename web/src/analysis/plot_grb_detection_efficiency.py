import matplotlib
matplotlib.use("Agg")
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
from scipy.interpolate import interp1d, interp2d
import scipy as sp
import matplotlib.patheffects as path_effects
print(sp.__version__)

import scipy.ndimage as ndimage
from astropy.table import Table
import csv
from collections import OrderedDict
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import math


def find_nearest_index(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


fong_2015_median_opening_angle = 16.0 # degrees
wu_2018_opening_angle = 5.0 # degrees

# additional_E_iso_scale = 2.0/(1.0 - np.cos(np.radians(wu_2018_opening_angle)/2.0))
# additional_E_iso_scale = 2.0/(1.0 - np.cos(np.radians(fong_2015_median_opening_angle)/2.0))

# DEBUG
additional_E_iso_scale = 1.0

grb_axis = "onaxis"
# grb_axis = "offaxis"

is_poster_plot = False

# Find matches to GRB170817A
sss17a_sgrb_props = (2.5, 0.3) # E_iso, n -- E_iso in units of 1e51 ergs
closest_sep = 9999
closest_model_props = None
closest_sub_dir = ""


# gw_id = "S190814bv"
gw_id = "S190425z"

# Unpack unordered results
E_keys = []
all_n = []

E_n = OrderedDict()
E_prob = OrderedDict()
all_probs = []
# for i in range(10):
for i in range(27):
    sub_dir = i+1

    results_table = Table.read("./web/events/%s/model_detection/grb/%s/Detection_Results_SGRB_%s_%s.prob" %
                               (gw_id, grb_axis, grb_axis, sub_dir), format='ascii.ecsv')

    # E0 in units of 1e50 ergs... Convert to FOE with additional scaling to E_iso
    E_result = np.asarray(list(results_table['E'])) * additional_E_iso_scale / 10.0
    n_result = list(results_table['n'])
    prob_result = list(results_table['Prob'])


    for j, E in enumerate(E_result):

        if E not in E_n:
            E_n[E] = []
        if E not in E_prob:
            E_prob[E] = []
        if E not in E_keys:
            E_keys.append(E)

        if True:
            all_n.append(n_result[j])
            all_probs.append(prob_result[j])

            E_n[E].append(n_result[j])
            E_prob[E].append(prob_result[j])


            model_sep = np.sqrt((sss17a_sgrb_props[0] - E_result[j])**2 +
                                 (sss17a_sgrb_props[1] - n_result[j])**2)
            if model_sep < closest_sep:
                closest_sep = model_sep
                closest_model_props = [E_result[j]*10, n_result[j]]
                closest_sub_dir = sub_dir

print(len(all_probs))
print("%s closest model" % grb_axis)
print("closest sep: %s" % closest_sep)
print("closest model: %s" % closest_model_props)
print("closest subdir: %s" % closest_sub_dir)


n_rows = OrderedDict()
prob_rows = OrderedDict()
for E, n_list in E_n.items():

    sorted_indices = np.argsort(E_n[E]) # sorted indices of n_list
    n_sorted = np.asarray(E_n[E])[sorted_indices]
    prob_sorted = np.asarray(E_prob[E])[sorted_indices]

    n_rows[E] = n_sorted
    prob_rows[E] = prob_sorted

Sorted_E = sorted(E_keys)
Sorted_E_n = OrderedDict()
Sorted_E_prob = OrderedDict()

for E in Sorted_E:
    Sorted_E_n[E] = n_rows[E]
    Sorted_E_prob[E] = prob_rows[E]

# print(min(all_n))
# print(max(all_n))
# print(min(E_keys))
# print(max(E_keys))


# From the data, get a log-spaced sampling of the ranges
arr_len = 100
model_E_input = np.logspace(np.min(np.log10(E_keys)) + 0.000001, np.max(np.log10(E_keys)) - 0.000001, arr_len)
model_n_input = np.logspace(np.min(np.log10(all_n)), np.max(np.log10(all_n)), arr_len)
print("Range in n: %s, %s" % (np.min(all_n), np.max(all_n)))
print("Range in E: %s, %s" % (np.min(E_keys), np.max(E_keys)))


# Interpolate rows
Interp_E_prob = OrderedDict()
for E in Sorted_E:
    n = Sorted_E_n[E]

    p = prob_rows[E]
    test_f = interp1d(n, p, kind="slinear") # , fill_value="extrapolate"
    Interp_E_prob[E] = test_f(model_n_input)

# # Stand-alone interpolate columns
# Interp_n_prob = OrderedDict()
# nlist = Sorted_E_n[E_keys[0]]
# for i, n in enumerate(nlist):
#     p = []
#     for E in E_keys:
#         p.append(Sorted_E_prob[E][i])
#
#     test_f = interp1d(E_keys, p, kind="slinear")
#     Interp_n_prob[n] = test_f(model_E_input)

# row-dependent interpolate columns
Interp_n_prob = OrderedDict()
for i, n in enumerate(model_n_input):

    column_probs = []
    sorted_Es = []

    for j, E in enumerate(Interp_E_prob.keys()):
        sorted_Es.append(E)
        column_probs.append(Interp_E_prob[E][i])

    test_f = interp1d(sorted_Es, column_probs, kind="slinear")
    interp_probs = test_f(model_E_input)
    Interp_n_prob[n] = interp_probs


# FLATTEN THE GRID
model_tuples = []
for i, (n, prob_col) in enumerate(Interp_n_prob.items()):
    for j, p in enumerate(prob_col):
        model_tuples.append((n, model_E_input[j], p))

X = []
Y = []
Z = []
for mt in model_tuples:
    X.append(mt[0])
    Y.append(mt[1])
    Z.append(mt[2])

# CHECK INITIAL PROB RANGE
print("min/max prob from files: %s, %s" % (np.min(Z), np.max(Z)))

# INTERPOLATE IN LOGSPACE
_2d_func = interp2d(np.log10(X), np.log10(Y), Z, kind="linear")
# min_plot_prob = _2d_func(min(X), min(Y))
# GET INTERPOLATED VALUES IN LOGSPACE
z_new = _2d_func(np.log10(model_n_input), np.log10(model_E_input))

# INTERPOLATED PROB RANGE (SHOULD MATCH INITIAL)
print("min/max prob from interp: %s, %s" % (np.min(z_new), np.max(z_new)))

# # GRB 170817A - convert Troja et al. 2017 to FOE using 13 degrees and Eq 1 from Wu 2018
# GRB170817A_E_k_iso_FOE = (3e+50 * 2/(1-np.cos(np.radians(13.0)/2.0))) / 1e+51
#
# # GRB 170817A - density from Troja et al. 2017
# GRB170817A_n = 0.07

# GRB 170817A - Ari (https://iopscience.iop.org/article/10.3847/2041-8213/aa91b3), Section 4.2 and Fig. 5
GRB170817A_E_k_iso_FOE = 2.5e+51 / 1e+51
# GRB 170817A - density Ari
GRB170817A_n = 3e-1



print("\n\nDone with interp\n\n")


min_prob_orig = np.min(all_probs)
max_prob_orig = np.max(all_probs)
print("Min Orig Prob: %s" % min_prob_orig)
print("Max Orig Prob: %s" % max_prob_orig)

min_prob = np.min(z_new)
max_prob = np.max(z_new)
print("Min Interp Prob: %s" % min_prob)
print("Max Interp Prob: %s" % max_prob)


min_prob = min_prob_orig
max_prob = max_prob_orig

# Convert to percent for colorbar...
norm = colors.LogNorm(min_prob_orig * 100, max_prob_orig * 100)
# test_norm = colors.LogNorm(min_plot_prob, max_prob_orig)
# norm = colors.LogNorm(min_plot_prob * 100, max_prob_orig * 100)
tks = np.logspace(np.log10(min_prob_orig * 100), np.log10(max_prob_orig * 100), 4)
# tks = np.logspace(np.log10(min_plot_prob[0] * 100), np.log10(max_prob_orig * 100), 4)
print(tks)

fig = plt.figure(figsize=(10,10), dpi=600)
ax = fig.add_subplot(111)


print("plotting...")

xx, yy = np.meshgrid(model_n_input, model_E_input)
test = ndimage.gaussian_filter(z_new, sigma=1.0, order=0)


# DEBUG
do_squares = False
do_no_interp = True
do_1D_interp = False

if do_squares:
    if do_no_interp:

        for E in Sorted_E:
            n = Sorted_E_n[E]
            p = prob_rows[E]

            for i in range(len(n)):
                n_i = n[i]
                p_i = p[i]
                clr = plt.cm.inferno(norm(p_i * 100))
                ax.plot(n_i, E, 's', color=clr, markeredgecolor=clr, markersize=17.5)


    else:
        if do_1D_interp:
            for i, n in enumerate(model_n_input):
                n_probs = Interp_n_prob[n]
                for j, prob in enumerate(n_probs):
                    clr = plt.cm.inferno(norm(prob * 100))
                    ax.plot(n, model_E_input[j], 's', color=clr, markeredgecolor=clr, markersize=17.5)
        else: # use 2D interp
            for n in model_n_input:
                for E in model_E_input:
                    p = _2d_func(n, E)
                    clr = plt.cm.inferno(norm(p * 100))
                    ax.plot(n, E, 's', color=clr, markeredgecolor=clr, markersize=17.5)
else:
    print("plotting contours with 2D interpolation")
    # DATA GRID IN LINEAR SPACE (PLOT AXES WILL BE LOG SCALED)
    from matplotlib import ticker
    CS_filled = ax.contourf(xx, yy, test, cmap=plt.cm.inferno, locator=ticker.LogLocator(),
                            levels=np.logspace(np.log10(min_prob), np.log10(max_prob), 1000))
                            # , levels=np.logspace(np.log10(min_prob), np.log10(max_prob), 1000))
    # CS_filled = ax.contourf(xx, yy, test, cmap=plt.cm.viridis, levels=np.logspace(np.log10(min_prob), np.log10(max_prob), 200))


if grb_axis == "onaxis":

    # 0425 on-axis
    CS_lines90 = ax.contour(xx, yy, test, colors="mediumturquoise", levels=[0.0, 0.10], linewidths=2.0)
    plt.setp(CS_lines90.collections, path_effects=[path_effects.withStroke(linewidth=3.0, foreground='black')])

    ##### GRB 170817A #####

    # 0814
    # fmt_dict = {
    #     0.1: "10%",
    #     0.5: "50%",
    #     0.90: "90%"
    # }
    # manual_locations1 = [(1e-1, 1e-2)]
    # 0425
    fmt_dict = {
        0.0: "",
        0.10: "10%"
    }
    clbls90 = ax.clabel(CS_lines90, inline=True, fontsize=24, fmt=fmt_dict)  # inline_spacing=130,
    plt.setp(clbls90, path_effects=[path_effects.withStroke(linewidth=3.0, foreground='black')], zorder=9900)

    # for c in clbls90:
    #     c.set_rotation(20)

    # # 0814
    # manual_locations2 = [(2e-1, 2e-3)]
    # clbls50 = ax.clabel(CS_lines50, inline=True, inline_spacing=250, fontsize=24, fmt=fmt_dict,
    #                     manual=manual_locations2)
    #
    # for c in clbls50:
    #     c.set_rotation(20)
    #
    # plt.setp(clbls50, path_effects=[path_effects.withStroke(linewidth=1.0, foreground='black')], zorder=9900)
    # plt.setp(clbls90, path_effects=[path_effects.withStroke(linewidth=1.0, foreground='black')], zorder=9900)

    if not is_poster_plot:
        pass
        # ax.text(2.0e-4, 2e-3, r'$\mathrm{\theta_{obs}=0\degree}$', fontsize=32, color="white",
        #         path_effects=[path_effects.withStroke(linewidth=2.0, foreground='black')], zorder=9999)

elif grb_axis == "offaxis":

    # 0814
    # fmt_dict = {
    #     0.1:"10%",
    #     0.5: "50%",
    #     0.90: "90%"
    # }
    # 0425
    fmt_dict = {
        0.10: "10%",
        0.20: "20%",
        0.30: "30%",
    }
    # manual_locations90 = [(2e-4, 6e-1)]
    # manual_locations50 = [(1e-3, 1.1e-1)]
    # manual_locations10 = [(2e-2, 2e-2)]

    # manual_locations90 = [(1e-2, 1.5e-1)]

    # # 0814
    # manual_locations90 = [(1e-2, 1.0e-1)]
    # manual_locations50 = [(4e-2, 4e-2)]
    # manual_locations10 = [(2e-1, 8e-3)]
    #
    # # clbls90 = ax.clabel(CS_lines90, inline=True, inline_spacing=-150, fontsize=24, fmt=fmt_dict,
    # #                     manual=manual_locations90)
    #
    # clbls90 = ax.clabel(CS_lines90, inline=True, inline_spacing=-150, fontsize=24, fmt=fmt_dict,
    #                     manual=manual_locations90)
    #
    # clbls50 = ax.clabel(CS_lines50, inline=True, inline_spacing=-50, fontsize=24, fmt=fmt_dict,
    #                     manual=manual_locations50)
    #
    # clbls10 = ax.clabel(CS_lines10, inline=True, inline_spacing=350, fontsize=24, fmt=fmt_dict,
    #                     manual=manual_locations10)

    # for c5, c10, c90 in zip(clbls50, clbls10, clbls90):
    #     c5.set_rotation(-40)
    #     c10.set_rotation(-40)
    #     c90.set_rotation(-40)
    #
    # plt.setp(clbls50, path_effects=[path_effects.withStroke(linewidth=1.0, foreground='black')], zorder=9900)
    # plt.setp(clbls10, path_effects=[path_effects.withStroke(linewidth=1.0, foreground='black')], zorder=9900)
    # plt.setp(clbls90, path_effects=[path_effects.withStroke(linewidth=1.0, foreground='black')], zorder=9900)



    # bbox = dict(boxstyle='square,pad=0') use_clabeltext=True,
    # plt.setp(clbls, path_effects=[path_effects.withStroke(linewidth=1.0, foreground='black')], zorder=9900,
    #          rotation=-40) # bbox=dict(boxstyle='square,pad=0', fc='none', ec='none')

    # manual_locations = [(1e-4, 1.1), (3e-4, 1.0), (6e-4, 7e-1), (1e-3, 6e-1), (4e-3, 1e-1)]
    # ax.clabel(CS_lines, inline=False, fontsize=24, fmt="%0.1f")
    # ax.text(2e-2, 3e-2, "20%", fontsize=24, color='red', rotation=-45)  # rotation is counter-clockwise
    # ax.text(1e-3, 1.1e-1, "60%", fontsize=24, color='red', rotation=-45)  # rotation is counter-clockwise
    # ax.text(1e-4, 5e-1, "90%", fontsize=24, color='red', rotation=-45)  # rotation is counter-clockwise

    # 0425
    # clbls10 = ax.clabel(CS_lines10, inline=True, inline_spacing=-175, fontsize=24, fmt=fmt_dict)  # inline_spacing=130,
    # clbls20 = ax.clabel(CS_lines20, inline=True, inline_spacing=-173, fontsize=24, fmt=fmt_dict)  # inline_spacing=130,
    # clbls30 = ax.clabel(CS_lines30, inline=True, inline_spacing=-155, fontsize=24, fmt=fmt_dict)  # inline_spacing=130,
    # for c in clbls90:
    #     c.set_rotation(20)
    # plt.setp(clbls10, path_effects=[path_effects.withStroke(linewidth=1.0, foreground='black')], zorder=9900)
    # plt.setp(clbls20, path_effects=[path_effects.withStroke(linewidth=1.0, foreground='black')], zorder=9900)
    # plt.setp(clbls30, path_effects=[path_effects.withStroke(linewidth=1.0, foreground='black')], zorder=9900)

    # ax.text(2.0e-4, 2e-3, r'$\mathrm{\theta_{obs}=17\degree}$', fontsize=32, color="white",
    #     path_effects=[path_effects.withStroke(linewidth=2.0, foreground='black')], zorder=9999)

    # ax.annotate(r'$\mathrm{\theta_{obs}=17\degree}$', xy=(3.0e-5, 2e-3), xycoords="data", fontsize=40, color='white')



# load GRBs
grbs = []
class sGRB:
    def __init__(self, name, Eiso, Eiso_err_upper, Eiso_err_lower, n, n_err_upper, n_err_lower):
        self.name = name
        self.Eiso = Eiso
        self.Eiso_err_upper = Eiso_err_upper
        self.Eiso_err_lower = Eiso_err_lower
        self.n = n
        self.n_err_upper = n_err_upper
        self.n_err_lower = n_err_lower


with open("./web/models/grb/fong_2015_table_3.txt", 'r') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',')
    next(csvreader)  # skip header

    for row in csvreader:
        name = row[0]

        # Eiso from Fong et al 2015 in units of 10^52 ergs -- convert to FOE
        eiso = float(row[1])*10.0
        eiso_err_upper = float(row[2])*10.0
        eiso_err_lower = float(row[3])*10.0

        n = float(row[4])
        n_err_upper = float(row[5])
        n_err_lower = float(row[6])

        grbs.append(sGRB(name, eiso, eiso_err_upper, eiso_err_lower, n, n_err_upper, n_err_lower))

# def make_error_boxes(ax, xdata, ydata, xerror, yerror, facecolor='k',
#                      edgecolor='k', alpha=0.15):
#
#     xerror = [[g.n_err_lower], [g.n_err_upper]]
#     yerror = [[g.Eiso_err_lower], [g.Eiso_err_upper]]
#
#     rect = Rectangle(xy=(xdata-xerror[0][0], ydata-yerror[0][0]), width=(xerror[0][0] + xerror[1][0]),
#                      height=(yerror[0][0] + yerror[1][0]))
#
#     # Create patch collection with specified colour/alpha
#     pc = PatchCollection([rect], facecolor=facecolor, alpha=alpha,
#                          edgecolor=edgecolor, linestyle=":")
#
#     # Add collection to axes
#     ax.add_collection(pc)
#
#     # Plot errorbars
#     artists = ax.errorbar(xdata, ydata, xerr=xerror, yerr=yerror,
#                           fmt='None', ecolor='white') # k
#     artists[2][0].set_path_effects([path_effects.Stroke(linewidth=3.0, foreground="black"), path_effects.Normal()])
#     artists[2][1].set_path_effects([path_effects.Stroke(linewidth=3.0, foreground="black"), path_effects.Normal()])
#
#
#
#     ax.errorbar(x=xdata, y=ydata, fmt='D', mfc='white', mec="None", ms=4)
#
#     return artists

for gi, g in enumerate(grbs):
    xerror = [[g.n_err_lower], [g.n_err_upper]]
    yerror = [[g.Eiso_err_lower], [g.Eiso_err_upper]]

    # Omit GRBs outside of our window...
    if (g.n - xerror[0][0]) <= 1.1e-6 or \
            (g.n + xerror[1][0]) >= 1.0 or \
            (g.Eiso - yerror[0][0]) <= 1e-3 or \
            (g.Eiso + yerror[1][0]) >= 10.0:
        continue

    # make_error_boxes(ax, g.n, g.Eiso, xerror, yerror)

    e = ax.errorbar(x=g.n, y=g.Eiso, xerr=xerror, yerr=yerror, fmt='D', mfc='black', ms=18,
                mec="white", ecolor="white", elinewidth=1.0, mew=1.5, capsize=5.0)

    e[1][0].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])
    e[1][1].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])
    e[1][2].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])
    e[1][3].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])

    e[2][0].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])
    e[2][1].set_path_effects([path_effects.Stroke(linewidth=5.0, foreground="black"), path_effects.Normal()])

ax.errorbar(x=1e-6, y=0, fmt='D', color='red', label="SGRBs from\nFong et al. 2015", ms=16)



##### GRB 170817A #####

if not is_poster_plot:
    ax.text(4.0e-2, 3.6, "GRB 170817A", fontsize=20, zorder=9999, color="white", #color="deepskyblue",
                path_effects = [path_effects.withStroke(linewidth=10.0, foreground='black')])
else:
    ax.text(1.6e-1, 3.5, "Short GRB", fontsize=22, ha="center", zorder=9999, color="deepskyblue",
            # color="deepskyblue",
            path_effects=[path_effects.withStroke(linewidth=2.0, foreground='black')])
ax.errorbar(x=GRB170817A_n, y=GRB170817A_E_k_iso_FOE, fmt='*', mew=1.5, mec="black", ms=24.0, color='black',
            label="GRB 170817A\nMurguia-berthier et al. 2017")


# x_test, y_test = np.log10(3e-1), np.log10(2.5)
# print("test_point_prob: %s" % _2d_func(x_test, y_test))
# e = ax.errorbar(x=x_test, y=y_test, fmt='D', mfc='red', ms=10,
#                 mec="black", ecolor="red", elinewidth=2.0, mew=1.0, capsize=5.0)
# for mt in model_tuples:
#     clr = norm(mt[2])
    # ax.errorbar(x=np.log10(mt[0]), y=np.log10(mt[1]), fmt='s', mfc=clr, ms=3,
    #             mec="black", ecolor="red", elinewidth=2.0, mew=1.0, capsize=5.0)
    # ax.plot(x=np.log10(mt[0]), y=np.log10(mt[1]), color=clr, marker="s")
    # X.append(mt[0])
    # Y.append(mt[1])
    # Z.append(mt[2])






ax.set_yscale('log')
ax.set_xscale('log')


sm = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.inferno)
sm.set_array([])  # can be an empty list

class ScalarFormatterClass(ScalarFormatter):
    def _set_format(self):
        self.format = "%2.2e"

cbformat = ScalarFormatterClass()  # create the

print("axis: %s" % grb_axis)
print("min/max prob: %s, %s" % (min_prob, max_prob))


# [1.87613862e-06 4.63174163e-04 1.14346724e-01 2.82294963e+01] # in percent
# tks_strings_on_axis = [
#     "0%",
#     "5e-4%",
#     "0.1%",
#     "28%",
#     # "0.0",
#     # "4.6e-6",
#     # "1.1e-3",
#     # "0.28",
#
# ]

#
# tks_strings_on_axis = [
#     "0%",
#     "0.4%",
#     "3.2%",
#     "28%"
# ]

tks_strings_on_axis = ["%s" % tk for tk in tks]

cb = fig.colorbar(sm, ax=ax, orientation='vertical', fraction=0.05, pad=0.02) #, alpha=0.80 , format=cbformat
cb.set_ticks(tks)
cb.ax.tick_params(length=6.0, labelsize=24) # width=2.0,

if grb_axis == "onaxis":
    cb.ax.set_yticklabels(tks_strings_on_axis, fontsize=24)
    # 2D Interpolation is on the LOGARITHM of the input...
    detect_prob = _2d_func(np.log10(GRB170817A_n), np.log10(GRB170817A_E_k_iso_FOE))
    print("detect of on-axis: %s" % detect_prob)
# else:
#     cb.ax.set_yticklabels(tks_strings_off_axis, fontsize=24)

if not is_poster_plot:
    cb.set_label("", fontsize=16, labelpad=9.0)
else:
    cb.set_label("Probability to Detect", fontsize=32, labelpad=11.0)
cb.ax.tick_params(length=8.0, width=2.0)
# cb.ax.locator_params(nbins=5)
cb.outline.set_linewidth(2.0)

ax.margins(x=0, y=0)
ax.tick_params(axis='both', which='major', labelsize=24, length=12.0, width=2)
ax.tick_params(axis='both', which='minor', labelsize=24, length=8.0, width=2)


# model_E_input2 = np.logspace(np.log10(1.0e-3), np.log10(10.0), arr_len)
# model_n_input2 = np.logspace(np.log10(1.0e-6), np.log10(1.0), arr_len)
# ax.set_xlim([1.0e-6, 1.0])
# ax.set_ylim([5.0e-2, 10.0])
# ax.set_ylim(1e-4, 1000)

if not is_poster_plot:
    plt.xlabel(r'n $\mathrm{\left(cm^{-3}\right)}$',fontsize=32)
else:
    plt.xlabel(r'Surrounding Gas Density $\mathrm{\left(cm^{-3}\right)}$', fontsize=32)

if not is_poster_plot:
    plt.ylabel(r'$\mathrm{E_{K,iso}}$ $\left(\times 10^{51} \mathrm{ergs}\right)}$',fontsize=32)
else:
    plt.ylabel(r'Burst Energy $\left(\times 10^{51} \mathrm{ergs}\right)}$', fontsize=32)
# if grb_axis == "onaxis":
#
#     ll =plt.legend(loc="upper left", framealpha=1.0, fontsize=24, borderpad=0.2, handletextpad=0.2)
#
#     # DEBUG LEGEND
#     # ll = plt.legend(loc="upper right", framealpha=1.0, fontsize=16, borderpad=0.2, handletextpad=0.2,
#     #                 bbox_to_anchor=(1.2, 1.5))
#     for text in ll.get_texts():
#         text.set_color("black")

# pe = [path_effects.withStroke(linewidth=2.0, foreground='white')]
# ax.grid(which='major', axis='both', linestyle=':', color="black", zorder=9999, path_effects=pe)


for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(2.0)

if not is_poster_plot:
    # fig.savefig('GW190814_GRB_Prob2Detect_%s.png' % grb_axis, bbox_inches='tight')
    fig.savefig('./web/events/S190425z/model_detection/grb/%s/0425_GRB_Prob2Detect_%s.png' % (grb_axis, grb_axis),
                bbox_inches='tight')
else:
    fig.savefig('./web/events/S190425z/model_detection/GW190814_GRB_poster.png', bbox_inches='tight', transparent=True)
plt.close('all')



#
#
# model_n_input = np.logspace(-6.0, 0.0, 100)
# test_f = interp1d(test_n_list, test_prob_list, kind="slinear")
# model_probs = test_f(model_n_input)
#
# fig = plt.figure(figsize=(10,10), dpi=1000)
# ax = fig.add_subplot(111)
# for i, n in enumerate(test_n_list):
#     ax.plot(n, test_prob_list[i], marker='s', markerfacecolor='None', markeredgecolor='k', markersize=17.5)
#
# for i, mp in enumerate(model_probs):
#     ax.plot(model_n_input[i], mp, 'b+', markersize=12.0)
#
#
#
#
# ax.margins(x=0,y=0)
# ax.tick_params(axis='both', which='major', labelsize=24, length=12.0, width=2)
# ax.set_xscale('log')
# plt.xlabel(r'n $\mathrm{[cm^{-3}]}$',fontsize=32)
# plt.ylabel('Probability',fontsize=32)
# fig.savefig('190814_GRB_slice.png', bbox_inches='tight')


plt.close('all')



do_debug = False
if do_debug:
    print("Plotting cuts...")
    fig = plt.figure(figsize=(20, 10), dpi=1000)



    ## DEBUG -- finding numbers for charlie
    log_model_n = np.log10(model_n_input) # X
    log_model_E = np.log10(model_E_input) # Y

    test_n = GRB170817A_n
    # test_n = 6.0e-5 # read value off plot
    fixed_n = np.full_like(log_model_E, test_n) # fixed n in the shape of E_iso
    prob_fixed_n = _2d_func(np.log10(fixed_n), log_model_E)

    test_E = GRB170817A_E_k_iso_FOE
    # test_E = 1.84e-1
    # test_E = 2.0e-2
    fixed_E = np.full_like(log_model_n, test_E) # fixed E_iso in the shape of n
    prob_fixed_E = _2d_func(log_model_n, np.log10(fixed_E))


    # given a column (i.e. an "n"), find E for 50%
    idx_E = find_nearest_index(prob_fixed_n[:, 0], 0.5)
    E_50 = model_E_input[idx_E]
    print("n = 3e-1, E = %0.3E, prob ~ 0.5" % (E_50))

    # given a row (i.e. an "E"), find n for 50%
    idx_n = find_nearest_index(prob_fixed_E[0, :], 0.5)
    print("E = 2.5, n = %0.3E, prob ~ 0.5" % (model_n_input[idx_n]))


    ax1 = fig.add_subplot(221)
    test_n_index = find_nearest_index(np.log10(model_n_input), np.log10(test_n))
    ax1.plot(np.power(10, log_model_E), z_new[:,test_n_index], label='n=%0.2E' % model_n_input[test_n_index])
    ax1.set_xlabel('E', fontsize=32)
    ax1.set_ylabel('Prob', fontsize=32)
    ax1.grid()
    ax1.set_xscale('log')
    ax1.legend()
    ax1.set_yticks(list(ax1.get_yticks()) + [0.5])
    ax1.vlines(E_50, 0, 1.0, colors='r')

    ax2 = fig.add_subplot(222)
    ax2.plot(np.power(10, log_model_E), prob_fixed_n[:,0], label='n=%0.2E' % test_n)
    ax2.set_xlabel('E', fontsize=32)
    ax2.set_ylabel('Prob', fontsize=32)
    ax2.grid()
    ax2.set_xscale('log')
    ax2.legend()
    ax2.set_yticks(list(ax2.get_yticks()) + [0.5])
    ax2.vlines(E_50, 0, 1.0, colors='r')



    ax3 = fig.add_subplot(223)
    test_E_index = find_nearest_index(np.log10(model_E_input), np.log10(test_E))
    ax3.plot(np.power(10, log_model_n), z_new[test_E_index, :], label='E=%0.2E' % model_E_input[test_E_index])
    ax3.set_xlabel('n', fontsize=32)
    ax3.set_ylabel('Prob', fontsize=32)
    ax3.grid()
    ax3.set_xscale('log')
    ax3.legend()
    if 0.4 in list(ax3.get_yticks()) and 0.6 in list(ax3.get_yticks()):
        ax3.set_yticks(list(ax3.get_yticks()) + [0.5])

    ax4 = fig.add_subplot(224)
    ax4.plot(np.power(10, log_model_n), prob_fixed_E[0,:], label='E=%0.2E' % test_E)
    ax4.set_xlabel('n', fontsize=32)
    ax4.set_ylabel('Prob', fontsize=32)
    ax4.grid()
    ax4.set_xscale('log')
    ax4.legend()

    if 0.4 in list(ax4.get_yticks()) and 0.6 in list(ax4.get_yticks()):
        ax4.set_yticks(list(ax4.get_yticks()) + [0.5])




    fig.savefig('./web/events/S190425z/model_detection/190814_GRB_test_%s.png' % grb_axis, bbox_inches='tight')



print("... Done.")