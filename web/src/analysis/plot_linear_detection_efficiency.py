import os.path

import matplotlib
matplotlib.use("Agg")

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

from astropy.table import Table
import csv
from collections import OrderedDict
import scipy.ndimage as ndimage


def find_nearest_index(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# model_X_input = np.linspace(-2.0, 10.0, 100)
# model_Y_input = np.linspace(-23.0, -11.0, 100)
# model_X_input = np.linspace(-2.0, 10.0, 100)
# model_Y_input = np.linspace(-23.0, -11.0, 100)
model_X_input = np.linspace(0.0015, 2.0, 100)
model_Y_input = np.linspace(-21.0, -12.0, 100)


# Unpack unordered results
M_list = []
M_dM = OrderedDict()
M_prob = OrderedDict()
test_prob_list = []
# for i in range(16):
# for i in range(373):


for i in range(100):
    sub_dir = i+1
    results_table = Table.read("./web/events/S190425z/model_detection/linear/Detection_Results_Linear_%s.prob" % sub_dir,
                               format='ascii.ecsv')

    M_list = list(results_table['M'])
    dm_list = list(results_table['dM'])
    prob_list = list(results_table['Prob'])
    test_prob_list += prob_list
    for j, M in enumerate(M_list):

        if M not in M_dM:
            M_dM[M] = []

        if M not in M_prob:
            M_prob[M] = []

        M_dM[M].append(dm_list[j])
        M_prob[M].append(prob_list[j])

print(np.min(test_prob_list), np.max(test_prob_list))

M_dM_probs = OrderedDict()
dMs = []
for key, value in M_dM.items():

    sorted_indices = np.argsort(M_dM[key])
    dM = np.asarray(M_dM[key])[sorted_indices]
    dMs = dM
    prob = np.asarray(M_prob[key])[sorted_indices]

    f_dM_prob = interp1d(dM, prob, kind="slinear")
    interp_probs = f_dM_prob(model_X_input)

    M_dM_probs[key] = interp_probs


column_interp = OrderedDict()
for i, dM in enumerate(model_X_input):
    dm_prob_list = []

    sorted_Ms = []
    for j, M in enumerate(sorted(M_dM_probs.keys())):
        sorted_Ms.append(M)
        dm_prob_list.append(M_dM_probs[M][i])

    f_M_prob = interp1d(sorted_Ms, dm_prob_list, kind="slinear")
    interp_probs = f_M_prob(model_Y_input)

    column_interp[dM] = interp_probs


model_tuples = []
for i, (dM, prob_col) in enumerate(column_interp.items()):
    for index, prob in enumerate(prob_col):
        model_tuples.append((dM, model_Y_input[index], prob))


X = []
Y = []
Z = []

for mt in model_tuples:

    X.append(mt[0])
    Y.append(mt[1])
    Z.append(mt[2])





z_test = ndimage.gaussian_filter(Z, sigma=5.0, order=0)
_2d_func = interp2d(X, Y, z_test, kind="linear")

# _2d_func = interp2d(X, Y, Z, kind="linear")
xx, yy = np.meshgrid(model_X_input, model_Y_input)
z_new = _2d_func(model_X_input, model_Y_input)





min_prob = np.min(z_new)
max_prob = np.max(z_new)

print(min_prob, max_prob)
norm = colors.Normalize(min_prob, max_prob)
# norm = colors.Normalize(0.0, 0.95)

fig = plt.figure(figsize=(10, 10), dpi=600)
ax = fig.add_subplot(111)

ax.invert_yaxis()

# 0814 contours
# test = ndimage.gaussian_filter(z_new, sigma=1.0, order=0)
test = z_new
# ax.contourf(xx, yy, z_new, levels=np.linspace(0.0, max_prob, 200), cmap=plt.cm.viridis)
ax.contourf(xx, yy, test, levels=np.linspace(0.0, max_prob, 200), cmap=plt.cm.inferno)

# CS = ax.contour(xx, yy, z_new, levels=[0.0, 0.10, 0.30, 0.5, 0.7, 0.9, max_prob], colors="red")
#
CS = ax.contour(xx, yy, test, levels=[0.0, 0.15, 0.25, 0.35, max_prob], linewidths=2.0, colors="mediumturquoise")
plt.setp(CS.collections, path_effects=[path_effects.withStroke(linewidth=3.0, foreground='black')])

# For Zoom out
# manual_locations = [(10.0, -17.0), (7.0, -20.0), (5.0, -21.0), (3.0, -21), (0.0, -22.0)]
# manual_locations = [(1e-1, -17.5), (2e-1, -17.0), (2e-1, -16.25), (3e-1, -16)] # settings used for 8/2020
manual_locations = [(8e-2, -16.5), (1.5e-1, -16.0), (2e-1, -15.75), (3e-1, -15.5)]  # settings used for 12/2020
fmt_dict = {
    0.0:"",
    0.15:"15%",
    0.25:"25%",
    0.35:"35%"
}
# manual=manual_locations
clbls = ax.clabel(CS, inline=True, fontsize=24, fmt=fmt_dict, inline_spacing=-100)
plt.setp(clbls, path_effects=[path_effects.withStroke(linewidth=3.0, foreground='black')], zorder=9900)

for c in clbls:
    c.set_rotation(40)

# For Zoom in
# ax.text(7e-3, -16.05, "90%", fontsize=24, color="red")
# ax.text(9e-3, -15.4, "70%", fontsize=24, color="red")
# ax.text(5e-2, -15.45, "50%", fontsize=24, color="red", rotation=35)
# ax.text(8.6e-2, -15.4, "30%", fontsize=24, color="red", rotation=55)
# ax.text(1.9e-1, -15.8, "10%", fontsize=24, color="red", rotation=55)

# for i in range(5):
#     p = CS.collections[i].get_paths()[1]
#     v = p.vertices
#     x = v[:,0]
#     y = v[:,1]
#
#     upper_bound = -19.0
#     # upper_bound = -17.8
#     # lower_bound = -15.0
#     lower_bound = -11.8
#
#     model_y1 = list(np.linspace(y.min(), upper_bound, 15))
#     # model_y2 = list(np.linspace(upper_bound, lower_bound, 8))
#     model_y2 = list(np.linspace(upper_bound, lower_bound, 12))
#     model_y3 = list(np.linspace(lower_bound, y.max(), 4))
#
#     y1 = y[np.where(y <= upper_bound)[0]]
#     x1 = x[np.where(y <= upper_bound)[0]]
#
#     y2 = y[(y <= lower_bound) & (y > upper_bound)]
#     x2 = x[(y <= lower_bound) & (y > upper_bound)]
#
#     y3 = y[np.where(y > lower_bound)[0]]
#     x3 = x[np.where(y > lower_bound)[0]]
#
#     tks_x1 = interpolate.splrep(y1, x1, s=1)
#     tks_x2 = interpolate.splrep(y2, x2, s=0.02)
#     tks_x3 = interpolate.splrep(y3, x3, s=1)
#
#     new_x1 = interpolate.splev(model_y1, tks_x1, der=0)
#     new_x2 = interpolate.splev(model_y2, tks_x2, der=0)
#     new_x3 = interpolate.splev(model_y3, tks_x3, der=0)
#
#     xx1 = list(new_x1)
#     xx2 = list(new_x2)
#     xx3 = list(new_x3)
#     xxx = []
#     xxx += xx1
#     xxx += xx2
#     xxx += xx3
#
#     yyy = []
#     yyy += model_y1
#     yyy += model_y2
#     yyy += model_y3
#
#     model_y = np.linspace(y.min(), y.max(), 1000)
#     tks_final = interpolate.splrep(yyy, xxx, s=50)
#     model_x = interpolate.splev(model_y, tks_final, der=0)
#
#     # ax.plot(test1, xnew1, 'b-', zorder=9999, linewidth=1.0)
#     # ax.plot(test2, xnew2, 'b-', zorder=9999, linewidth=1.0)
#     # ax.plot(test3, xnew3, 'g-', zorder=9999, linewidth=1.0)
#     # ax.plot(model_x, model_y, 'k-', zorder=9999, linewidth=1.0)
#     ax.plot(xxx,yyy,'r', linewidth=3.0) #zorder=9999,

# These transients are only for the zoom in
# transient_colors = OrderedDict([
#     ("Ia", "red"),
#     ("Ia-91bg", "darkblue"),
#     ("II", "mediumspringgreen"),
#     ("Ia-91T", "mediumslateblue"),
#     ("SLSN-I", "mediumturquoise"),
#     ("IIb", "orange"),
#     ("Ic", "grey"),
#     ("Ib", "orchid"),
#     ("IIn", "dodgerblue")]
# )

transient_types = [
    "Ia",
    "Ia-91bg",
    "II",
    "SLSN-I",
    "IIb",
    "Ic",
    "Ib",
    "IIn",
    "Ia-91T"
]

transient_colors_id =OrderedDict([(val, i) for (i, val) in enumerate(transient_types)])
cmap = plt.get_cmap('Spectral_r')
transient_color_map = cmap(np.linspace(0, 1, len(transient_types)))

transients = []
with open('./web/models/linear/dm7_i_all_transients.txt') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=' ')
    next(csv_reader)
    for row in csv_reader:
        transient_type = row[0]
        dm7 = float(row[1])
        abs_peak = float(row[2])

        clr = transient_color_map[transient_colors_id[transient_type]]

        transients.append((abs_peak, dm7/7.0, clr, transient_type))

transient_legend = {
    "Ia":False,
    "Ia-91bg":False,
    "II":False,
    "SLSN-I":False,
    "IIb":False,
    "Ic":False,
    "Ib":False,
    "IIn":False,
    "Ia-91T":False,
}

mkr_sz = 25
_mew=1.0
for t in transients:
    if np.abs((t[1] - 0.003)) > 0.001: # cut out transients close the left edge of the plot
        if not transient_legend[t[3]]:
            transient_legend[t[3]] = True
            # ax.plot(t[1], t[0], "*", markerfacecolor=t[2], markeredgecolor=t[2], markersize=20, alpha=1.0, label=t[3])
            ax.plot(t[1], t[0], "*", markerfacecolor=t[2], markeredgecolor='k', markersize=mkr_sz, alpha=1.0, label=t[3],
                    mew=_mew)
        else:
            # ax.plot(t[1], t[0], "*", markerfacecolor=t[2], markeredgecolor=t[2], markersize=20, alpha=1.0)
            ax.plot(t[1], t[0], "*", markerfacecolor=t[2], markeredgecolor='k', markersize=mkr_sz, alpha=1.0, mew=_mew)

sm = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.inferno)
sm.set_array([]) # can be an empty list
# tks = np.linspace(min_prob, max_prob, 5)
tks = np.asarray([0.0, 0.15, 0.25, 0.35, max_prob])

cb = fig.colorbar(sm, ax=ax, orientation='vertical', fraction=0.05, pad=0.02, alpha=0.80) # fraction=0.04875

# cb.set_ticks([0.0, 0.944] + list(tks))
# cb.set_ticks(tks)
# tks_strings = ["0%", "25%", "50%", "70%", "95%"]
# tks_strings = ["0%", "15%", "25%", "37%"]
tks_strings = ["0%", "15%", "25%", "35%"]
cb.ax.set_yticklabels(tks_strings, fontsize=24)
cb.set_label("", fontsize=16, labelpad=9.0)
cb.ax.locator_params(nbins=5)
cb.ax.tick_params(length=8.0, width=2.0)

# WRONG - SSS17a
# dm7 = 0.56
# Mabs = -16.56

# ADJUSTED FROM KASEN
# SDSS_r dm1 = 0.690769300139
# Mabs = -16.1218

# ax.text(0.55, -15.8, "SSS17a", fontsize=22, ha="center", zorder=9999, color="yellow",
#             path_effects = [path_effects.withStroke(linewidth=2.5, foreground='black')])
ax.text(0.135, -16.06, "AT 2017gfo", fontsize=20, zorder=9999, color="white" #ha="center",
            , path_effects = [path_effects.withStroke(linewidth=10.0, foreground='black')])
# )
# ax.plot(0.690769300139, -16.1218, '*', color='black', markeredgecolor='black',
ax.plot(0.690769300139,     -16.1218, '*', color='black', markeredgecolor='black',
        markersize=mkr_sz, mew=_mew) #markeredgecolor="black" , alpha=0.25
print("Prob of Kasen model for SSS17a, r-band: %s" % _2d_func(0.690769300139, -16.1218))



# ax.annotate('0.9', xy=(0.4, -18.1), fontsize=24, color='red')
# ax.annotate('0.7', xy=(0.4, -17.4), fontsize=24, color='red')
# ax.annotate('0.5', xy=(0.4, -16.8), fontsize=24, color='red')
# ax.annotate('0.3', xy=(0.4, -16.4), fontsize=24, color='red')
# ax.annotate('0.1', xy=(0.4, -15.9), fontsize=24, color='red')
# ax.annotate('.9', xy=(-0.4, -13.0), fontsize=24, color='red')
# ax.annotate('.7', xy=(-0.16, -13.8), fontsize=24, color='red')
# ax.annotate('.5', xy=(-0.02, -14.4), fontsize=24, color='red')
# ax.annotate('.3', xy=(0.08, -14.8), fontsize=24, color='red')
# ax.annotate('.1', xy=(0.18, -15.1), fontsize=24, color='red')

ax.tick_params(axis='both', which='major', labelsize=24)

# These limits and legend are only for the zoom in
# ax.set_ylim([-15, -20])
ax.set_ylim([-14, -20])
# ax.set_xlim([3e-3, 9e-1])
ax.set_xlim([3e-3, 1.5])
ax.set_xscale('log')
leg = ax.legend(loc="upper right", fontsize=18, borderpad=0.35, handletextpad=0.0, labelspacing=0.4, framealpha=1.0,
                edgecolor='k', ncol=2, columnspacing=0.0)
leg.get_frame().set_linewidth(2.0)

plt.xlabel(r'$\Delta$M (mag $\mathrm{day^{-1}}$)',fontsize=32)
plt.ylabel(r'$\mathrm{M_{0}}$ (mag)',fontsize=32)

# Only for zoom-in
# plt.ylabel(r'$\mathrm{M_{peak}}$ [mag]',fontsize=32)

for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(2.0)
cb.outline.set_linewidth(2.0)

ax.tick_params(axis='both', which='major', labelsize=24, length=12.0, width=2)
ax.tick_params(axis='both', which='minor', labelsize=24, length=8.0, width=2)

fig.savefig('./web/events/S190425z/model_detection/linear/0425_Linear_Prob2Detect.png', bbox_inches='tight')
# fig.savefig('190814_Linear_Prob2Detect_all_Zoom.png', bbox_inches='tight')
plt.close('all')



# For the dM of 0.69, when do we get to 50%
# For the M of -16.12, when do we get to 50%
fixed_M = np.full_like(model_X_input, -16.1218) # fixed M in the shape of dM
fixed_dM = np.full_like(model_Y_input, 0.690769300139) # fixed dM in the shape of M

prob_fixed_M = _2d_func(model_X_input, fixed_M)
prob_fixed_dM = _2d_func(fixed_dM, model_Y_input)


# given a column (i.e. an "n"), find E for 50%
idx_M = find_nearest_index(prob_fixed_dM[:, 0], 0.5)
M_50 = model_Y_input[idx_M]
print("dM = 0.69, M_0 = %0.3E, prob ~ 0.5" % (M_50))

# given a row (i.e. an "E"), find n for 50%
idx_dM = find_nearest_index(prob_fixed_M[0, :], 0.5)
print("M_0 = -16.12, dM = %0.3E, prob ~ 0.5" % (model_X_input[idx_dM]))



# given a column (i.e. an "n"), find E for 50%
# idx_M= find_nearest_index(prob_fixed_dM[:, 0], 0.5)
# E_50 = model__input[idx_E]
# print("n = 3e-1, E = %0.3E, prob ~ 0.5" % (E_50))
#
# # given a row (i.e. an "E"), find n for 50%
# idx_n = find_nearest_index(prob_fixed_M[0, :], 0.5)
# print("E = 2.5, n = %0.3E, prob ~ 0.5" % (model_n_input[idx_n]))

# For the dM of 0.690769300139, when do we get to 50%
# For the M of -16.1218, when do we get to 50%

# model_X_input = np.linspace(-2.0, 10.0, 100)
# model_Y_input = np.linspace(-23.0, -11.0, 100)

# fixed_M = np.full_like((1,100), -16.1218) # fixed M in the shape of dM
# fixed_dM = np.full_like((1,100), 0.690769300139) # fixed dM in the shape of M

# prob_fixed_M = _2d_func(model_Y_input, fixed_M)
# prob_fixed_dM = _2d_func(fixed_dM, model_X_input)
#
# fig = plt.figure(figsize=(20, 10), dpi=1000)
#
# ax1 = fig.add_subplot(121)
# ax1.plot(model_Y_input, prob_fixed_dM[:, 0])
# ax1.set_xlabel('M', fontsize=32)
# ax1.set_ylabel('Prob', fontsize=32)
#
# ax2 = fig.add_subplot(122)
# ax2.plot(model_X_input, prob_fixed_M[0, :])
# ax2.set_xlabel('dM', fontsize=32)
# ax2.set_ylabel('Prob', fontsize=32)
#
# fig.savefig('190814_Linear_test.png', bbox_inches='tight')
#
#
# # along col
# idx_M = find_nearest_index(prob_fixed_dM[:, 0], 0.5)
# print("dM = 0.691, M = %0.2f, prob ~ 0.5" % (model_Y_input[idx_M]))
#
# # along row
# idx_dM =

# print("dM = %0.2f, M = -16.121, prob ~ 0.5" % (model_X_input[idx_dM]))
#
# # fixed_M = np.full_like(model_X_input, -16.1218) # fixed M in the shape of dM
# # fixed_dM = np.full_like(model_Y_input, 0.690769300139) # fixed dM in the shape of M
# # prob_fixed_M = _2d_func(model_X_input, fixed_M)
# # prob_fixed_dM = _2d_func(fixed_dM, model_Y_input)
#
# print("... Done.")
