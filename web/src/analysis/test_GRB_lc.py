import matplotlib
matplotlib.use("Agg")
from matplotlib.ticker import ScalarFormatter

import matplotlib as mpl
import matplotlib.pyplot as plt
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

model_1 = "./web/models/grb/onaxis/5/grb_0.0000_100.0000_0.0000100.dat"
model_2 = "./web/models/grb/onaxis/4/grb_0.0000_100.0000_0.0015849.dat"
model_3 = "./web/models/grb/onaxis/7/grb_0.0000_100.0000_3.9810717.dat"

model_4 = "./web/models/grb/onaxis/7/grb_0.0000_100.0000_0.0000063.dat"
model_5 = "./web/models/grb/onaxis/7/grb_0.0000_100.0000_0.0000158.dat"
model_6 = "./web/models/grb/onaxis/7/grb_0.0000_100.0000_0.1000000.dat"
model_7 = "./web/models/grb/onaxis/7/grb_0.0000_100.0000_0.0025119.dat"

model_table1 = Table.read(model_1, format='ascii.ecsv')
model_table2 = Table.read(model_2, format='ascii.ecsv')
model_table3 = Table.read(model_3, format='ascii.ecsv')

model_table4 = Table.read(model_4, format='ascii.ecsv')
model_table5 = Table.read(model_5, format='ascii.ecsv')
model_table6 = Table.read(model_6, format='ascii.ecsv')
model_table7 = Table.read(model_7, format='ascii.ecsv')

t = model_table1["time"]
sdss_r1 = model_table1["sdss_r"]
sdss_r2 = model_table2["sdss_r"]
sdss_r3 = model_table3["sdss_r"]

sdss_r4 = model_table4["sdss_r"]
sdss_r5 = model_table5["sdss_r"]
sdss_r6 = model_table6["sdss_r"]
sdss_r7 = model_table7["sdss_r"]

print("Plotting...")
fig = plt.figure(figsize=(10, 10), dpi=300)
ax = fig.add_subplot(111)
ax.plot(t, sdss_r1, linestyle="-", color="blue", linewidth=3.0, label="Ek,iso=10; n=1e-5")
ax.plot(t, sdss_r2, linestyle="-", color="red", linewidth=3.0, label="Ek,iso=10; n=1.5e-3")

ax.plot(t, sdss_r3, linestyle="-", color="magenta", linewidth=3.0, label="Ek,iso=10; n=3.98")
ax.plot(t, sdss_r4, linestyle="-", color="green", linewidth=3.0, label="Ek,iso=10; n=1e-6")

ax.plot(t, sdss_r5, linestyle="-", color="orange", linewidth=3.0, label="Ek,iso=10; n=1.5e-5")
ax.plot(t, sdss_r6, linestyle="-", color="black", linewidth=3.0, label="Ek,iso=10; n=0.1")

ax.plot(t, sdss_r7, linestyle="-", color="purple", linewidth=3.0, label="Ek,iso=10; n=2.5e-3")

# ax.set_ylim([-18 ,-10])
ax.set_xlim([0, 2])
ax.invert_yaxis()

ax.set_ylabel("Absolute Mag. (AB)", fontsize=32, labelpad=9.0)
ax.set_xlabel("Days From Merger", fontsize=32)
ax.legend()

fig.savefig('./web/events/S190425z/model_detection/grb/onaxis/grb_lc_test.png', bbox_inches='tight')



print("... Done.")