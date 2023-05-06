import os
import pickle
import time

import healpy as hp
import numpy as np
from astropy import units as u
import astropy.coordinates as coord
from dustmaps.config import config as dustmaps_config
from dustmaps.sfd import SFDQuery

from web.src.objects.Pixel_Element import *

utilities_base_dir = "./web/src/utilities"
dustmaps_config["data_dir"] = utilities_base_dir

nside128 = 128

t1 = time.time()

print("Generate all sky coords")
theta, phi = hp.pix2ang(nside=nside128, ipix=np.arange(hp.nside2npix(nside128)))
pix_coord = coord.SkyCoord(ra=np.rad2deg(phi), dec=np.rad2deg(0.5 * np.pi - theta), unit=(u.deg, u.deg))

print("Retrieving dust info")
sfd = SFDQuery()
ebv = sfd(pix_coord)

print("E(B-V) data length: %s" % len(ebv))
min_ebv = np.min(ebv)
max_ebv = np.max(ebv)
print("Max E(B-V): %s" % min_ebv)
print("Min E(B-V): %s" % max_ebv)

t2 = time.time()

print("\n********* start DEBUG ***********")
print("EBV query - execution time: %s" % (t2 - t1))
print("********* end DEBUG ***********\n")

