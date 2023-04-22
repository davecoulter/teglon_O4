import time
import csv
import copy
import healpy as hp
import numpy as np
from astropy.cosmology import z_at_value, FlatLambdaCDM
import os

# from src.utilities.Database_Helpers import *
from web.src.utilities.Database_Helpers import *

glade_catalog_in = "./web/src/utilities/galaxy_catalog_files/GLADE_2.4.dat"
glade_catalog_mid = "./web/src/utilities/galaxy_catalog_files/GLADE_2.4.mid.dat"
glade_catalog_out = "./web/src/utilities/galaxy_catalog_files/GLADE_2.4.out.dat"
# reset output
glade_H0 = 70.0
glade_Om0 = 0.27
z_err = 300./3e5
glade_cosmo = FlatLambdaCDM(H0=glade_H0, Om0=glade_Om0)
def get_z_dist_and_err(z):
    glade_d = glade_cosmo.luminosity_distance(z)
    glade_d_symmetric_err = copy.deepcopy(glade_d) # glade_cosmo.luminosity_distance(z)

    idx = np.where(glade_d.value > 5.0)[0]
    if len(idx) > 0:
        z_high = z[idx] + z_err
        z_low = z[idx] - z_err
        glade_d_close = glade_cosmo.luminosity_distance(z_low)
        glade_d_far = glade_cosmo.luminosity_distance(z_high)
        _glade_d_symmetric_err = (glade_d_far - glade_d_close) / 2.0
        glade_d_symmetric_err[idx] = _glade_d_symmetric_err

    return glade_d.value, glade_d_symmetric_err.value

def get_B_luminosity_proxy(b_mag, z_dist):
    proxy = (z_dist**2)*(10.0**(-0.4*b_mag))
    return proxy

transform_GLADE = True
if transform_GLADE:
    print("Processing GLADE catalog...")
    if os.path.exists(glade_catalog_mid):
        with open(glade_catalog_mid, 'w') as reset_output:
            reset_output.write("")
    if os.path.exists(glade_catalog_out):
        with open(glade_catalog_out, 'w') as reset_output:
            reset_output.write("")

    t1 = time.time()
    row_count = 0
    with open(glade_catalog_in, 'r') as fileObject:
        csv_reader = csv.reader(fileObject, delimiter=' ', skipinitialspace=True)
        row_count = sum(1 for row in csv_reader)

    glade_z_to_transform = []
    glade_B_to_transform = []
    glade_ra_to_transform = []
    glade_dec_to_transform = []
    with open(glade_catalog_in, 'r') as csvfile_in:
        with open(glade_catalog_mid, 'a') as csvfile_out:
            # Read CSV lines
            csvreader = csv.reader(csvfile_in, delimiter=' ', skipinitialspace=True)
            csvwriter = csv.writer(csvfile_out, delimiter=',') # , escapechar=' ', quoting=csv.QUOTE_NONE

            for i, row in enumerate(csvreader):
                if i % 1e5 == 0:
                    frac = float(i)/float(row_count)*100
                    print("\tpoint transform: %0.0f%% complete... " % frac)

                # Only import galaxies with "good" distances and B magnitudes...
                z_str = row[10]
                app_B_mag_str = row[11]
                if row[5] == 'G' and z_str != 'null' and app_B_mag_str != 'null':
                    ra = float(row[6])
                    dec = float(row[7])
                    coord_str = "POINT(%s %s)" % (dec, ra - 180.0)  # Dec, RA order due to MySQL convention for lat/lon

                    glade_B_to_transform.append(float(app_B_mag_str))
                    glade_z_to_transform.append(float(z_str))
                    glade_ra_to_transform.append(np.radians(ra))
                    glade_dec_to_transform.append(np.radians(dec))

                    output_row = row
                    output_row += [coord_str]
                    csvwriter.writerow(output_row)

    transformed_z_dist, transformed_z_dist_err = get_z_dist_and_err(np.asarray(glade_z_to_transform))
    transformed_B_lum_proxy = get_B_luminosity_proxy(np.asarray(glade_B_to_transform), transformed_z_dist)
    B_lum_proxy_norm = np.sum(transformed_B_lum_proxy)
    b_weights = transformed_B_lum_proxy/B_lum_proxy_norm
    nside_2_pixl_indices = hp.ang2pix(nside=2, theta=0.5 * np.pi - np.asarray(glade_dec_to_transform),
                                      phi=np.asarray(glade_ra_to_transform))
    nside_4_pixl_indices = hp.ang2pix(nside=4, theta=0.5 * np.pi - np.asarray(glade_dec_to_transform),
                                      phi=np.asarray(glade_ra_to_transform))
    nside_8_pixl_indices = hp.ang2pix(nside=8, theta=0.5 * np.pi - np.asarray(glade_dec_to_transform),
                                      phi=np.asarray(glade_ra_to_transform))
    nside_16_pixl_indices = hp.ang2pix(nside=16, theta=0.5 * np.pi - np.asarray(glade_dec_to_transform),
                                      phi=np.asarray(glade_ra_to_transform))
    nside_32_pixl_indices = hp.ang2pix(nside=32, theta=0.5 * np.pi - np.asarray(glade_dec_to_transform),
                                      phi=np.asarray(glade_ra_to_transform))
    nside_64_pixl_indices = hp.ang2pix(nside=64, theta=0.5 * np.pi - np.asarray(glade_dec_to_transform),
                                      phi=np.asarray(glade_ra_to_transform))
    nside_128_pixl_indices = hp.ang2pix(nside=128, theta=0.5 * np.pi - np.asarray(glade_dec_to_transform),
                                      phi=np.asarray(glade_ra_to_transform))
    nside_256_pixl_indices = hp.ang2pix(nside=256, theta=0.5 * np.pi - np.asarray(glade_dec_to_transform),
                                      phi=np.asarray(glade_ra_to_transform))
    nside_512_pixl_indices = hp.ang2pix(nside=512, theta=0.5 * np.pi - np.asarray(glade_dec_to_transform),
                                      phi=np.asarray(glade_ra_to_transform))
    nside_1024_pixl_indices = hp.ang2pix(nside=1024, theta=0.5 * np.pi - np.asarray(glade_dec_to_transform),
                                      phi=np.asarray(glade_ra_to_transform))
    nside_2048_pixl_indices = hp.ang2pix(nside=2048, theta=0.5 * np.pi - np.asarray(glade_dec_to_transform),
                                         phi=np.asarray(glade_ra_to_transform))

    print("\n")
    with open(glade_catalog_mid, 'r') as csvfile_in:
        with open(glade_catalog_out, 'a') as csvfile_out:
            # Read CSV lines
            csvreader = csv.reader(csvfile_in, delimiter=',', skipinitialspace=True)
            csvwriter = csv.writer(csvfile_out, delimiter=',') # , escapechar=' ', quoting=csv.QUOTE_NONE

            for i, row in enumerate(csvreader):
                if i % 1e5 == 0:
                    frac = float(i)/float(row_count)*100
                    print("\tz transform: %0.0f%% complete... " % frac)

                z_dist = transformed_z_dist[i]
                z_dist_err = transformed_z_dist_err[i]
                nside_2_idx = nside_2_pixl_indices[i]
                nside_4_idx = nside_4_pixl_indices[i]
                nside_8_idx = nside_8_pixl_indices[i]
                nside_16_idx = nside_16_pixl_indices[i]
                nside_32_idx = nside_32_pixl_indices[i]
                nside_64_idx = nside_64_pixl_indices[i]
                nside_128_idx = nside_128_pixl_indices[i]
                nside_256_idx = nside_256_pixl_indices[i]
                nside_512_idx = nside_512_pixl_indices[i]
                nside_1024_idx = nside_1024_pixl_indices[i]
                nside_2048_idx = nside_2048_pixl_indices[i]
                b_lum_proxy = transformed_B_lum_proxy[i]
                b_weight = b_weights[i]
                output_row = row
                output_row += [z_dist, z_dist_err, nside_2_idx, nside_4_idx, nside_8_idx, nside_16_idx, nside_32_idx,
                               nside_64_idx, nside_128_idx, nside_256_idx, nside_512_idx, nside_1024_idx,
                               nside_2048_idx, b_lum_proxy, b_weight]
                csvwriter.writerow(output_row)
    t2 = time.time()
    print("\n********* start DEBUG ***********")
    print("GLADE extension execution time: %s" % (t2 - t1))

do_data_upload = True
if do_data_upload:
    print("Bulk uploading GLADE catalog...")
    t1 = time.time()
    upload_sql = """LOAD DATA LOCAL INFILE '%s'
            INTO TABLE Galaxy
            FIELDS TERMINATED BY ','
            LINES TERMINATED BY '\n'
            (PGC, Name_GWGC, Name_HyperLEDA, Name_2MASS, Name_SDSS_DR12, flag1, RA, _Dec, dist, dist_err, z,
            B, B_err, B_abs, J, J_err, H, H_err, K, K_err, flag2, flag3, @point, z_dist, z_dist_err, N2_Pixel_Index, 
            N4_Pixel_Index, N8_Pixel_Index, N16_Pixel_Index, N32_Pixel_Index, N64_Pixel_Index, N128_Pixel_Index,
            N256_Pixel_Index, N512_Pixel_Index, N1024_Pixel_Index, N2048_Pixel_Index, LuminosityProxy, Bweight)
            SET Coord := ST_PointFromText(@point, 4326);"""

    success = bulk_upload(upload_sql % glade_catalog_out)
    if not success:
        print("\nUnsuccessful bulk upload. Exiting...")
    else:
        print("Successful upload!")

    t2 = time.time()
    print("\n********* start DEBUG ***********")
    print("CSV upload execution time: %s" % (t2 - t1))