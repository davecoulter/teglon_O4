# region Imports
import os
import sys
import pickle
from scipy.interpolate import interp1d

from collections import OrderedDict
from astropy.table import Table
from astropy.time import Time
import pytz

from web.src.objects.Detector import *
from web.src.objects.Tile import *
from web.src.objects.Pixel_Element import *
from web.src.utilities.Database_Helpers import *

from astropy.cosmology import LambdaCDM
from astropy.coordinates import Distance
from astropy.coordinates.angles import Angle
from astropy.cosmology import z_at_value
from configparser import RawConfigParser

from scipy.special import erf
from multiprocessing import Pool
# endregion


configFile = './Settings.ini'
config = RawConfigParser()
config.read(configFile)

H0 = config.get('cosmology', 'H0')
O_M = config.get('cosmology', 'OMEGA_M')
O_DE = config.get('cosmology', 'OMEGA_DE')

def initial_z(pixel_data):
    cosmo = LambdaCDM(H0=H0, Om0=O_M, Ode0=O_DE)

    pix_index = pixel_data[0]
    mean_dist = pixel_data[1]

    d = Distance(mean_dist, u.Mpc)
    z = d.compute_z(cosmology=cosmo)

    return (pix_index, float(z))


def get_healpix_pixel_id(galaxy_info):
    phi, theta = np.radians(float(galaxy_info[8])), 0.5 * np.pi - np.radians(float(galaxy_info[9]))

    # map NSIDE is last argument of galaxy_info
    # return the galaxy_id with the pixel index in the NSIDE of the healpix map
    return (galaxy_info[0], hp.ang2pix(int(galaxy_info[-1]), theta, phi))

"""
pixel_data = (
    pix_key, 
    mean_dist, 
    dist_sigma, 
    mwe, 
    prob_2D, 
    band, 
    mjd, 
    pix_synopsis.lim_mags[band][mjd],
    model_func_dict,
    delta_mjd
)
"""
def integrate_pixel(pixel_data):
    band_mapping = {
        "sdss_u": "SDSS u",
        "sdss_g": "SDSS g",
        "sdss_r": "SDSS r",
        "sdss_i": "SDSS i",
        "sdss_z": "SDSS z",
        "johnson_B": "Landolt B",
        "johnson_V": "Landolt V",
        "johnson_R": "Landolt R",
        "johnson_I": "Landolt I",
        "ukirt_J": "UKIRT J",
        "ukirt_H": "UKIRT H",
        "ukirt_K": "UKIRT K",
        "Clear": "Clear"
    }
    reverse_band_mapping = {v: k for k, v in band_mapping.items()}

    pix_key = pixel_data[0]
    pix_index = pix_key[0]
    model_key = pix_key[1]

    mean_dist = pixel_data[1]
    dist_sigma = pixel_data[2]
    mwe = pixel_data[3]
    prob_2D = pixel_data[4]
    band = pixel_data[5]
    mjd = pixel_data[6]
    lim_mag = pixel_data[7]
    model_func_dict = pixel_data[8]
    delta_mjd = pixel_data[9]


    f_band = model_func_dict[reverse_band_mapping[band]]
    abs_in_band = f_band(delta_mjd)

    # compute distance upper bound, given the limiting magnitude:
    pwr = (lim_mag - abs_in_band + 5.0 - mwe) / 5.0 - 6.0
    d_upper = 10 ** pwr
    prob_to_detect = 0.0

    if d_upper > 0.0:
        # cdf = lambda d: norm.cdf(d, mean_dist, dist_sigma)
        # dist_norm = cdf(np.inf) - cdf(0.0) # truncating integral at d = 0 Mpc
        # prob_to_detect = (prob_2D/dist_norm) * (cdf(d_upper) - cdf(0.0))

        prob_to_detect = prob_2D * (0.5 * erf((d_upper - mean_dist) / (np.sqrt(2) * dist_sigma)) - \
                                    0.5 * erf(-mean_dist / (np.sqrt(2) * dist_sigma)))

    new_pix_key = (pix_index, model_key, band, mjd)
    return (new_pix_key, prob_to_detect)


class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--gw_id', default='', type="str",
                          help='LIGO superevent name, e.g. `S190425z`.')

        parser.add_option('--healpix_dir', default='./web/events/{GWID}', type="str",
                          help='Working directory.')

        parser.add_option('--healpix_file', default='bayestar.fits.gz', type="str",
                          help='Healpix filename. Default `bayestar.fits.gz`.')

        parser.add_option('--model_output_dir', default="./web/events/{GWID}/model_detection", type="str",
                          help='Directory for where to output processed models.')

        parser.add_option('--num_cpu', default="5", type="int",
                          help='Number of CPUs to use for multiprocessing')

        parser.add_option('--model_base_path', default="./web/models/{MODEL_TYPE}", type="str",
                          help='Path to models to process')

        parser.add_option('--model_type', default='kne', type="str",
                          help="Type of model to process. Values: linear, grb, kne")

        parser.add_option('--sub_dir', default="1", type="str",
                          help='Sub directory (for batching)')

        parser.add_option('--batch_dir', default="", type="str",
                          help='GRB Model sub sub directory (for batching... lol)')

        parser.add_option('--merger_time_MJD', default="58598.346134259256", type="float",
                          help='''Time of the merger in MJD. This is used to compute where on the light curve we are. 
                                  Default to GW190425''')

        parser.add_option('--clobber', action="store_true", default=False,
                            help='''Remove per event pickles and start with fresh database queries''')

        return (parser)

    def main(self):

        # region Sanity checks
        is_error = False

        formatted_healpix_dir = self.options.healpix_dir
        formatted_model_base_input_dir = self.options.model_base_path
        formatted_model_output_dir = self.options.model_output_dir
        if self.options.gw_id == "":
            is_error = True
            print("GWID is required.")
        else:
            if "{GWID}" in formatted_healpix_dir:
                formatted_healpix_dir = formatted_healpix_dir.replace("{GWID}", self.options.gw_id)

            if "{GWID}" in formatted_model_output_dir:
                formatted_model_output_dir = formatted_model_output_dir.replace("{GWID}", self.options.gw_id)

        if "{MODEL_TYPE}" in formatted_model_base_input_dir:
            formatted_model_base_input_dir = formatted_model_base_input_dir.replace("{MODEL_TYPE}", self.options.model_type)

        healpix_map_select = "SELECT id, NSIDE FROM HealpixMap WHERE GWID = '%s' and Filename = '%s'"
        map_check = query_db([healpix_map_select % (self.options.gw_id, self.options.healpix_file)])[0]
        healpix_map_id = None
        healpix_map_nside = None
        if len(map_check) != 1:
            is_error = True
            print("Map does not exist for combination of `%s` and `%s`!" %
                  (self.options.gw_id, self.options.healpix_file))
        else:
            healpix_map_id = int(map_check[0][0])
            healpix_map_nside = int(map_check[0][1])


        if is_error:
            print("Exiting...")
            return 1
        # endregion

        # region Setup Pickles
        gw_event_pickle_dir = "%s/%s" % (formatted_healpix_dir, "pickles")

        map_pickle_file = "%s/%s" % (gw_event_pickle_dir, "map_pix.pkl")
        resolved_z_pickle_file = "%s/%s" % (gw_event_pickle_dir, "resolved_z.pkl")
        detectors_pickle_file = "%s/%s" % (gw_event_pickle_dir, "detectors.pkl")
        obs_tiles_pickle_file = "%s/%s" % (gw_event_pickle_dir, "obs_tiles.pkl")
        enclosed_pix_pickle_file = "%s/%s" % (gw_event_pickle_dir, "enc_pix.pkl")
        if self.options.clobber:
            if os.path.exists(map_pickle_file):
                os.remove(map_pickle_file)
            if os.path.exists(resolved_z_pickle_file):
                os.remove(resolved_z_pickle_file)
            if os.path.exists(detectors_pickle_file):
                os.remove(detectors_pickle_file)
            if os.path.exists(obs_tiles_pickle_file):
                os.remove(obs_tiles_pickle_file)
            if os.path.exists(enclosed_pix_pickle_file):
                os.remove(enclosed_pix_pickle_file)

        utilities_base_dir = "/app/web/src/utilities"
        global_pickle_dir = "%s/pickles/" % utilities_base_dir

        # LOADING NSIDE 128 SKY PIXELS AND EBV INFORMATION
        print("\nLoading NSIDE 128 pixels...")
        nside128 = 128
        N128_dict = None
        with open('%s/N128_dict.pkl' % global_pickle_dir, 'rb') as handle:
            N128_dict = pickle.load(handle)
        del handle

        print("\nLoading existing EBV...")
        ebv = None
        with open('%s/ebv.pkl' % global_pickle_dir, 'rb') as handle:
            ebv = pickle.load(handle)
        # endregion

        prep_start = time.time()
        print("Processes to use: %s" % self.options.num_cpu)

        merger_time = self.options.merger_time_MJD
        t = Time(merger_time, format='mjd')
        print("** Running models for MJD: %s; UTC: %s **" % (merger_time, t.to_datetime(timezone=pytz.utc)))

        # region Load Models from Disk
        model_path = None
        if self.options.model_type == "grb":
            model_path = "%s/%s/%s" % (formatted_model_base_input_dir, self.options.batch_dir, self.options.sub_dir)
        else:
            model_path = "%s/%s" % (formatted_model_base_input_dir, self.options.sub_dir)

        model_files = []
        for file_index, file in enumerate(os.listdir(model_path)):
            if file.endswith(".dat"):
                model_files.append("%s/%s" % (model_path, file))

        # Band abbreviation, band_id mapping
        band_mapping = {
            "sdss_u": "SDSS u",
            "sdss_g": "SDSS g",
            "sdss_r": "SDSS r",
            "sdss_i": "SDSS i",
            "sdss_z": "SDSS z",
            "johnson_B": "Landolt B",
            "johnson_V": "Landolt V",
            "johnson_R": "Landolt R",
            "johnson_I": "Landolt I",
            "ukirt_J": "UKIRT J",
            "ukirt_H": "UKIRT H",
            "ukirt_K": "UKIRT K",
            "Clear": "Clear"
        }
        reverse_band_mapping = {v: k for k, v in band_mapping.items()}
        models = {}
        for index, mf in enumerate(model_files):
            model_table = Table.read(mf, format='ascii.ecsv')
            model_props = model_table.meta['comment']
            model_key = None

            if self.options.model_type == 'linear':
                M = float(model_props[0].split("=")[1])
                dM = float(model_props[1].split("=")[1])
                model_key = (M, dM)
            elif self.options.model_type == 'grb':
                E = float(model_props[1].split("=")[1])
                n = float(model_props[2].split("=")[1])
                theta_obs = float(model_props[3].split("=")[1])
                model_key = (E, n, theta_obs)
            elif self.options.model_type == 'kne':
                vej = float(model_props[0].split('=')[-1].strip())
                mej = float(model_props[1].split('=')[-1].strip())
                ye = float(model_props[2].split('=')[-1].strip())
                model_key = (vej, mej, ye)
            else:
                raise Exception("Unkown model type! `%s`" % self.options.model_type)

            base_name = os.path.basename(mf)
            print("Loading `%s`" % base_name)

            available_model_bands = []
            # Create MaskedTable from Table
            mt_masked = Table(model_table, masked=True)

            # Mask out NaN. Generate composite mask from union of all column masks
            for col in mt_masked.columns:
                if col not in available_model_bands:
                    available_model_bands.append(col)

                nan_indices = np.where(np.isnan(mt_masked[col]))
                mt_masked.mask[nan_indices] = True

            # Temporal mask: mask times > 15 days
            late_time_indices = np.where(mt_masked['time'] > 15.0)
            mt_masked.mask[late_time_indices] = True

            # mask will be set for all cols; just take this one
            composite_mask = mt_masked['time'].mask
            model_time = np.asarray(mt_masked['time'][~composite_mask])

            # Spin up interpolation funcs
            models[model_key] = {}
            for col in mt_masked.columns:
                y_vals = np.asarray(mt_masked[col][~composite_mask])
                light_curve_func = interp1d(model_time, y_vals, fill_value="extrapolate")
                models[model_key][col] = light_curve_func
        # endregion

        band_dict_by_name = {
            "SDSS u": (1, "SDSS u", 4.239),
            "SDSS g": (2, "SDSS g", 3.303),
            "SDSS r": (3, "SDSS r", 2.285),
            "SDSS i": (4, "SDSS i", 1.698),
            "SDSS z": (5, "SDSS z", 1.263),
            "Landolt B": (6, "Landolt B", 3.626),
            "Landolt V": (7, "Landolt V", 2.742),
            "Landolt R": (8, "Landolt R", 2.169),
            "Landolt I": (9, "Landolt I", 1.505),
            "UKIRT J": (10, "UKIRT J", 0.709),
            "UKIRT H": (11, "UKIRT H", 0.449),
            "UKIRT K": (12, "UKIRT K", 0.302),
            "Clear": (13, "Clear", 0.91)
        }
        band_dict_by_id = {
            1: (1, "SDSS u", 4.239),
            2: (2, "SDSS g", 3.303),
            3: (3, "SDSS r", 2.285),
            4: (4, "SDSS i", 1.698),
            5: (5, "SDSS z", 1.263),
            6: (6, "Landolt B", 3.626),
            7: (7, "Landolt V", 2.742),
            8: (8, "Landolt R", 2.169),
            9: (9, "Landolt I", 1.505),
            10: (10, "UKIRT J", 0.709),
            11: (11, "UKIRT H", 0.449),
            12: (12, "UKIRT K", 0.302),
            13: (13, "Clear", 0.91)
        }

        map_nside = 256  # 0425
        map_pixels = None
        if not os.path.exists(map_pickle_file):
            map_pix_select = '''
                SELECT 
                    hp.id,
                    hp.Pixel_Index,
                    hp.Mean,
                    hp.Stddev,
                    hp.Prob as '2D_Prob',
                    hpc.NetPixelProb as '4D_Prob',
                    sp.Pixel_Index as N128_Pixel_Index
                FROM 
                    HealpixPixel hp
                JOIN 
                    HealpixPixel_Completeness hpc on hpc.HealpixPixel_id = hp.id
                JOIN 
                    SkyPixel sp on sp.id = hp.N128_SkyPixel_id
                JOIN 
                    ObservedTile_HealpixPixel ot_hp on ot_hp.HealpixPixel_id = hp.id
                JOIN 
                    ObservedTile ot on ot.id = ot_hp.ObservedTile_id
                WHERE
                    hp.HealpixMap_id = %s AND
                    sp.NSIDE = 128;
            '''
            map_pixels = query_db([map_pix_select % healpix_map_id])[0]
            with open(map_pickle_file, 'wb') as handle:
                pickle.dump(map_pixels, handle, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open(map_pickle_file, 'rb') as handle:
                map_pixels = pickle.load(handle)

        print("Retrieved %s map pixels..." % len(map_pixels))

        # Initialize map pix dict for later access
        map_pixel_dict = OrderedDict()

        class Pixel_Synopsis():
            def __init__(self, mean_dist, dist_sigma, prob_2D, pixel_index, N128_index, pix_ebv):  # z
                self.mean_dist = mean_dist
                self.dist_sigma = dist_sigma
                self.forced_norm = 0.0

                self.prob_2D = prob_2D
                self.pixel_index = pixel_index
                self.N128_index = N128_index
                self.pix_ebv = pix_ebv
                self.z = 0.0

                # From the tiles that contain this pixel
                # band:value
                self.measured_bands = []
                self.lim_mags = OrderedDict()
                self.delta_mjds = OrderedDict()

                # From the model (only select the bands that have been imaged)
                self.A_lambda = OrderedDict()  # band:value

                # Final calculation
                # model:band:value
                self.best_integrated_probs = OrderedDict()

            def __str__(self):
                return str(self.__dict__)

        count_bad_pixels = 0

        initial_integrands = []
        for p in map_pixels:
            mean_dist = float(p[2])
            dist_sigma = float(p[3])
            prob_2D = float(p[5]) # this is prob 4D for now
            pixel_index = int(p[1])
            N128_pixel_index = int(p[6])
            pix_ebv = ebv[N128_pixel_index]

            if mean_dist == 0.0:
                # distance did not converge for this pixel. pass...
                print("Bad Index: %s" % pixel_index)
                count_bad_pixels += 1
                continue

            p_new = Pixel_Synopsis(
                mean_dist,
                dist_sigma,
                prob_2D,
                pixel_index,
                N128_pixel_index,
                pix_ebv)

            initial_integrands.append((p_new.pixel_index, p_new.mean_dist, p_new.dist_sigma))
            map_pixel_dict[pixel_index] = p_new

        resolved_z = None
        if not os.path.exists(resolved_z_pickle_file):
            print("Starting z-cosmo pool (%s distances)..." % len(initial_integrands))
            it1 = time.time()

            # Calculate everything based on num_cpu...
            chunks = int(len(initial_integrands) // self.options.num_cpu) + 1
            print("\tNum CPUs: %s" % self.options.num_cpu)
            print("\t\tTasks per CPU: %s" % chunks)

            with Pool(processes=self.options.num_cpu) as pool:
                resolved_z = list(pool.imap_unordered(initial_z, initial_integrands, chunksize=chunks))

            it2 = time.time()
            print("... finished z-cosmo pool: %s [seconds]" % (it2 - it1))

            with open(resolved_z_pickle_file, 'wb') as handle:
                pickle.dump(resolved_z, handle, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open(resolved_z_pickle_file, 'rb') as handle:
                resolved_z = pickle.load(handle)

        for rz in resolved_z:
            map_pixel_dict[rz[0]].z = rz[1]

        print("\nMap pixel dict complete. %s bad pixels." % count_bad_pixels)

        # Get serialized detectors and tiles
        detector_dictionary = None
        print("\nLoading Detectors...")
        if not os.path.exists(detectors_pickle_file):
            detector_select = '''
                            SELECT id, Name, ST_AsText(Poly) FROM Detector; 
                        '''
            detector_result = query_db([detector_select])[0]
            detector_dictionary = {int(dr[0]): Detector(dr[1], Detector.get_detector_vertices_from_teglon_db(dr[2]),
                                              detector_id=int(dr[0])) for dr in detector_result}

            with open(detectors_pickle_file, 'wb') as handle:
                pickle.dump(detector_dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open(detectors_pickle_file, 'rb') as handle:
                detector_dictionary = pickle.load(handle)

        print("\nLoading All Tiles...")
        ot_result = None
        if not os.path.exists(obs_tiles_pickle_file):
            observed_tile_select = '''
                SELECT 
                    id, 
                    Detector_id,
                    FieldName,
                    RA,
                    _Dec,
                    Band_id,
                    MJD,
                    Mag_Lim, 
                    PositionAngle
                FROM
                    ObservedTile
                WHERE 
                    HealpixMap_id = %s;                      
            '''
            ot_result = query_db([observed_tile_select % healpix_map_id])[0]

            with open(obs_tiles_pickle_file, 'wb') as handle:
                pickle.dump(ot_result, handle, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open(obs_tiles_pickle_file, 'rb') as handle:
                ot_result = pickle.load(handle)

        enc_pix_by_tile_id_dict = None
        if not os.path.exists(enclosed_pix_pickle_file):
            enc_pix_select = '''
                SELECT 
                    ot.id as Tile_id,
                    hp.Pixel_Index
                FROM
                    ObservedTile_HealpixPixel ot_hp
                JOIN
                    HealpixPixel hp on hp.id = ot_hp.HealpixPixel_id
                JOIN 
                    ObservedTile ot on ot.id = ot_hp.ObservedTile_id
                WHERE
                    ot.HealpixMap_id = %s;
            '''
            pix_result = query_db([enc_pix_select % healpix_map_id])[0]
            enc_pix_by_tile_id_dict = {}
            for pr in pix_result:
                tile_id = int(pr[0])
                pix_index = int(pr[1])
                if tile_id not in enc_pix_by_tile_id_dict:
                    enc_pix_by_tile_id_dict[tile_id] = []

                if pix_index not in enc_pix_by_tile_id_dict[tile_id]:
                    enc_pix_by_tile_id_dict[tile_id].append(pix_index)

            with open(enclosed_pix_pickle_file, 'wb') as handle:
                pickle.dump(enc_pix_by_tile_id_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open(enclosed_pix_pickle_file, 'rb') as handle:
                enc_pix_by_tile_id_dict = pickle.load(handle)

        observed_tiles = []
        for ot in ot_result:
            id = int(ot[0])
            ra = float(ot[3])
            dec = float(ot[4])
            pa = float(ot[8])
            detect_id = int(ot[1])
            fieldname = ot[2]
            mjd = float(ot[6])
            mag_lim = float(ot[7])
            band_id = int(ot[5])

            # Hack
            band_name = band_dict_by_id[band_id][1]
            if band_name not in reverse_band_mapping:
                continue
            model_band_name = reverse_band_mapping[band_name]
            if model_band_name not in available_model_bands:
                continue

            t = Tile(ra, dec, nside=map_nside, position_angle_deg=pa, tile_id=id,
                     detector=detector_dictionary[detect_id])

            t.field_name = fieldname
            t.mjd = mjd
            t.mag_lim = mag_lim
            t.band_id = band_id

            observed_tiles.append(t)

        print("Observed Tiles: %s" % len(observed_tiles))

        print("\nUpdating pixel `delta_mjds` and `lim_mags`...")
        # region Initialize models
        # For each tile:
        #   we want the MJD of observation, and add that to the list of a pixels' MJD collection.
        #   we want the limiting mag, add that to the list of a pixel's lim mag collection
        for t in observed_tiles:
            # pix_indices = t.enclosed_pixel_indices
            pix_indices = enc_pix_by_tile_id_dict[t.id]

            for i in pix_indices:
                # get band from id...
                band = band_dict_by_id[t.band_id]
                band_name = band[1]

                # Some pixels are omitted because their distance information did not converge
                if i not in map_pixel_dict:
                    continue
                pix_synopsis_new = map_pixel_dict[i]

                if band_name not in pix_synopsis_new.measured_bands:
                    pix_synopsis_new.measured_bands.append(band_name)
                    pix_synopsis_new.delta_mjds[band_name] = {}
                    pix_synopsis_new.lim_mags[band_name] = {}

                # DEBUG
                if t.mjd in pix_synopsis_new.lim_mags[band_name]:
                    if pix_synopsis_new.lim_mags[band_name][t.mjd] != t.mag_lim:
                        print("Collision! Tile: %s" % t)
                        print("MJD: %s" % t.mjd)
                        print("Previous mag lim: %s" % pix_synopsis_new.lim_mags[band_name][t.mjd])
                        print("New mag lim: %s" % t.mag_lim)
                # DEBUG

                # time-dilate the delta_mjd, which is used to get the Abs Mag LC point
                pix_synopsis_new.delta_mjds[band_name][t.mjd] = (t.mjd - merger_time) / (1.0 + pix_synopsis_new.z)
                pix_synopsis_new.lim_mags[band_name][t.mjd] = (t.mag_lim)

        print("\nInitializing %s models..." % len(models))
        # endregion

        print("\nUpdating pixel `A_lambda` and `model_observer_time_arr`...")
        # region Update pixel MWE and time-dilation
        # Set pixel `model_observer_time_arr` and `A_lambda`
        for pix_index, pix_synopsis in map_pixel_dict.items():
            for model_param_tuple, model_dict in models.items():
                for model_col in model_dict.keys():
                    if model_col in band_mapping:

                        band = band_dict_by_name[band_mapping[model_col]]
                        band_id = band[0]
                        band_name = band[1]
                        band_coeff = band[2]

                        if band_name in pix_synopsis.measured_bands:
                            if band_name not in pix_synopsis.A_lambda:
                                pix_a_lambda = pix_synopsis.pix_ebv * band_coeff
                                pix_synopsis.A_lambda[band_name] = pix_a_lambda

        prep_end = time.time()
        print("Prep time: %s [seconds]" % (prep_end - prep_start))
        # endregion

        # region Integrate
        compute_start = time.time()
        # NEW Do the calculation...
        print("\nUpdating `map_pixel_dict`...")
        count = 0

        for pix_index, pix_synopsis in map_pixel_dict.items():
            for band in pix_synopsis.measured_bands:
                for model_param_tuple, model_dict in models.items():
                    pixel_delta_mjd = pix_synopsis.delta_mjds[band]

                    for i, (mjd, delta_mjd) in enumerate(pixel_delta_mjd.items()):
                        if model_param_tuple not in pix_synopsis.best_integrated_probs:
                            pix_synopsis.best_integrated_probs[model_param_tuple] = {}

                        if band not in pix_synopsis.best_integrated_probs[model_param_tuple]:
                            pix_synopsis.best_integrated_probs[model_param_tuple][band] = {mjd: 0.0}
            count += 1
            if count % 1000 == 0:
                print("Processed: %s" % count)

        compute_end = time.time()
        print("Compute time: %s [seconds]" % (compute_end - compute_start))

        # Compose integrands
        print("Building integrands...")
        integrands_start = time.time()
        pixels_to_integrate = []
        for model_key, model_func_dict in models.items():
            for pix_index, pix_synopsis in map_pixel_dict.items():
                for band in pix_synopsis.measured_bands:
                    mwe = pix_synopsis.A_lambda[band]
                    prob_2D = pix_synopsis.prob_2D
                    mean_dist = pix_synopsis.mean_dist
                    dist_sigma = pix_synopsis.dist_sigma

                    for i, (mjd, delta_mjd) in enumerate(pix_synopsis.delta_mjds[band].items()):
                        pix_key = (pix_index, model_key)
                        pixels_to_integrate.append((
                            pix_key,
                            mean_dist,
                            dist_sigma,
                            mwe,
                            prob_2D,
                            band,
                            mjd,
                            pix_synopsis.lim_mags[band][mjd],
                            model_func_dict,
                            delta_mjd
                        ))
        integrands_end = time.time()
        print("Building integrands time: %s [seconds]" % (integrands_end - integrands_start))

        print("Integrating %s total model/pixel combinations..." % len(pixels_to_integrate))
        mp_start = time.time()

        integrated_pixels = None
        # Calculate everything based on num_cpu...

        chunks = int(len(pixels_to_integrate) // self.options.num_cpu // 3) + 1
        print("\tNum CPUs: %s" % self.options.num_cpu)
        print("\t\tTasks per CPU: %s" % chunks)

        with Pool(processes=self.options.num_cpu, maxtasksperchild=3) as pool:
            integrated_pixels = list(pool.imap_unordered(integrate_pixel, pixels_to_integrate, chunksize=chunks))

        for ip in integrated_pixels:
            pix_key = ip[0]
            pix_index = pix_key[0]
            model_param_tuple = pix_key[1]
            band = pix_key[2]
            mjd = pix_key[3]
            prob = ip[1]
            map_pixel_dict[pix_index].best_integrated_probs[model_param_tuple][band][mjd] = prob
        mp_end = time.time()

        print("Integration time: %s [seconds]" % (mp_end - mp_start))
        # endregion

        # probability by model by band by pixel, for all epochs per pixel
        probs_by_model_band_pix_index = {}

        # summed pixel probability by model and band
        probs_by_model_band = {}

        # final probability by model
        net_prob_by_model = {}

        for model_param_tuple, model_dict in models.items():

            # Initialize the model key for each dictionary
            if model_param_tuple not in net_prob_by_model:
                net_prob_by_model[model_param_tuple] = {}

            if model_param_tuple not in probs_by_model_band:
                probs_by_model_band[model_param_tuple] = {}

            if model_param_tuple not in probs_by_model_band_pix_index:
                probs_by_model_band_pix_index[model_param_tuple] = {}

            for pix_index, pix_synopsis in map_pixel_dict.items():
                for band in pix_synopsis.measured_bands:

                    # Initialize the band key for each dictionary
                    if band not in probs_by_model_band_pix_index[model_param_tuple]:
                        probs_by_model_band_pix_index[model_param_tuple][band] = {}

                    if band not in probs_by_model_band[model_param_tuple]:
                        probs_by_model_band[model_param_tuple][band] = {}

                    # sanity -- if the pixel prob is itself 0.0, we can just skip
                    if pix_synopsis.prob_2D <= 0:
                        probs_by_model_band_pix_index[model_param_tuple][band][pix_index] = 0.0
                        continue

                    # For each pixel:
                    #   take the compliment from each epoch's integrated probability
                    #   take the product of these compliments (the product of the individual marginal probabilities
                    #       of a non-detection)
                    #   take the compliment of this join marginal probability to find the prob of at least one detection

                    # test = 1.0 - reduce(lambda y, z: y * z,
                    #                  map(lambda x: (1.0 - x),
                    #                      [p[1] for p in pix_synopsis.best_integrated_probs[model_param_tuple][band].items()]))

                    product_of_compliments = 1.0

                    test_probs = []
                    for mjd, prob in pix_synopsis.best_integrated_probs[model_param_tuple][band].items():
                        # product_of_compliments *= (1.0 - prob)
                        # product_of_compliments *= (pix_synopsis.prob_2D - prob)
                        product_of_compliments *= (1.0 - prob / pix_synopsis.prob_2D)
                        test_probs.append(prob)

                    # compliment_of_product_of_compliments = 1.0 - product_of_compliments
                    # compliment_of_product_of_compliments = pix_synopsis.prob_2D - product_of_compliments
                    compliment_of_product_of_compliments = pix_synopsis.prob_2D * (1 - product_of_compliments)
                    # print(np.max(test_probs), compliment_of_product_of_compliments)

                    probs_by_model_band_pix_index[model_param_tuple][band][
                        pix_index] = compliment_of_product_of_compliments

            # Sum all pixels for a given band and model, that makes up the per model, per band prob
            for band, pix_prob_dict in probs_by_model_band_pix_index[model_param_tuple].items():

                # test2 = reduce(lambda x, y: x + y, [p[1] for p in pix_prob_dict.items()])

                summed_prob = 0.0
                for pix_index, prob in pix_prob_dict.items():
                    summed_prob += prob

                probs_by_model_band[model_param_tuple][band] = summed_prob

                # print("Summed prob for band [%s]: %s" % (band, summed_prob))

            # Finally, as above, take the compliment of the product of the compliments for each band
            # to compute the probability of at least 1 detection in any band
            # test3 = 1.0 - reduce(lambda y, z: y * z,
            #                  map(lambda x: (1.0 - x),
            #                      [p[1] for p in probs_by_model_band[model_param_tuple].items()]))

            total_covered_prob = 0.0
            for pix_index, pix_synopsis in map_pixel_dict.items():
                total_covered_prob += pix_synopsis.prob_2D

            # print("Total covered prob: %s" % total_covered_prob)

            net_product_of_compliments = 1.0
            for band, sum_prob in probs_by_model_band[model_param_tuple].items():
                # net_product_of_compliments *= (total_covered_prob - sum_prob)
                net_product_of_compliments *= (1.0 - sum_prob / total_covered_prob)

            # net_compliment_of_product_of_compliments = total_covered_prob - net_product_of_compliments
            net_compliment_of_product_of_compliments = total_covered_prob * (1.0 - net_product_of_compliments)

            # print("Total model prob: %s" % net_compliment_of_product_of_compliments)

            net_prob_by_model[model_param_tuple] = net_compliment_of_product_of_compliments

            # test4 = 1

        # Build ascii.ecsv formatted output
        cols = None
        dtype = None
        output_filename = None

        if self.options.model_type == 'linear':
            cols = ['M', 'dM', 'Prob']
            dtype = ['f8', 'f8', 'f8']
            output_filename = 'Detection_Results_Linear_%s.prob' % self.options.sub_dir
        elif self.options.model_type == 'grb':
            cols = ['E', 'n', 'theta_obs', 'Prob']
            dtype = ['f8', 'f8', 'f8', 'f8']
            output_filename = 'Detection_Results_SGRB_%s_%s.prob' % (self.options.batch_dir, self.options.sub_dir)
        elif self.options.model_type == 'kne':
            cols = ['vej', 'mej', 'ye', 'Prob']
            dtype = ['f8', 'f8', 'f8', 'f8']
            output_filename = 'Detection_KNe_%s.prob' % self.options.sub_dir

        result_table = Table(dtype=dtype, names=cols)

        if self.options.model_type == 'linear':
            for model_param_tuple, prob in net_prob_by_model.items():
                result_table.add_row([model_param_tuple[0], model_param_tuple[1], prob])
        elif self.options.model_type == 'grb' or self.options.model_type == 'kne':
            for model_param_tuple, prob in net_prob_by_model.items():

                try:
                    result_table.add_row([model_param_tuple[0], model_param_tuple[1], model_param_tuple[2], prob])
                except Exception as e:
                    print("Here")
                    test = 1


                    print(e)
                    print("\nExiting")
                    return 1

        result_table.write("%s/%s" % (formatted_model_output_dir, output_filename),
                           overwrite=True, format='ascii.ecsv')


if __name__ == "__main__":
    useagestring = """python ComputeModelDetection.py [options]


Assumes .dat files with Astropy ECSV format: 
https://docs.astropy.org/en/stable/api/astropy.io.ascii.Ecsv.html

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
    print("Teglon `ComputeModelDetection` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")
