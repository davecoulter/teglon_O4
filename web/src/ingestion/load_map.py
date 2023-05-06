import os
import multiprocessing as mp
from scipy.special import erf
import healpy as hp
from ligo.skymap import distance

from web.src.objects.Detector import *
from web.src.objects.Tile import *
from web.src.objects.Pixel_Element import *
from web.src.utilities.Database_Helpers import *

import psutil
import urllib.request
from astropy.time import Time

from ligo.gracedb.rest import GraceDb

initialize_start = time.time()

def initialize_tile(tile):
    tile.enclosed_pixel_indices
    # tile.query_polygon_string
    return tile

class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--gw_id', default="", type="str",
                          help='LIGO superevent name, e.g. `S190425z` ')

        parser.add_option('--healpix_dir', default='./web/events/{GWID}', type="str",
                          help='Directory for where to look for the healpix file (Default = ../Events/{GWID})')

        parser.add_option('--healpix_file', default="bayestar.fits.gz", type="str", help='healpix filename')

        parser.add_option('--analysis_mode', action="store_true", default=False,
                          help='''Upload the healpix file at the native resolution (Default = False), with all 
                            pixels. Otherwise, upload a rescaled version, and only the 90th percentile pixels.''')

        parser.add_option('--clobber', action="store_true", default=False,
                          help='''Replace the existing (superevent, healpix_file) combination if it already exists. 
                          If `clobber` not specified and combination exists, return error message.''')

        parser.add_option('--skip_swope', action="store_true", default=False,
                          help='''Do not register Swope tiles for this map (Default = False)''')

        parser.add_option('--skip_thacher', action="store_true", default=False,
                          help='''Do not register Thacher tiles for this map (Default = False)''')

        parser.add_option('--skip_t80', action="store_true", default=False,
                          help='''Do not register T80S_T80S-Cam tiles for this map (Default = False)''')

        return (parser)

    def main(self):

        api_endpoint = "https://gracedb.ligo.org/api/"
        utilities_base_dir = "/app/web/src/utilities"
        pickle_output_dir = "%s/pickles/" % utilities_base_dir

        isDEBUG = False
        build_map = True
        build_pixels = True
        build_tile_pixel_relation = True
        build_galaxy_pixel_relation = True
        build_completeness_func = False
        build_galaxy_weights = True

        healpix_map_select = "SELECT id, NSIDE FROM HealpixMap WHERE GWID = '%s' and Filename = '%s'"

        is_error = False

        # Parameter checks
        if self.options.gw_id == "":
            is_error = True
            print("GWID is required.")

        formatted_healpix_dir = self.options.healpix_dir
        formatted_candidates_dir = self.options.healpix_dir + "/candidates"
        formatted_obs_tiles_dir = self.options.healpix_dir + "/observed_tiles"
        formatted_model_dir = self.options.healpix_dir + "/model_detection"
        if "{GWID}" in formatted_healpix_dir:
            formatted_healpix_dir = formatted_healpix_dir.replace("{GWID}", self.options.gw_id)
            formatted_candidates_dir = formatted_candidates_dir.replace("{GWID}", self.options.gw_id)
            formatted_obs_tiles_dir = formatted_obs_tiles_dir.replace("{GWID}", self.options.gw_id)
            formatted_model_dir = formatted_model_dir.replace("{GWID}", self.options.gw_id)

        hpx_path = "%s/%s" % (formatted_healpix_dir, self.options.healpix_file)

        if self.options.healpix_file == "":
            is_error = True
            print("You must specify which healpix file to process.")

        map_check = query_db([healpix_map_select % (self.options.gw_id, self.options.healpix_file)])[0]
        if len(map_check) > 0:

            current_map_id = map_check[0][0]

            if self.options.clobber:
                print("\n************ The combination of GWID `%s` and healpix file `%s` already exists in the db.************ \nClobbering...")
                # Delete the map in the db...

                create_baks = 'CALL BackupTables(%s);' % current_map_id
                lock_tables = '''
                            LOCK TABLE HealpixMap WRITE, HealpixPixel WRITE, HealpixPixel_Completeness WRITE, HealpixPixel_Galaxy_Weight WRITE, 
                            ObservedTile WRITE, ObservedTile_HealpixPixel WRITE, StaticTile_HealpixPixel WRITE, HealpixMap_bak WRITE, HealpixPixel_bak WRITE, 
                            HealpixPixel_Completeness_bak WRITE, HealpixPixel_Galaxy_Weight_bak WRITE, ObservedTile_bak WRITE, ObservedTile_HealpixPixel_bak WRITE, 
                            StaticTile_HealpixPixel_bak WRITE;
                        '''
                delete_map = 'CALL DeleteMap(%s);' % current_map_id
                unlock_tables = 'UNLOCK TABLES;'

                query_db([create_baks], commit=True)
                query_db([lock_tables], commit=True)
                query_db([delete_map], commit=True)
                query_db([unlock_tables], commit=True)

                print("\n************ Map ID `%s` has been deleted, proceeding with load_map...************\n" % current_map_id)

            else:
                is_error = True
                print("The combination of GWID `%s` and healpix file `%s` already exists in the db. Please choose a "
                      "unique combination." % (self.options.gw_id, self.options.healpix_file))

        if self.options.skip_swope:
            print("**** Not registering Swope tiles for this event! ****")

        if self.options.skip_thacher:
            print("**** Not registering Thacher tiles for this event! ****")

        # Only build tile-pixel relations if we don't skip both Swope AND Thacher
        if self.options.skip_swope and self.options.skip_thacher:
            build_tile_pixel_relation = False

        if not build_tile_pixel_relation:
            print("**** Skipping tile-pixel relations entirely! ****")

        if is_error:
            print("Exiting...")
            return 1

        ### Quantities that are required ###
        nside128 = 128
        nside256 = 256
        prob = None
        distmu = None
        distsigma = None
        distnorm = None
        header_gen = None
        map_nside = None
        healpix_map_id = None
        N128_dict = None
        map_pixel_dict = None

        # Can I remove this?
        if not build_map:
            print("Skipping build map...")
            print("\tLoading existing map from disk...")
            orig_prob, orig_distmu, orig_distsigma, orig_distnorm, header_gen = hp.read_map(hpx_path, field=(0, 1, 2, 3), h=True)

            orig_npix = len(orig_prob)
            print("\tOriginal number of pix in '%s': %s" % (self.options.healpix_file, orig_npix))

            orig_map_nside = hp.npix2nside(orig_npix)

            # Check if NSIDE is > 256. If it is, and the --orig_res flag is not specified, rescale map to NSIDE 256
            if orig_map_nside > nside256 and not self.options.analysis_mode:
                print("Rescaling map to NSIDE = 256")
                rescaled_prob, rescaled_distmu, rescaled_distsigma, rescaled_distnorm = hp.ud_grade([orig_prob,
                                                                                                    orig_distmu,
                                                                                                    orig_distsigma,
                                                                                                    orig_distnorm],
                                                                                                    nside_out=nside256,
                                                                                                    order_in="RING",
                                                                                                    order_out="RING")
                rescaled_npix = len(rescaled_prob)
                print("\tRescaled number of pix in '%s': %s" % (self.options.healpix_file, rescaled_npix))

                rescaled_map_nside = hp.npix2nside(rescaled_npix)
                print("\tRescaled resolution (nside) of '%s': %s\n" % (self.options.healpix_file, rescaled_map_nside))

                original_pix_per_rescaled_pix = orig_npix / rescaled_npix
                print("Original pix per rescaled pix for %s" % original_pix_per_rescaled_pix)

                print("Renormalizing and initializing rescaled map...")
                rescaled_prob = rescaled_prob * original_pix_per_rescaled_pix

                prob = rescaled_prob
                distmu = rescaled_distmu
                distsigma = rescaled_distsigma
                distnorm = rescaled_distnorm
                map_nside = rescaled_map_nside
            else:
                print("##### Loading ORIGINAL RESOLUTION MAP @ NSIDE = %s #####" % orig_map_nside)
                prob = orig_prob
                distmu = orig_distmu
                distsigma = orig_distsigma
                distnorm = orig_distnorm
                map_nside = orig_map_nside

            print("Getting map id")
            healpix_map_id = query_db([healpix_map_select % (self.options.gw_id, self.options.healpix_file)])[0][0][0]
            print("Done map id")
        else:
            print("Building map...")
            # Create target directory if it doesn't already exist...
            if not os.path.exists(formatted_healpix_dir):
                os.mkdir(formatted_healpix_dir)
                print("\n\nDirectory ", formatted_healpix_dir, " Created ")
            else:
                print("\n\nDirectory ", formatted_healpix_dir, " already exists")

            if not os.path.exists(formatted_candidates_dir):
                os.mkdir(formatted_candidates_dir)
                print("\n\nDirectory ", formatted_candidates_dir, " Created ")
            else:
                print("\n\nDirectory ", formatted_candidates_dir, " already exists")

            if not os.path.exists(formatted_obs_tiles_dir):
                os.mkdir(formatted_obs_tiles_dir)
                print("\n\nDirectory ", formatted_obs_tiles_dir, " Created ")
            else:
                print("\n\nDirectory ", formatted_obs_tiles_dir, " already exists")

            if not os.path.exists(formatted_model_dir):
                os.mkdir(formatted_model_dir)
                print("\n\nDirectory ", formatted_model_dir, " Created ")
            else:
                print("\n\nDirectory ", formatted_model_dir, " already exists")

            # try:
            #     os.mkdir(formatted_healpix_dir)
            #     print("\n\nDirectory ", formatted_healpix_dir, " Created ")
            #
            #     # Create the `ObservedTiles` and `Candidates` subdirectories
            #     os.mkdir(formatted_candidates_dir)
            #     os.mkdir(formatted_obs_tiles_dir)
            #     os.mkdir(formatted_model_dir)
            # except FileExistsError:
            #     print("\n\nDirectory ", formatted_healpix_dir, " already exists")

            # if not self.options.load_local_file:
            # Get file -- ADD check and only get if you need to...

            t_0 = ""
            t_0_datetime_str = ""
            gw_url = ""

            try:
                t1 = time.time()

                gdb_client = GraceDb(api_endpoint)
                event_json = gdb_client.superevent(self.options.gw_id).json()
                far = event_json["far"]
                t_0 = event_json["t_0"]
                t_0_time_obj = Time(t_0, scale="utc", format="gps")
                t_0_datetime_str = t_0_time_obj.to_value(format="iso")
                gw_url = event_json["links"]["self"]

                # Download and save map file.
                # The default self.options.healpix_file - "bayestar.fits.gz" always has the most up-to-date flat file
                #   in a search context
                file_response = gdb_client.files(self.options.gw_id, self.options.healpix_file)
                with open(hpx_path, "wb") as f:
                    f.write(file_response.data)

                # print("Downloading `%s`..." % self.options.healpix_file)
                # t1 = time.time()
                # gw_file_formatter = "https://gracedb.ligo.org/apiweb/superevents/%s/files/%s"
                # gw_file = gw_file_formatter % (self.options.gw_id, self.options.healpix_file)
                #
                # with urllib.request.urlopen(gw_file) as response:
                #     with open(hpx_path, 'wb') as out_file:
                #         shutil.copyfileobj(response, out_file)


                t2 = time.time()

                print("\n********* start DEBUG ***********")
                print("Downloading `%s` - execution time: %s" % (self.options.healpix_file, (t2 - t1)))
                print("********* end DEBUG ***********\n")

            except urllib.error.HTTPError as e:
                print("\nError:")
                print(e.code, self.options.healpix_file)
                print("\n\tExiting...")
                return 1
            except urllib.error.URLError as e:
                print("\nError:")
                if hasattr(e, 'reason'):
                    print(e.reason, self.options.healpix_file)
                elif hasattr(e, 'code'):
                    print(e.code, self.options.healpix_file)
                print("\n\tExiting...")
                return 1
            except Error as e:
                print("\nError:")
                print(e)
                print("\n\tExiting...")
                return 1

            # Get event details from GraceDB page
            # try:
            #     print("Scraping `%s` details..." % self.options.gw_id)
            #     t1 = time.time()
            #     gw_url_formatter = "https://gracedb.ligo.org/superevents/%s/view/"
            #     gw_url = gw_url_formatter % self.options.gw_id
            #
            #     r = requests.get(gw_url)
            #     data = r.text
            #     soup = BeautifulSoup(data, "lxml")
            #
            #     ## Current HTML is different. HACK: adding hard coded values so that this just works
            #     # tds = soup.find("table", {"class": "superevent"}).find_all('tr')[1].find_all('td')
            #     # FAR_hz = tds[3].text
            #     # FAR_per_yr = tds[4].text
            #     # t_start = tds[5].text
            #     # t_0 = tds[6].text
            #     # t_end = tds[7].text
            #     # UTC_submission_time = parse(tds[8].text)
            #
            #     # 0814 values
            #     # t_0 = 1249852257.012957
            #     # UTC_submission_time = parse("2019-08-14 21:11:18")
            #
            #     # 0425 values
            #     t_0 = 1240215503.01
            #     UTC_submission_time = parse("2019-04-25 08:18:26")
            #
            #
            #     t2 = time.time()
            #     print("\n********* start DEBUG ***********")
            #     print("Scraping `%s` details - execution time: %s" % (self.options.gw_id, (t2 - t1)))
            #     print("********* end DEBUG ***********\n")
            #
            # except urllib.error.HTTPError as e:
            #     print("\nError:")
            #     print(e.code, self.options.healpix_file)
            #     print("\n\tExiting...")
            #     return 1
            # except urllib.error.URLError as e:
            #     print("\nError:")
            #     if hasattr(e, 'reason'):
            #         print(e.reason, self.options.healpix_file)
            #     elif hasattr(e, 'code'):
            #         print(e.code, self.options.healpix_file)
            #     print("\n\tExiting...")
            #     return 1
            # except Error as e:
            #     print("\nError:")
            #     print(e)
            #     print("\n\tExiting...")
            #     return 1

            # Unpack healpix file and insert map into db
            print("Unpacking '%s':%s..." % (self.options.gw_id, self.options.healpix_file))
            t1 = time.time()

            (orig_prob, orig_distmu, orig_distsigma, orig_distnorm), header_gen = hp.read_map(hpx_path, field=(0, 1, 2, 3), h=True)

            header = dict(header_gen)
            orig_npix = len(orig_prob)
            print("\tOriginal number of pix in '%s': %s" % (self.options.healpix_file, orig_npix))

            sky_area = 4 * 180 ** 2 / np.pi
            area_per_orig_px = sky_area / orig_npix
            print("\tSky Area per original pix in '%s': %s [sq deg]" % (self.options.healpix_file, area_per_orig_px))

            orig_map_nside = hp.npix2nside(orig_npix)
            print("\tOriginal resolution (nside) of '%s': %s\n" % (self.options.healpix_file, orig_map_nside))


            # Check if NSIDE is > 256. If it is, and the --orig_res flag is not specified, rescale map to NSIDE 256
            if orig_map_nside > nside256 and not self.options.analysis_mode:
                print("Rescaling map to NSIDE = 256")
                rescaled_prob, rescaled_distmu, rescaled_distsigma, rescaled_distnorm = hp.ud_grade([orig_prob,
                                                                                                    orig_distmu,
                                                                                                    orig_distsigma,
                                                                                                    orig_distnorm],
                                                                                                    nside_out=nside256,
                                                                                                    order_in="RING",
                                                                                                    order_out="RING")
                rescaled_npix = len(rescaled_prob)
                print("\tRescaled number of pix in '%s': %s" % (self.options.healpix_file, rescaled_npix))

                area_per_rescaled_px = sky_area / rescaled_npix
                print("\tSky Area per rescaled pix in '%s': %s [sq deg]" % (self.options.healpix_file, area_per_rescaled_px))

                rescaled_map_nside = hp.npix2nside(rescaled_npix)
                print("\tRescaled resolution (nside) of '%s': %s\n" % (self.options.healpix_file, rescaled_map_nside))

                original_pix_per_rescaled_pix = orig_npix / rescaled_npix
                print("Original pix per rescaled pix for %s" % original_pix_per_rescaled_pix)

                print("Renormalizing and initializing rescaled map...")
                rescaled_prob = rescaled_prob * original_pix_per_rescaled_pix

                prob = rescaled_prob
                distmu = rescaled_distmu
                distsigma = rescaled_distsigma
                distnorm = rescaled_distnorm
                map_nside = rescaled_map_nside
            else:
                print("##### Loading ORIGINAL RESOLUTION MAP @ NSIDE = %s #####" % orig_map_nside)

                prob = orig_prob
                distmu = orig_distmu
                distsigma = orig_distsigma
                distnorm = orig_distnorm
                map_nside = orig_map_nside

            # Initialize with 0.0 prob to galaxies
            healpix_map_insert = '''INSERT INTO HealpixMap 
            (GWID, URL, Filename, NSIDE, t_0, SubmissionTime, NetProbToGalaxies, RescaledNSIDE) 
            VALUES (%s,%s,%s,%s,%s,%s,%s,%s);'''

            # healpix_map_data = [(self.options.gw_id, gw_url, self.options.healpix_file, orig_map_nside, t_0,
            #                      UTC_submission_time.strftime('%Y-%m-%d %H:%M:%S'), 0.0, map_nside)]
            healpix_map_data = [(self.options.gw_id, gw_url, self.options.healpix_file, orig_map_nside, t_0,
                                 t_0_datetime_str, 0.0, map_nside)]

            print("\nInserting %s Healpix Map: %s ..." % (self.options.gw_id, self.options.healpix_file))
            insert_records(healpix_map_insert, healpix_map_data)
            healpix_map_id = query_db([healpix_map_select % (self.options.gw_id, self.options.healpix_file)])[0][0][0]
            print("...Done")

        if not build_pixels:
            print("Skipping pixels...")
            print("\tLoading NSIDE 128 pixels...")
            with open(pickle_output_dir + 'N128_dict.pkl', 'rb') as handle:
                N128_dict = pickle.load(handle)
                del handle

            t1 = time.time()

            print("Get map pixel")
            map_pixel_select = '''SELECT 
                id, HealpixMap_id, Pixel_Index, Prob, Distmu, Distsigma, 
                Distnorm, Mean, Stddev, Norm, N128_SkyPixel_id 
            FROM HealpixPixel WHERE HealpixMap_id = %s;'''

            map_pixels = query_db([map_pixel_select % healpix_map_id])[0]

            map_pixel_dict = {}
            for p in map_pixels:
                map_pixel_dict[int(p[2])] = p

            # clean up
            print("freeing `map_pixels`...")
            del map_pixels

            process = psutil.Process(os.getpid())
            mb = process.memory_info().rss / 1e+6
            print("Total mem usage: %0.3f [MB]" % mb)
            print("\n*****************************")

            t2 = time.time()
            print("\n********* start DEBUG ***********")
            print("Pixel select execution time: %s" % (t2 - t1))
            print("********* end DEBUG ***********\n")
        else:
            print("Building pixels...")
            # Process healpix pixels and bulk insert into db. To do this, `LOAD DATA LOCAL INFILE` must be enabled in MySQL.
            # The pixels are unpacked, associated with NSIDE 128 pixels, their LIGO distance distributions are resolved,
            # and finally they are written to the working directory as a CSV that is bulk uploaded to the db. Clean up the file
            # afterward.

            healpix_map_select = "SELECT id FROM HealpixMap WHERE GWID = '%s' and Filename = '%s'"
            healpix_map_id = query_db([healpix_map_select % (self.options.gw_id, self.options.healpix_file)])[0][0][0]

            # Associate these healpix pixels with our N128 Sky Pixels
            t1 = time.time()
            N128_select = "SELECT id, Pixel_Index FROM SkyPixel WHERE NSIDE = 128;"
            N128_result = query_db([N128_select])[0]
            print("Number of NSIDE 128 sky pixels: %s" % len(N128_result))

            N128_dict = {}
            for n128 in N128_result:
                N128_dict[int(n128[1])] = int(n128[0])

            # N128_dict = {}
            # with open(pickle_output_dir + 'N128_dict.pkl', 'rb') as handle:
            #     N128_dict = pickle.load(handle)
            #     del handle

            # Clean up N128_result
            print("freeing `N128_result`...")
            del N128_result

            t2 = time.time()
            print("\n********* start DEBUG ***********")
            print("N128 select execution time: %s" % (t2 - t1))
            print("********* end DEBUG ***********\n")

            # Resolve pixel distance distribution. Clean output arrays prior to serialization...
            # MySQL cannot store the np.inf value returned by healpy. Get a value that's close but not too close!!
            t1 = time.time()
            max_double_value = 1.5E+308

            from ligo.skymap.postprocess.util import find_greedy_credible_levels
            percentiles = find_greedy_credible_levels(prob)
            _90th_indices = np.where(percentiles <= 0.91)

            windowed_prob = prob[_90th_indices]
            windowed_distmu = distmu[_90th_indices]
            windowed_distsigma = distsigma[_90th_indices]
            windowed_distnorm = distnorm[_90th_indices]

            # distmu[distmu > max_double_value] = max_double_value
            # distsigma[distsigma > max_double_value] = max_double_value
            # distnorm[distnorm > max_double_value] = max_double_value
            #
            # mean, stddev, norm = distance.parameters_to_moments(distmu, distsigma)
            #
            # mean[mean > max_double_value] = max_double_value
            # stddev[stddev > max_double_value] = max_double_value
            # norm[norm > max_double_value] = max_double_value

            windowed_distmu[windowed_distmu > max_double_value] = max_double_value
            windowed_distsigma[windowed_distsigma > max_double_value] = max_double_value
            windowed_distnorm[windowed_distnorm > max_double_value] = max_double_value

            mean, stddev, norm = distance.parameters_to_moments(windowed_distmu, windowed_distsigma)

            mean[mean > max_double_value] = max_double_value
            stddev[stddev > max_double_value] = max_double_value
            norm[norm > max_double_value] = max_double_value

            # theta, phi = hp.pix2ang(map_nside, range(len(prob)))
            theta, phi = hp.pix2ang(map_nside, _90th_indices)
            N128_indices = hp.ang2pix(nside128, theta, phi)

            healpix_pixel_data = []
            # for i, n128_i in enumerate(N128_indices):
            for i, (pix_i, n128_i) in enumerate(zip(_90th_indices[0], N128_indices[0])):
                N128_pixel_id = N128_dict[n128_i]
                # healpix_pixel_data.append((healpix_map_id, i, prob[i], distmu[i], distsigma[i], distnorm[i], mean[i],
                #                            stddev[i], norm[i], N128_pixel_id))
                healpix_pixel_data.append((healpix_map_id, pix_i, windowed_prob[i], windowed_distmu[i],
                                           windowed_distsigma[i], windowed_distnorm[i], mean[i],
                                           stddev[i], norm[i], N128_pixel_id))

            # Clean up
            print("freeing `theta`...")
            print("freeing `phi`...")
            print("freeing `N128_indices`...")
            del theta
            del phi
            del N128_indices

            t2 = time.time()
            print("\n********* start DEBUG ***********")
            print("Distance resolution execution time: %s" % (t2 - t1))
            print("********* end DEBUG ***********\n")

            # Create CSV, upload, and clean up CSV
            upload_csv = "%s/%s_bulk_upload.csv" % (formatted_healpix_dir, self.options.gw_id)
            try:
                t1 = time.time()
                print("Creating `%s`" % upload_csv)
                with open(upload_csv, 'w') as csvfile:
                    csvwriter = csv.writer(csvfile)
                    for data in healpix_pixel_data:
                        csvwriter.writerow(data)

                t2 = time.time()
                print("\n********* start DEBUG ***********")
                print("CSV creation execution time: %s" % (t2 - t1))
                print("********* end DEBUG ***********\n")
            except Error as e:
                print("Error in creating CSV:\n")
                print(e)
                print("\nExiting")
                return 1

            t1 = time.time()
            upload_sql = """LOAD DATA LOCAL INFILE '%s' 
                    INTO TABLE HealpixPixel 
                    FIELDS TERMINATED BY ',' 
                    LINES TERMINATED BY '\n' 
                    (HealpixMap_id, Pixel_Index, Prob, Distmu, Distsigma, Distnorm, Mean, Stddev, Norm, N128_SkyPixel_id);"""

            success = bulk_upload(upload_sql % upload_csv)
            if not success:
                print("\nUnsuccessful bulk upload. Exiting...")
                return 1

            t2 = time.time()
            print("\n********* start DEBUG ***********")
            print("CSV upload execution time: %s" % (t2 - t1))

            try:
                print("Removing `%s`..." % upload_csv)
                os.remove(upload_csv)

                # Clean up
                print("freeing `healpix_pixel_data`...")
                del healpix_pixel_data

                print("... Done")
            except Error as e:
                print("Error in file removal")
                print(e)
                print("\nExiting")
                return 1

            t1 = time.time()
            map_pixel_select = '''SELECT 
                id, HealpixMap_id, Pixel_Index, Prob, Distmu, Distsigma, 
                Distnorm, Mean, Stddev, Norm, N128_SkyPixel_id 
            FROM HealpixPixel WHERE HealpixMap_id = %s'''

            map_pixels = query_db([map_pixel_select % healpix_map_id])[0]
            print(len(map_pixels))
            map_pixel_dict = {}
            for p in map_pixels:
                map_pixel_dict[int(p[2])] = p

            # Clean up
            print("freeing `map_pixels`...")
            del map_pixels

            t2 = time.time()
            print("\n********* start DEBUG ***********")
            print("Pixel select execution time: %s" % (t2 - t1))
            print("********* end DEBUG ***********\n")

            with open(pickle_output_dir + 'N128_dict.pkl', 'wb') as handle:
                pickle.dump(N128_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

        if not build_tile_pixel_relation:
            print("Skipping tile-pixel relation...")
        else:
            print("Building tile-pixel relation...")

            # select_detector_id = "SELECT id, Name, Deg_width, Deg_height, Deg_radius, Area FROM Detector WHERE Name='%s'"
            select_detector = "SELECT id, Name, ST_AsText(Poly) FROM Detector WHERE Name='%s';"
            tile_select = "SELECT id, Detector_id, FieldName, RA, _Dec, Coord, Poly, EBV, N128_SkyPixel_id FROM StaticTile WHERE Detector_id = %s;"
            tile_pixel_upload_csv = "%s/%s_tile_pixel_upload.csv" % (formatted_healpix_dir, self.options.gw_id)

            tile_pixel_data = []
            if not self.options.skip_swope:
                ##### DO SWOPE ######
                # Get detector -> static tile rows
                swope_detector_result = query_db([select_detector % "SWOPE"])[0][0]
                swope_id = swope_detector_result[0]
                swope_name = swope_detector_result[1]
                swope_poly = swope_detector_result[2]
                swope_detector_vertices = Detector.get_detector_vertices_from_teglon_db(swope_poly)
                swope_detector = Detector(swope_name, swope_detector_vertices, detector_id=swope_id)
                swope_static_tile_rows = query_db([tile_select % swope_id])[0]

                swope_tiles = []
                for r in swope_static_tile_rows:

                    t = Tile(central_ra_deg=float(r[3]), central_dec_deg=float(r[4]), detector=swope_detector,
                             nside=map_nside, tile_id=int(r[0]), tile_mwe=float(r[7]))
                    swope_tiles.append(t)

                # clean up
                print("freeing `swope_static_tile_rows`...")
                del swope_static_tile_rows

                t1 = time.time()
                initialized_swope_tiles = None
                with mp.Pool() as pool:
                    initialized_swope_tiles = pool.map(initialize_tile, swope_tiles)

                # clean up
                print("freeing `swope_tiles`...")
                del swope_tiles
                t2 = time.time()

                print("\n********* start DEBUG ***********")
                print("Swope Tile initialization execution time: %s" % (t2 - t1))
                print("********* end DEBUG ***********\n")

                # Insert Tile/Healpix pixel relations
                for t in initialized_swope_tiles:
                    for p in t.enclosed_pixel_indices:
                        if p in map_pixel_dict:
                            tile_pixel_data.append((t.id, map_pixel_dict[p][0], healpix_map_id))

                # Create CSV
                try:
                    t1 = time.time()
                    print("Appending `%s`" % tile_pixel_upload_csv)
                    with open(tile_pixel_upload_csv, 'a') as csvfile:
                        csvwriter = csv.writer(csvfile)
                        for data in tile_pixel_data:
                            csvwriter.writerow(data)

                    t2 = time.time()
                    print("\n********* start DEBUG ***********")
                    print("Tile-Pixel CSV creation execution time: %s" % (t2 - t1))
                    print("********* end DEBUG ***********\n")
                except Error as e:
                    print("Error in creating Tile-Pixel CSV:\n")
                    print(e)
                    print("\nExiting")
                    return 1

                # clean up
                print("freeing `tile_pixel_data`...")
                del tile_pixel_data

            tile_pixel_data = []
            if not self.options.skip_thacher:
                ##### DO THACHER ######
                # Get detector -> static tile rows
                thacher_detector_result = query_db([select_detector % "THACHER"])[0][0]
                thacher_id = thacher_detector_result[0]
                thacher_name = thacher_detector_result[1]
                thacher_poly = thacher_detector_result[2]
                thacher_detector_vertices = Detector.get_detector_vertices_from_teglon_db(thacher_poly)
                thacher_detector = Detector(thacher_name, thacher_detector_vertices, detector_id=thacher_id)
                thacher_static_tile_rows = query_db([tile_select % thacher_id])[0]

                thacher_tiles = []
                for r in thacher_static_tile_rows:
                    t = Tile(central_ra_deg=float(r[3]), central_dec_deg=float(r[4]), detector=thacher_detector,
                             nside=map_nside, tile_id=int(r[0]), tile_mwe=float(r[7]))
                    thacher_tiles.append(t)

                # clean up
                print("freeing `thacher_static_tile_rows`...")
                del thacher_static_tile_rows

                t1 = time.time()
                initialized_thacher_tiles = None
                with mp.Pool() as pool:
                    initialized_thacher_tiles = pool.map(initialize_tile, thacher_tiles)

                # clean up
                print("freeing `thacher_tiles`...")
                del thacher_tiles
                t2 = time.time()

                print("\n********* start DEBUG ***********")
                print("Thacher Tile initialization execution time: %s" % (t2 - t1))
                print("********* end DEBUG ***********\n")

                # Insert Tile/Healpix pixel relations
                for t in initialized_thacher_tiles:
                    for p in t.enclosed_pixel_indices:
                        if p in map_pixel_dict:
                            tile_pixel_data.append((t.id, map_pixel_dict[p][0], healpix_map_id))

                # Append to existing CSV, upload, and clean up CSV
                try:
                    t1 = time.time()
                    print("Appending `%s`" % tile_pixel_upload_csv)
                    with open(tile_pixel_upload_csv, 'a') as csvfile:
                        csvwriter = csv.writer(csvfile)
                        for data in tile_pixel_data:
                            csvwriter.writerow(data)

                    t2 = time.time()
                    print("\n********* start DEBUG ***********")
                    print("Tile-Pixel CSV append execution time: %s" % (t2 - t1))
                    print("********* end DEBUG ***********\n")
                except Error as e:
                    print("Error in creating Tile-Pixel CSV:\n")
                    print(e)
                    print("\nExiting")
                    return 1

            tile_pixel_data = []
            if not self.options.skip_t80:
                ##### DO T80 ######
                # Get detector -> static tile rows
                t80_detector_result = query_db([select_detector % "T80S_T80S-Cam"])[0][0]
                t80_id = t80_detector_result[0]
                t80_name = t80_detector_result[1]
                t80_poly = t80_detector_result[2]
                t80_detector_vertices = Detector.get_detector_vertices_from_teglon_db(t80_poly)
                t80_detector = Detector(t80_name, t80_detector_vertices, detector_id=t80_id)
                t80_static_tile_rows = query_db([tile_select % t80_id])[0]

                t80_tiles = []
                for r in t80_static_tile_rows:
                    t = Tile(central_ra_deg=float(r[3]), central_dec_deg=float(r[4]), detector=t80_detector,
                             nside=map_nside, tile_id=int(r[0]), tile_mwe=float(r[7]))
                    t80_tiles.append(t)

                # clean up
                print("freeing `t80_static_tile_rows`...")
                del t80_static_tile_rows

                t1 = time.time()
                initialized_t80_tiles = None
                with mp.Pool() as pool:
                    initialized_t80_tiles = pool.map(initialize_tile, t80_tiles)

                # clean up
                print("freeing `t80_tiles`...")
                del t80_tiles
                t2 = time.time()

                print("\n********* start DEBUG ***********")
                print("T80 Tile initialization execution time: %s" % (t2 - t1))
                print("********* end DEBUG ***********\n")

                # Insert Tile/Healpix pixel relations
                for t in initialized_t80_tiles:
                    for p in t.enclosed_pixel_indices:
                        if p in map_pixel_dict:
                            tile_pixel_data.append((t.id, map_pixel_dict[p][0], healpix_map_id))

                # Append to existing CSV, upload, and clean up CSV
                try:
                    t1 = time.time()
                    print("Appending `%s`" % tile_pixel_upload_csv)
                    with open(tile_pixel_upload_csv, 'a') as csvfile:
                        csvwriter = csv.writer(csvfile)
                        for data in tile_pixel_data:
                            csvwriter.writerow(data)

                    t2 = time.time()
                    print("\n********* start DEBUG ***********")
                    print("Tile-Pixel CSV append execution time: %s" % (t2 - t1))
                    print("********* end DEBUG ***********\n")
                except Error as e:
                    print("Error in creating Tile-Pixel CSV:\n")
                    print(e)
                    print("\nExiting")
                    return 1

        print("Bulk uploading Tile-Pixel...")
        t1 = time.time()
        st_hp_upload_sql = """LOAD DATA LOCAL INFILE '%s' 
                            INTO TABLE StaticTile_HealpixPixel 
                            FIELDS TERMINATED BY ',' 
                            LINES TERMINATED BY '\n' 
                            (StaticTile_id, HealpixPixel_id, HealpixMap_id);"""

        success = bulk_upload(st_hp_upload_sql % tile_pixel_upload_csv)
        if not success:
            print("\nUnsuccessful bulk upload. Exiting...")
            return 1

        t2 = time.time()
        print("\n********* start DEBUG ***********")
        print("Tile-Pixel CSV upload execution time: %s" % (t2 - t1))

        try:
            print("Removing `%s`..." % tile_pixel_upload_csv)
            os.remove(tile_pixel_upload_csv)

            # clean up
            print("freeing `tile_pixel_data`...")
            del tile_pixel_data

            print("... Done")
        except Error as e:
            print("Error in file removal")
            print(e)
            print("\nExiting")
            return 1

        if not build_galaxy_pixel_relation:
            print("Skipping galaxy-pixel relation...")
        else:
            pass
            # print("Building galaxy-pixel relation...")
            # # 1. Select all Galaxies from GD2
            # # 2. Resolve all Galaxies to a HealpixPixel index/id
            # # 3. Store HealpixPixel_Galaxy association
            # t1 = time.time()
            # galaxy_select = '''
            #     SELECT 	id,
            #             PGC,
            #             Name_GWGC,
            #             Name_HyperLEDA,
            #             Name_2MASS,
            #             Name_SDSS_DR12,
            #             RA,
            #             _Dec,
            #             Coord,
            #             dist,
            #             dist_err,
            #             z_dist,
            #             z_dist_err,
            #             z,
            #             B,
            #             B_err,
            #             B_abs,
            #             J,
            #             J_err,
            #             H,
            #             H_err,
            #             K,
            #             K_err,
            #             flag1,
            #             flag2,
            #             flag3,
            #             N2_Pixel_Index,
            #             N4_Pixel_Index,
            #             N8_Pixel_Index,
            #             N16_Pixel_Index,
            #             N32_Pixel_Index,
            #             N64_Pixel_Index,
            #             N128_Pixel_Index,
            #             N256_Pixel_Index,
            #             N512_Pixel_Index,
            #             N1024_Pixel_Index,
            #             %s
            #     FROM Galaxy
            #     WHERE z_dist < 1206.0
            # '''
            # galaxy_result = query_db([galaxy_select % map_nside])[0]
            # print("Number of Galaxies: %s" % len(galaxy_result))
            # t2 = time.time()
            # print("\n********* start DEBUG ***********")
            # print("Galaxy select execution time: %s" % (t2 - t1))
            # print("********* end DEBUG ***********\n")
            #
            # t1 = time.time()
            #
            # with mp.Pool() as pool:
            #     galaxy_pixel_relations = pool.map(get_healpix_pixel_id, galaxy_result)
            #
            # galaxy_pixel_data = []
            # for gp in galaxy_pixel_relations:
            #     _pixel_index = gp[1]
            #     _pixel_id = map_pixel_dict[_pixel_index][0]
            #     _galaxy_id = gp[0]
            #     galaxy_pixel_data.append((_pixel_id, _galaxy_id))
            #
            # t2 = time.time()
            # print("\n********* start DEBUG ***********")
            # print("Galaxy-pixel creation execution time: %s" % (t2 - t1))
            # print("********* end DEBUG ***********\n")
            #
            # # Create CSV, upload, and clean up CSV
            # galaxy_pixel_upload_csv = "%s/%s_gal_pix_upload.csv" % (formatted_healpix_dir, self.options.gw_id)
            # try:
            #     t1 = time.time()
            #     print("Creating `%s`" % galaxy_pixel_upload_csv)
            #     with open(galaxy_pixel_upload_csv, 'w') as csvfile:
            #         csvwriter = csv.writer(csvfile)
            #         for data in galaxy_pixel_data:
            #             csvwriter.writerow(data)
            #
            #     t2 = time.time()
            #     print("\n********* start DEBUG ***********")
            #     print("CSV creation execution time: %s" % (t2 - t1))
            #     print("********* end DEBUG ***********\n")
            # except Error as e:
            #     print("Error in creating CSV:\n")
            #     print(e)
            #     print("\nExiting")
            #     return 1
            #
            # t1 = time.time()
            # upload_sql = """LOAD DATA LOCAL INFILE '%s'
            #         INTO TABLE HealpixPixel_Galaxy
            #         FIELDS TERMINATED BY ','
            #         LINES TERMINATED BY '\n'
            #         (HealpixPixel_id, Galaxy_id);"""
            #
            # success = bulk_upload(upload_sql % galaxy_pixel_upload_csv)
            # if not success:
            #     print("\nUnsuccessful bulk upload. Exiting...")
            #     return 1
            #
            # t2 = time.time()
            # print("\n********* start DEBUG ***********")
            # print("CSV upload execution time: %s" % (t2 - t1))
            #
            # try:
            #     print("Removing `%s`..." % galaxy_pixel_upload_csv)
            #     os.remove(galaxy_pixel_upload_csv)
            #
            #     # Clean up
            #     print("freeing `galaxy_pixel_data`...")
            #     del galaxy_pixel_data
            #
            #     print("... Done")
            # except Error as e:
            #     print("Error in file removal")
            #     print(e)
            #     print("\nExiting")
            #     return 1

        # clean up
        print("freeing `prob`...")
        print("freeing `distmu`...")
        print("freeing `distsigma`...")
        print("freeing `distnorm`...")
        print("freeing `header_gen`...")
        # print("freeing `map_pixel_dict`...")
        del prob
        del distmu
        del distsigma
        del distnorm
        del header_gen
        # del map_pixel_dict

        if not build_completeness_func:
            print("Skipping completeness func...")
        # else:
        elif False:
            print("Building completeness func...")

            pixel_completeness_upload_csv = "%s/%s_pixel_completeness_upload.csv" % (
            formatted_healpix_dir, self.options.gw_id)

            completeness_select = '''
                SELECT 
                    sp1.id as N128_id, 0.5*(sd.D1+sd.D2) as Dist, sc.SmoothedCompleteness
                FROM SkyPixel sp1 
                JOIN SkyPixel sp2 on sp2.id = sp1.Parent_Pixel_id 
                JOIN SkyPixel sp3 on sp3.id = sp2.Parent_Pixel_id 
                JOIN SkyPixel sp4 on sp4.id = sp3.Parent_Pixel_id 
                JOIN SkyPixel sp5 on sp5.id = sp4.Parent_Pixel_id 
                JOIN SkyPixel sp6 on sp6.id = sp5.Parent_Pixel_id 
                JOIN SkyPixel sp7 on sp7.id = sp6.Parent_Pixel_id 
                JOIN SkyCompleteness sc on sc.SkyPixel_id in (sp1.id, sp2.id, sp3.id, sp4.id, sp5.id, sp6.id, sp7.id) 
                JOIN SkyDistance sd on sd.id = sc.SkyDistance_id 
                WHERE sp1.id in (%s)
                ORDER BY sp1.id, sd.D1 
            '''

            healpix_pixel_completeness_select = '''
                SELECT 
                    id, 
                    HealpixMap_id, 
                    Pixel_Index, 
                    Prob, 
                    Distmu, 
                    Distsigma, 
                    Distnorm, 
                    Mean, 
                    Stddev, 
                    Norm, 
                    N128_SkyPixel_id
                FROM HealpixPixel hp
                WHERE HealpixMap_id = %s and hp.N128_SkyPixel_id in (%s)
            '''

            # get clusters of n128 by id
            # get completeness where n128 id in (above cluster)
            # match healpix pixels to n128 ids, and compute completeness and push changes
            # repeat until done with clusters

            # Create list from dict so we can sort
            n128_pixel_id_list = [pixel_id for pixel_index, pixel_id in
                                  N128_dict.items()]  # pixel_index, pixel_id are both INT

            # clean up
            print("freeing `N128_dict`...")
            del N128_dict

            batch_size = 10000
            i = 0
            j = batch_size
            k = len(n128_pixel_id_list)

            process = psutil.Process(os.getpid())
            mb = process.memory_info().rss / 1e+6
            print("\nTotal mem usage: %0.3f [MB]\n" % mb)
            print("\nLength of N128 pixels: %s" % k)
            print("Batch size: %s" % batch_size)
            print("Starting loop...")

            number_of_iters = k // batch_size
            if k % batch_size > 0:
                number_of_iters += 1

            iter_num = 1
            while j < k:

                t1 = time.time()

                print("%s:%s" % (i, j))
                batch = n128_pixel_id_list[i:j]
                batch_id_string = ','.join(map(str, batch))

                completeness_result = query_db([completeness_select % batch_id_string])[0]
                print("Number of completeness records: %s" % len(completeness_result))
                healpix_pixel_completeness_result = \
                query_db([healpix_pixel_completeness_select % (healpix_map_id, batch_id_string)])[0]
                print("Number of healpix pixel records: %s" % len(healpix_pixel_completeness_result))

                i = j
                j += batch_size

                # DO WORK HERE
                completeness_values_dict = {}
                for c in completeness_result:
                    N128_id = int(c[0])

                    if N128_id not in completeness_values_dict:
                        completeness_values_dict[N128_id] = [[0.0], [1.0]]

                    comp = float(c[2]) if float(c[2]) <= 1.0 else 1.0
                    completeness_values_dict[N128_id][0].append(float(c[1]))
                    completeness_values_dict[N128_id][1].append(comp)

                # Add point at infinity
                for key, val in completeness_values_dict.items():
                    completeness_values_dict[key][0].append(np.inf)
                    completeness_values_dict[key][1].append(0.0)

                # clean up
                print("freeing `completeness_result`...")
                del completeness_result

                pixel_completeness_records = []
                for p in healpix_pixel_completeness_result:
                    pix_id = p[0]
                    pix_prob = float(p[3])
                    pix_mean_dist = float(p[7])
                    N128_id = int(p[10])

                    pix_completeness = np.interp(pix_mean_dist,
                                                 completeness_values_dict[N128_id][0],
                                                 completeness_values_dict[N128_id][1])

                    renorm2dprob = pix_prob * (1.0 - pix_completeness)
                    pixel_completeness_records.append((pix_id, pix_completeness, renorm2dprob, -1.0, healpix_map_id))

                # Append to data to CSV
                try:
                    t1 = time.time()

                    print("Appending `%s`" % pixel_completeness_upload_csv)
                    with open(pixel_completeness_upload_csv, 'a') as csvfile:
                        csvwriter = csv.writer(csvfile)
                        for data in pixel_completeness_records:
                            csvwriter.writerow(data)

                    t2 = time.time()
                    print("\n********* start DEBUG ***********")
                    print("Pixel Completeness CSV append execution time: %s" % (t2 - t1))
                    print("********* end DEBUG ***********\n")
                except Error as e:
                    print("Error in creating Pixel Completeness CSV:\n")
                    print(e)
                    print("\nExiting")
                    return 1

                # batch_insert(pixel_update, pixel_update_records)

                # clean up
                print("freeing `completeness_values_dict`...")
                print("freeing `pixel_completeness_records`...")
                print("freeing `healpix_pixel_completeness_result`...")
                del completeness_values_dict
                del pixel_completeness_records
                del healpix_pixel_completeness_result

                t2 = time.time()
                print("\n********* start DEBUG ***********")
                print("SELECT %s/%s complete - execution time: %s" % (iter_num, number_of_iters, (t2 - t1)))
                print("********* end DEBUG ***********\n")

                iter_num += 1

                process = psutil.Process(os.getpid())
                mb = process.memory_info().rss / 1e+6
                print("\nTotal mem usage: %0.3f [MB]\n" % mb)

            print("Out of loop...")

            # FINISH WORK HERE
            t1 = time.time()

            print("\n%s:%s" % (i, k))
            batch = n128_pixel_id_list[i:k]
            batch_id_string = ','.join(map(str, batch))

            completeness_result = query_db([completeness_select % batch_id_string])[0]
            healpix_pixel_completeness_result = \
            query_db([healpix_pixel_completeness_select % (healpix_map_id, batch_id_string)])[0]

            completeness_values_dict = {}
            for c in completeness_result:
                N128_id = int(c[0])

                if N128_id not in completeness_values_dict:
                    completeness_values_dict[N128_id] = [[0.0], [1.0]]

                comp = float(c[2]) if float(c[2]) <= 1.0 else 1.0
                completeness_values_dict[N128_id][0].append(float(c[1]))
                completeness_values_dict[N128_id][1].append(comp)

            # Add point at infinity
            for key, val in completeness_values_dict.items():
                completeness_values_dict[key][0].append(np.inf)
                completeness_values_dict[key][1].append(0.0)

            # clean up
            print("freeing `completeness_result`...")
            del completeness_result

            pixel_completeness_records = []
            for p in healpix_pixel_completeness_result:
                pix_id = p[0]
                pix_prob = float(p[3])
                pix_mean_dist = float(p[7])
                N128_id = int(p[10])

                pix_completeness = np.interp(pix_mean_dist,
                                             completeness_values_dict[N128_id][0],
                                             completeness_values_dict[N128_id][1])

                renorm2dprob = pix_prob * (1.0 - pix_completeness)
                pixel_completeness_records.append((pix_id, pix_completeness, renorm2dprob, -1.0, healpix_map_id))

            t2 = time.time()
            print("\n********* start DEBUG ***********")
            print("SELECT %s/%s complete - execution time: %s" % (iter_num, number_of_iters, (t2 - t1)))
            print("********* end DEBUG ***********\n")

            # Append to data to CSV
            try:
                t1 = time.time()
                print("Appending `%s`" % pixel_completeness_upload_csv)
                with open(pixel_completeness_upload_csv, 'a') as csvfile:
                    csvwriter = csv.writer(csvfile)
                    for data in pixel_completeness_records:
                        csvwriter.writerow(data)

                t2 = time.time()
                print("\n********* start DEBUG ***********")
                print("Pixel Completeness CSV append execution time: %s" % (t2 - t1))
                print("********* end DEBUG ***********\n")
            except Error as e:
                print("Error in creating Pixel Completeness CSV:\n")
                print(e)
                print("\nExiting")
                return 1

            t1 = time.time()
            upload_sql = """LOAD DATA LOCAL INFILE '%s' 
                    INTO TABLE HealpixPixel_Completeness 
                    FIELDS TERMINATED BY ',' 
                    LINES TERMINATED BY '\n' 
                    (HealpixPixel_id, PixelCompleteness, Renorm2DProb, NetPixelProb, HealpixMap_id);"""

            success = bulk_upload(upload_sql % pixel_completeness_upload_csv)
            if not success:
                print("\nUnsuccessful bulk upload. Exiting...")
                return 1

            t2 = time.time()
            print("\n********* start DEBUG ***********")
            print("CSV upload execution time: %s" % (t2 - t1))

            try:
                print("Removing `%s`..." % pixel_completeness_upload_csv)
                os.remove(pixel_completeness_upload_csv)

                # clean up
                print("freeing `completeness_values_dict`...")
                print("freeing `pixel_completeness_records`...")
                print("freeing `healpix_pixel_completeness_result`...")
                del completeness_values_dict
                del pixel_completeness_records
                del healpix_pixel_completeness_result

                print("... Done")
            except Error as e:
                print("Error in file removal")
                print(e)
                print("\nExiting")
                return 1

        if not build_completeness_func:
            print("Alternate completeness func...")
            t1 = time.time()
            pixel_completeness_upload_csv = "%s/%s_pixel_completeness_upload.csv" % (
                formatted_healpix_dir, self.options.gw_id)
            composed_completeness_dict = None
            with open(pickle_output_dir + 'composed_completeness_dict.pkl', 'rb') as handle:
                composed_completeness_dict = pickle.load(handle)

            pixel_completeness_records = []
            for p_index, p_tuple in map_pixel_dict.items():
                pix_id = p_tuple[0]
                pix_prob = float(p_tuple[3])
                pix_mean_dist = float(p_tuple[7])
                N128_id = int(p_tuple[10])

                pix_completeness = np.interp(pix_mean_dist,
                                             np.asarray(composed_completeness_dict[N128_id][0]),
                                             np.asarray(composed_completeness_dict[N128_id][1]))

                renorm2dprob = pix_prob * (1.0 - pix_completeness)
                pixel_completeness_records.append((pix_id, pix_completeness, renorm2dprob, -1.0, healpix_map_id))

            t2 = time.time()
            print("\n********* start DEBUG ***********")
            print("Alternate Completeness Calc complete - execution time: %s" % (t2 - t1))
            print("********* end DEBUG ***********\n")

            # Append to data to CSV
            try:
                t1 = time.time()
                print("Appending `%s`" % pixel_completeness_upload_csv)
                with open(pixel_completeness_upload_csv, 'a') as csvfile:
                    csvwriter = csv.writer(csvfile)
                    for data in pixel_completeness_records:
                        csvwriter.writerow(data)

                t2 = time.time()
                print("\n********* start DEBUG ***********")
                print("Pixel Completeness CSV append execution time: %s" % (t2 - t1))
                print("********* end DEBUG ***********\n")
            except Error as e:
                print("Error in creating Pixel Completeness CSV:\n")
                print(e)
                print("\nExiting")
                return 1

            t1 = time.time()
            upload_sql = """LOAD DATA LOCAL INFILE '%s' 
                                INTO TABLE HealpixPixel_Completeness 
                                FIELDS TERMINATED BY ',' 
                                LINES TERMINATED BY '\n' 
                                (HealpixPixel_id, PixelCompleteness, Renorm2DProb, NetPixelProb, HealpixMap_id);"""

            success = bulk_upload(upload_sql % pixel_completeness_upload_csv)
            if not success:
                print("\nUnsuccessful bulk upload. Exiting...")
                return 1

            t2 = time.time()
            print("\n********* start DEBUG ***********")
            print("CSV upload execution time: %s" % (t2 - t1))

            try:
                print("Removing `%s`..." % pixel_completeness_upload_csv)
                os.remove(pixel_completeness_upload_csv)

                # clean up
                print("freeing `completeness_values_dict`...")
                print("freeing `pixel_completeness_records`...")
                print("freeing `healpix_pixel_completeness_result`...")
                # del completeness_values_dict
                del composed_completeness_dict
                del pixel_completeness_records
                # del healpix_pixel_completeness_result

                print("... Done")
            except Error as e:
                print("Error in file removal")
                print(e)
                print("\nExiting")
                return 1


        if not build_galaxy_weights:
            print("Skipping galaxy weights ...")
        else:
            print("Building galaxy weights ...")

            t1 = time.time()
            # Set & Retrieve NetProbToGalaxies
            # healpix_map_NetProbToGalaxies_update = '''
            #     UPDATE HealpixMap hm
            #     SET NetProbToGalaxies = (SELECT 1-SUM(hpc.Renorm2DProb)
            #                              FROM HealpixPixel_Completeness hpc
            #                              WHERE HealpixMap_id = %s
            #                             )
            #     WHERE hm.id = %s;
            # '''

            total_pixel_prob_select = '''
                SELECT SUM(Prob) FROM HealpixPixel WHERE HealpixMap_id = %s
            '''
            total_pixel_prob = float(query_db([total_pixel_prob_select % healpix_map_id])[0][0][0])


            healpix_map_NetProbToGalaxies_update = '''
                UPDATE HealpixMap hm
                SET NetProbToGalaxies = (SELECT %s-SUM(hpc.Renorm2DProb) 
                                         FROM HealpixPixel_Completeness hpc 
                                         WHERE HealpixMap_id = %s
                                        )
                WHERE hm.id = %s;
            '''

            query_db([healpix_map_NetProbToGalaxies_update % (total_pixel_prob, healpix_map_id, healpix_map_id)], commit=True)

            select_NetProbToGalaxies = "SELECT NetProbToGalaxies FROM HealpixMap WHERE id = %s"
            select_NetProbToGalaxies_result = query_db([select_NetProbToGalaxies % healpix_map_id])[0]
            net_prob_to_galaxies = select_NetProbToGalaxies_result[0][0]
            print("Net probability to galaxies: %0.5f" % net_prob_to_galaxies)

            t2 = time.time()
            print("\n********* start DEBUG ***********")
            print("healpix_map_NetProbToGalaxies execution time: %s" % (t2 - t1))
            print("********* end DEBUG ***********\n")

            tt1 = time.time()

            galaxy_4D_select = '''
            SELECT 
                g.id,
                g.z_dist,
                g.z_dist_err,
                g.Bweight,
                hp.id, 
                hp.Pixel_Index,
                hp.Prob,
                hp.Mean,
                hp.Stddev,
                ABS(g.z_dist - hp.Mean)/SQRT(POW(hp.Stddev, 2) + POW(g.z_dist_err, 2)) as SigmaTotal 
            FROM Galaxy g
            JOIN HealpixPixel hp on hp.Pixel_Index = g.%s_Pixel_Index
            WHERE hp.HealpixMap_id = %s;
            '''
            precompute_result = query_db([galaxy_4D_select % ("N%s" % map_nside, healpix_map_id)])[0]
            t2 = time.time()

            print("\n********* start DEBUG ***********")
            print("Precompute Select execution time: %s" % (t2 - t1))
            print("********* end DEBUG ***********\n")

            galaxy_attributes = []
            four_d_prob_norm = 0.0
            for r in precompute_result:
                g_id = int(r[0])
                hp_id = int(r[4])

                lum_prob = float(r[3])
                two_d_prob = float(r[6])
                sigma_total = float(r[9])
                z_prob = 1.0 - erf(sigma_total)

                four_d_prob = z_prob * two_d_prob * lum_prob
                four_d_prob_norm += four_d_prob

                galaxy_attributes.append([g_id, hp_id, lum_prob, z_prob, two_d_prob, four_d_prob])

            # clean up
            print("freeing `precompute_result`...")
            del precompute_result

            galaxy_attribute_data = []
            for g in galaxy_attributes:
                # norm_4d_weight = g[4] / four_d_prob_norm

                g_id = g[0]
                hp_id = g[1]
                lum_prob = g[2]
                z_prob = g[3]
                two_d_prob = g[4]
                four_d_prob = g[5]

                # norm_4d_weight = g[5] / four_d_prob_norm
                norm_4d_weight = four_d_prob / four_d_prob_norm
                net_gal_prob = norm_4d_weight * net_prob_to_galaxies

                galaxy_attribute_data.append(
                    # (g[0], g[1], g[2], g[3], g[4], norm_4d_weight, norm_4d_weight * net_prob_to_galaxies))
                    (g_id, hp_id, lum_prob, z_prob, two_d_prob, norm_4d_weight, net_gal_prob, healpix_map_id))

            # clean up
            print("freeing `galaxy_attributes`...")
            del galaxy_attributes

            # Create CSV, upload, and clean up CSV
            galaxy_attributes_upload_csv = "%s/%s_galaxy_attributes_upload.csv" % (
            formatted_healpix_dir, self.options.gw_id)
            try:
                t1 = time.time()
                print("Creating `%s`" % galaxy_attributes_upload_csv)
                with open(galaxy_attributes_upload_csv, 'w') as csvfile:
                    csvwriter = csv.writer(csvfile)
                    for data in galaxy_attribute_data:
                        csvwriter.writerow(data)

                t2 = time.time()
                print("\n********* start DEBUG ***********")
                print("CSV creation execution time: %s" % (t2 - t1))
                print("********* end DEBUG ***********\n")
            except Error as e:
                print("Error in creating CSV:\n")
                print(e)
                print("\nExiting")
                return 1

            t1 = time.time()
            upload_sql = """LOAD DATA LOCAL INFILE '%s' 
                    INTO TABLE HealpixPixel_Galaxy_Weight 
                    FIELDS TERMINATED BY ',' 
                    LINES TERMINATED BY '\n' 
                    (Galaxy_id, HealpixPixel_id, LumWeight, zWeight, Prob2DWeight, Norm4DWeight, GalaxyProb, HealpixMap_id);"""

            success = bulk_upload(upload_sql % galaxy_attributes_upload_csv)
            if not success:
                print("\nUnsuccessful bulk upload. Exiting...")
                return 1

            t2 = time.time()
            print("\n********* start DEBUG ***********")
            print("CSV upload execution time: %s" % (t2 - t1))

            try:
                print("Removing `%s`..." % galaxy_attributes_upload_csv)
                os.remove(galaxy_attributes_upload_csv)

                # Clean up
                print("freeing `galaxy_attribute_data`...")
                del galaxy_attribute_data

                print("... Done")
            except Error as e:
                print("Error in file removal")
                print(e)
                print("\nExiting")
                return 1

            tt2 = time.time()
            print("\n********* start DEBUG ***********")
            print("HealpixPixel_Galaxy Galaxy Attribute execution time: %s" % (tt2 - tt1))
            print("********* end DEBUG ***********\n")

            # update the healpix pixel net prob column...
            print("Updating healpix pixel net prob...")
            t1 = time.time()

            healpix_pixel_net_prob_update = '''
                UPDATE HealpixPixel_Completeness hpc1 
                JOIN 
                ( 
                    SELECT 
                        hpc2.id, 
                        (hpc2.Renorm2DProb + IFNULL(SUM(hp_g_w.GalaxyProb),0.0)) as NetPixelProb 
                    FROM HealpixPixel_Completeness hpc2 
                    LEFT JOIN HealpixPixel_Galaxy_Weight hp_g_w on hp_g_w.HealpixPixel_id = hpc2.HealpixPixel_id
                    WHERE hpc2.HealpixMap_id = %s
                    GROUP BY hpc2.id
                ) temp on hpc1.id = temp.id 
                SET hpc1.NetPixelProb = temp.NetPixelProb 
                WHERE hpc1.HealpixMap_id = %s;
            '''




            query_db([healpix_pixel_net_prob_update % (healpix_map_id, healpix_map_id)], commit=True)
            t2 = time.time()
            print("\n********* start DEBUG ***********")
            print("Update HealpixPixel NetProb execution time: %s" % (t2 - t1))
            print("********* end DEBUG ***********\n")


if __name__ == "__main__":
    useagestring = """python Generate_Tiles.py [options]

Example with healpix_dir defaulted to 'Events/<gwid>':
python LoadMap.py --gw_id <gwid> --healpix_file <filename>

Example with healpix_dir specified:
python LoadMap.py --gw_id <gwid> --healpix_dir Events/<directory name> --healpix_file <filename>
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
    print("Teglon `LoadMap` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")
