import time
import pickle
import os
from collections import OrderedDict

from dustmaps.config import config as dustmaps_config
from dustmaps.sfd import SFDQuery

from web.src.utilities.HEALPix_Helpers import *
from web.src.utilities.Database_Helpers import *

from web.src.objects.Detector import *
from web.src.objects.Tile import *
from web.src.objects.SQL_Polygon import *
from web.src.objects.Pixel_Element import *
from web.src.objects.Completeness_Objects import *

class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--build_skypixels_pickle', action="store_true", default=False,
                          help='Build sky pixel table; HACK, remove table FK on parent pixel id first, then restore '
                               'after update (Default = False)')

        parser.add_option('--build_ebv_pickle', action="store_true", default=False,
                          help='Build EBV table (Default = False)')

        parser.add_option('--build_composed_completeness_pickle', action="store_true", default=False,
                          help='Compose completness records into function input (Default = False)')

        return (parser)

    def main(self):

        utilities_base_dir = "./web/src/utilities"
        pickle_output_dir = "%s/pickles/" % utilities_base_dir
        if not os.path.exists(pickle_output_dir):
            os.makedirs(pickle_output_dir)

        # Set up dustmaps config
        # DC - Why is this directory set this way?
        dustmaps_config["data_dir"] = utilities_base_dir

        nside2 = 2
        nside2_npix = hp.nside2npix(nside2)
        frac_sky_nside2 = (1 / nside2_npix)

        nside4 = 4
        nside4_npix = hp.nside2npix(nside4)
        frac_sky_nside4 = (1 / nside4_npix)

        nside8 = 8
        nside8_npix = hp.nside2npix(nside8)
        frac_sky_nside8 = (1 / nside8_npix)

        nside16 = 16
        nside16_npix = hp.nside2npix(nside16)
        frac_sky_nside16 = (1 / nside16_npix)

        nside32 = 32
        nside32_npix = hp.nside2npix(nside32)
        frac_sky_nside32 = (1 / nside32_npix)

        nside64 = 64
        nside64_npix = hp.nside2npix(nside64)
        frac_sky_nside64 = (1 / nside64_npix)

        nside128 = 128
        nside128_npix = hp.nside2npix(nside128)
        frac_sky_nside128 = (1 / nside128_npix)

        nside_npix = [nside2_npix, nside4_npix, nside8_npix, nside16_npix, nside32_npix, nside64_npix, nside128_npix]
        nsides = [nside2, nside4, nside8, nside16, nside32, nside64, nside128]

        sky_pixels = None

        if self.options.build_skypixels_pickle:

            sky_pixel_select = "SELECT id, Pixel_Index FROM SkyPixel WHERE NSIDE = %s;"
            sky_pixel_select2 = "SELECT id, NSIDE, Pixel_Index FROM SkyPixel;"

            sky_pixels = OrderedDict()

            for i, nside in enumerate(nsides):
                print("Initializing NSIDE=%s pix..." % nside)
                sky_pixels[nside] = {}

                for j in range(nside_npix[i]):
                    pe = Pixel_Element(j, nside, 0.0)
                    pe.query_polygon_string # initialize
                    sky_pixels[nside][j] = pe

                if i > 0:
                    for j, (pi, pe) in enumerate(sky_pixels[nside].items()):
                        # Get Parent Pixel
                        theta, phi = hp.pix2ang(pe.nside, pe.index)
                        parent_index = hp.ang2pix(nsides[i-1], theta, phi)

                        if not parent_index in sky_pixels[nsides[i-1]]:
                            raise("Orphaned NSIDE=%s pixel! Orphaned index, theta, phi: (%s, %s, %s)" % (nside, pe.index, theta, phi))
                        else:
                            parent_pixel = sky_pixels[nsides[i-1]][parent_index]
                            sky_pixels[nside][j].parent_pixel = parent_pixel
                print("... created %s NSIDE=%s pix." % (len(sky_pixels[nside]), nside))

            for i, (nside, pixel_dict) in enumerate(sky_pixels.items()):

                if i == 0:
                    continue
                else: # Associate parents
                    parent_ids = query_db([sky_pixel_select % (nsides[i - 1])])

                    for row in parent_ids[0]:
                        sky_pixels[nsides[i-1]][int(row[1])].id = int(row[0]) # Set pixel id in parent...

            pixel_tuples = query_db([sky_pixel_select2])[0]
            for tup in pixel_tuples:
                db_id = int(tup[0])
                p_nside = int(tup[1])
                p_index = int(tup[2])

                sky_pixels[p_nside][p_index].id = db_id

            with open(pickle_output_dir + "sky_pixels.pkl", 'wb') as handle:
                pickle.dump(sky_pixels, handle, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            print("Skipping Sky Pixels...")
            print("\tLoading existing pixels...")
            if os.path.exists(pickle_output_dir + "sky_pixels.pkl"):
                with open(pickle_output_dir + "sky_pixels.pkl", 'rb') as handle:
                    sky_pixels = pickle.load(handle)
            else:
                print('sky_pixels.pkl does not exist!')

        if self.options.build_ebv_pickle:
            print("\n\nBuilding MWE...")

            t1 = time.time()

            theta, phi = hp.pix2ang(nside=nside128, ipix=np.arange(hp.nside2npix(nside128)))
            pix_coord = coord.SkyCoord(ra=np.rad2deg(phi), dec=np.rad2deg(0.5 * np.pi - theta), unit=(u.deg, u.deg))

            print("Retrieving dust info")
            sfd = SFDQuery()
            ebv = sfd(pix_coord)

            t2 = time.time()

            print("\n********* start DEBUG ***********")
            print("EBV query - execution time: %s" % (t2 - t1))
            print("********* end DEBUG ***********\n")

            print("E(B-V) data length: %s" % len(ebv))
            min_ebv = np.min(ebv)
            max_ebv = np.max(ebv)
            print("Max E(B-V): %s" % min_ebv)
            print("Min E(B-V): %s" % max_ebv)

            with open(pickle_output_dir + 'ebv.pkl', 'wb') as handle:
                pickle.dump(ebv, handle, protocol=pickle.HIGHEST_PROTOCOL)

        if self.options.build_composed_completeness_pickle:
            print("Composing Completeness...")
            composed_completeness_dict = OrderedDict()

            n128_ids = []
            for pix_index, pix_ele in sky_pixels[nside128].items():
                n128_ids.append(pix_ele.id)

            select_composed_completeness = '''
                SELECT id, MeanDistance, SmoothedCompleteness, N128_SkyPixel_id 
                FROM CompletenessGrid 
                ORDER BY N128_SkyPixel_id, MeanDistance;
            '''

            composed_completeness_result = query_db([select_composed_completeness])[0]
            for ccr in composed_completeness_result:
                mean_dist = float(ccr[1])
                smoothed_completeness = float(ccr[2])
                n128_pix_id = int(ccr[3])

                if n128_pix_id not in composed_completeness_dict:
                    composed_completeness_dict[n128_pix_id] = [[], []]

                composed_completeness_dict[n128_pix_id][0].append(mean_dist)
                composed_completeness_dict[n128_pix_id][1].append(smoothed_completeness)

            # append a point at infinity
            for key, val in composed_completeness_dict.items():
                composed_completeness_dict[key][0].append(np.inf)
                composed_completeness_dict[key][1].append(0.0)

            with open(pickle_output_dir + 'composed_completeness_dict.pkl', 'wb') as handle:
                pickle.dump(composed_completeness_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    useagestring = """python InitializeTeglon.py [options]

Example non-box command line:
python InitializeTeglon.py  --is_debug True --build_skydistances True --build_skypixels True --build_detectors True
                            --build_bands True --build_MWE True --plot_MWE True 
                            --build_galaxy_skypixel_associations True --build_completeness True 
                            --plot_completeness True --build_static_grids True

All args default to true, to only install/restore certain tables, set non-used args to False
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
    print("Teglon `InitializeTeglon` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")