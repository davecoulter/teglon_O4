import os

import matplotlib
matplotlib.use("Agg")
from matplotlib import colors
from matplotlib import pyplot as plt

from datetime import datetime, timezone, timedelta

from astropy.time import Time
from astropy.coordinates import get_sun
from dustmaps.config import config as dustmaps_config
from dustmaps.sfd import SFDQuery
from ligo.skymap.postprocess.util import find_greedy_credible_levels
import ligo.skymap.plot
from astroplan import Observer
from astropy.coordinates import AltAz, EarthLocation
from astropy import units as u
import astropy.coordinates as coord
from astropy.table import Table

from web.src.objects.Detector import *
from web.src.objects.Tile import *
from web.src.objects.Pixel_Element import *
from web.src.utilities.Database_Helpers import *
from web.src.utilities.filesystem_utilties import *


class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--gw_id', default="", type="str",
                          help='LIGO superevent name, e.g. `S190425z` ')

        parser.add_option('--healpix_dir', default='./web/events/{GWID}', type="str",
                          help='Directory for where to look for the healpix file.')

        parser.add_option('--healpix_file', default="bayestar.fits.gz", type="str",
                          help='Healpix filename. Default `bayestar.fits.gz`.')

        parser.add_option('--tele', default="a", type="str",
                          help='''Telescope abbreviation for the telescope to extract files for. Default: `a`. 
                          Available: (s:Swope, t:Thacher, n:Nickel, t80: T80S_T80S-Cam, a:All). If `All` selected, 
                          combination plot is generated (vs single plot).''')

        parser.add_option('--band', default="r", type="str",
                          help='''Only used when plotting single telescope, otherwise defaults per telescope are used.
                               For single telescope, default: `r`. Available: (g, r, i, z, I)''')

        parser.add_option('--extinct', default="0.5", type="float",
                          help='''Extinction in mags in the specified band to be less than. 
                          Default: 0.5 mag. Must be > 0.0''')

        parser.add_option('--tile_file', default="{FILENAME}", type="str",
                          help='''Only used when plotting single telescope. Filename of tiles to plot. Otherwise, will 
                          default to tile files of the form: `{GWID}_{TELE_NAME}_4D_0.9_bayestar.fits.gz.txt`.''')

        parser.add_option('--num_tiles', default="-1", type="int",
                          help='''Number of tiles to plot for each tile file. If < 0, plot all. Otherwise, plot indices
                          [0:num_tiles].''')

        parser.add_option('--cum_prob_outer', default="0.9", type="float",
                          help='''Cumulative prob to cover in tiles. Default 0.9. Must be > 0.2 and < 0.95 
                          and > `cum_prob_inner`''')

        parser.add_option('--cum_prob_inner', default="0.5", type="float",
                          help='''Cumulative prob to cover in tiles. Default 0.5. Must be > 0.2 and < 0.95 
                          and < `cum_prob_inner`''')

        return (parser)

    def main(self):

        is_error = False
        nside128 = 128

        # Set up external resource directories
        utilities_base_dir = "./web/src/utilities"
        pickle_output_dir = "%s/pickles/" % utilities_base_dir

        print("Resetting dustmaps config...")
        dustmaps_config.reset()
        dustmaps_config["data_dir"] = utilities_base_dir

        detector_mapping = {
            "s": "SWOPE",
            "t": "THACHER",
            "n": "NICKEL",
            "t80": "T80S_T80S-Cam",
            "a": "All"
        }

        detector_geography = {
            "SWOPE": EarthLocation(lat=-29.0182 * u.deg, lon=-70.6926 * u.deg, height=2402 * u.m),
            "THACHER": EarthLocation(lat=34.46479 * u.deg, lon=-121.6431 * u.deg, height=630.0 * u.m),
            "NICKEL": EarthLocation(lat=37.3414 * u.deg, lon=-121.6429 * u.deg, height=1283 * u.m),
            "T80S_T80S-Cam": EarthLocation(lat=-70.8035 * u.deg, lon=-70.8035 * u.deg, height=2207.0 * u.m),
        }

        band_mapping = {
            "g": "SDSS g",
            "r": "SDSS r",
            "i": "SDSS i",
            "z": "SDSS z",
            "I": "Landolt I"
        }

        # Parameter checks
        print("\n\n**************************\nChecking Parameters...")
        if self.options.gw_id == "":
            is_error = True
            print("GWID is required.")

        if self.options.band not in band_mapping:
            is_error = True
            print("Invalid band selection. Available bands: %s" % band_mapping.keys())

        if self.options.tele not in detector_mapping:
            is_error = True
            print("Invalid telescope selection. Available telescopes: %s" % detector_mapping.keys())

        if is_error:
            print("\n\nErrors! See above output.")
            print("Exiting...")
            print("\n**************************")
            return 1
        else:
            print("Checked!\n**************************")


        formatted_healpix_dir = self.options.healpix_dir
        formatted_healpix_dir = formatted_healpix_dir.replace("{GWID}", self.options.gw_id)

        plot_all = True
        if self.options.tele != "a":
            plot_all = False

        tile_files = {}
        if plot_all:
            # get each default tile file for each instrument.
            default_tile_file_formatted = "{GWID}_{TELE_NAME}_4D_0.9_bayestar.fits.gz.txt"

            for detect_key, detect_val in detector_mapping.items():
                if detect_key != "a":
                    tile_files[detect_val] = "%s/%s" % (formatted_healpix_dir, default_tile_file_formatted.format(
                        GWID=self.options.gw_id, TELE_NAME=detect_val
                    ))
        else:
            tile_files[detector_mapping[self.options.tele]] = "%s/%s" % (formatted_healpix_dir, self.options.tile_file)


        # Create canvas
        plot_axes = {}
        ax_list = []
        if plot_all:
            fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2,
                                                                     subplot_kw={'projection': 'astro hours mollweide'})
            fig.delaxes(ax6)
            plt.tight_layout()

            ax_list += [ax1, ax2, ax3, ax4, ax5]

            for i, (detector_name, tile_file) in enumerate(tile_files.items()):
                 _ax = ax_list[i]
                 _ax.grid()
                 _ax.tick_params(axis='both', labelsize=5)
                 plot_axes[detector_name] = _ax

            ax_list[-1].grid()
            ax_list[-1].tick_params(axis='both', labelsize=5)
        else:
            fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'astro hours mollweide'})
            ax.grid()
            plot_axes[detector_mapping[self.options.tele]] = ax


        # Generate all plottable assets
        healpix_map_select = "SELECT id, RescaledNSIDE, t_0 FROM HealpixMap WHERE GWID = '%s' and Filename = '%s'"
        healpix_map_result = query_db([healpix_map_select % (self.options.gw_id, self.options.healpix_file)])[0][0]
        healpix_map_id = int(healpix_map_result[0])
        healpix_map_event_time = Time(healpix_map_result[2], format="gps")
        healpix_map_nside = int(healpix_map_result[1])

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

        map_2d_pix_dict = {int(mpr[2]): float(mpr[3]) for mpr in map_pix_result}
        map_4d_pix_dict = {int(mpr[2]): float(mpr[4]) for mpr in map_pix_result}

        num_pix = hp.nside2npix(healpix_map_nside)
        map_2d_pix = np.zeros(num_pix)
        map_4d_pix = np.zeros(num_pix)

        for i, mp in enumerate(map_2d_pix): # same number of pix, re-use enumeration
            if i in map_2d_pix_dict:
                map_2d_pix[i] = map_2d_pix_dict[i]
            if i in map_4d_pix_dict:
                map_4d_pix[i] = map_4d_pix_dict[i]

        _90_50_levels_2d = find_greedy_credible_levels(np.asarray(map_2d_pix))
        _90_50_levels_4d = find_greedy_credible_levels(np.asarray(map_4d_pix))

        # Plot Probability
        for detector_name, plot_ax in plot_axes.items():
            plot_ax.title.set_text(detector_name)
            plot_ax.contourf_hpx(_90_50_levels_2d, cmap='OrRd', levels=[0.0, 0.5, 0.9], linewidths=0.2, alpha=0.3)

        # if plotting All telescopes
        if plot_all:
            ax_list[-1].title.set_text("Resampled 4D Probability")
            ax_list[-1].contourf_hpx(_90_50_levels_4d, cmap='OrRd_r', levels=[0.0, 0.5, 0.9], linewidths=0.2, alpha=0.75)

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

        for detector_name, plot_ax in plot_axes.items():
            plot_ax.contour_hpx(ebv, colors='dodgerblue', levels=[0.5, np.max(ebv)], linewidths=0.2, alpha=0.5)
            plot_ax.contourf_hpx(ebv, cmap='Blues', levels=[0.5, np.max(ebv)], linewidths=0.2, alpha=0.3)

        # Sun
        time_of_trigger = Time(healpix_map_event_time.to_datetime(), scale="utc")
        theta, phi = hp.pix2ang(nside=64, ipix=np.arange(hp.nside2npix(64)))
        ra = np.rad2deg(phi)
        dec = np.rad2deg(0.5 * np.pi - theta)
        all_sky_coords = coord.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))

        # Pointings
        detector_select_by_name = '''
            SELECT
                id, 
                Name,
                Deg_width,
                Deg_height,
                Deg_radius,
                Area,
                MinDec,
                MaxDec,
                ST_AsText(Poly)
            FROM Detector 
            WHERE Name='%s'
        '''
        for detector_name, tile_file in tile_files.items():

            plot_ax = plot_axes[detector_name]

            detector_result = query_db([detector_select_by_name % detector_name])[0][0]
            detector_id = int(detector_result[0])
            detector_name = detector_result[1]
            detector_poly = detector_result[8]
            detector_vertices = Detector.get_detector_vertices_from_teglon_db(detector_poly)

            detector_obj = Detector(detector_name, detector_vertices, detector_id=detector_id)
            tiles = Table.read(tile_file, format='ascii.ecsv')
            tiles_to_plot = [Tile(ra, dec, detector=detector_obj, nside=healpix_map_nside,
                                  net_prob=prob) for ra, dec, prob in zip(tiles["RA"], tiles["Dec"], tiles["Prob"])]

            min_prob = np.min(tiles["Prob"])
            max_prob = np.max(tiles["Prob"])
            print("min prob for `%s`: %s" % (detector_name, min_prob))
            print("max prob for `%s`: %s" % (detector_name, max_prob))
            clr_norm = colors.LogNorm(min_prob, max_prob)
            for i, t in enumerate(tiles_to_plot[0:self.options.num_tiles]):
                if t.ra_deg < 358 and t.ra_deg > 2:
                    t.plot2(plot_ax, edgecolor='k',
                            facecolor=plt.cm.Greens(clr_norm(t.net_prob)),
                            linewidth=0.05, alpha=1.0, zorder=9999)

            # Sun Contours
            observatory = Observer(location=detector_geography[detector_name])
            sun_set = observatory.sun_set_time(time_of_trigger, which='nearest', horizon=-18 * u.deg)
            sun_rise = observatory.sun_rise_time(sun_set, which='next', horizon=-18 * u.deg)
            sun_coord = get_sun(sun_set)
            hours_of_the_night = int((sun_rise - sun_set).to_value('hr'))
            time_next = sun_set

            sun_separations = sun_coord.separation(all_sky_coords)
            deg_sep = np.asarray([sp.deg for sp in sun_separations])

            plot_axes[detector_name].contour_hpx(deg_sep, colors='darkorange', levels=[0, 60.0], alpha=0.5, linewidths=0.2)
            plot_axes[detector_name].contourf_hpx(deg_sep, colors='gold', levels=[0, 60.0], linewidths=0.2, alpha=0.3)
            plot_axes[detector_name].plot(sun_coord.ra.degree, sun_coord.dec.degree,
                                       transform=plot_axes[detector_name].get_transform('world'), marker="o",
                      markersize=8, markeredgecolor="darkorange", markerfacecolor="gold", linewidth=0.05)

            # Airmass constraints
            net_airmass = all_sky_coords.transform_to(AltAz(obstime=sun_set, location=detector_geography[detector_name])).secz
            for i in np.arange(hours_of_the_night):
                time_next += timedelta(1/24.)
                temp = all_sky_coords.transform_to(AltAz(obstime=time_next, location=detector_geography[detector_name])).secz
                temp_index = np.where((temp >= 1.0) & (temp <= 2.0))
                net_airmass[temp_index] = 1.5

            plot_axes[detector_name].contourf_hpx(net_airmass, colors='gray', levels=[np.min(net_airmass), 1.0], linewidths=0.2, alpha=0.3)
            plot_axes[detector_name].contourf_hpx(net_airmass, colors='gray', levels=[2.0, np.max(net_airmass)], linewidths=0.2, alpha=0.3)

        output_file = "all_telescopes_4D_0.9_bayestar.fits.gz.svg"
        if not plot_all:
            output_file = self.options.tile_file.replace(".txt", ".svg")

        output_path = "%s/%s" % (formatted_healpix_dir, output_file)
        fig.savefig(output_path, bbox_inches='tight', format="svg")  # ,dpi=840
        chmod_outputfile(output_path)
        plt.close('all')
        print("... Done.")




if __name__ == "__main__":
    useagestring = """python Plot_Teglon.py [options]

python Plot_Teglon.py --gw_id <gwid> --healpix_file <filename> --extinct 0.50 --prob_type 2D --cum_prob_outer 0.90 
--cum_prob_inner 0.50 --band r --tile_file <FILENAME> --tele s --galaxies_file <FILENAME> 
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
    print("Teglon `Plot_Teglon` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")


