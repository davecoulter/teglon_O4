import argparse
import time

from web.src.services.teglon import *


if __name__ == "__main__":
    start = time.time()

    teglon = Teglon()

    parser = argparse.ArgumentParser()
    parser.add_argument('--gw_id', default="", type=str,
                      help='LIGO superevent name, e.g. `S190425z` ')
    parser.add_argument('--healpix_dir', default='./web/events/{GWID}', type=str,
                      help='Directory for where to look for the healpix file.')
    parser.add_argument('--healpix_file', default="bayestar.fits.gz", type=str,
                      help='Healpix filename. Default `bayestar.fits.gz`.')
    parser.add_argument('--tele', default="a", type=str,
                      help='''Telescope abbreviation for the telescope to extract files for. Default: `a`.
                      Available: (s:Swope, t:Thacher, n:Nickel, t80: T80S_T80S-Cam, a:All). If `All` selected,
                      combination plot is generated (vs single plot).''')
    parser.add_argument('--band', default="r", type=str,
                      help='''Only used when plotting single telescope, otherwise defaults per telescope are used.
                           For single telescope, default: `r`. Available: (g, r, i, z, I)''')
    parser.add_argument('--extinct', default=0.5, type=float,
                      help='''Extinction in mags in the specified band to be less than.
                      Default: 0.5 mag. Must be > 0.0''')
    parser.add_argument('--tile_file', default="{FILENAME}", type=str,
                      help='''Only used when plotting single telescope. Filename of tiles to plot. Otherwise, will
                      default to tile files of the form: `{GWID}_{TELE_NAME}_4D_0.9_bayestar.fits.gz.txt`.''')
    parser.add_argument('--num_tiles', default=-1, type=int,
                      help='''Number of tiles to plot for each tile file. If < 0, plot all. Otherwise, plot indices
                      [0:num_tiles].''')
    parser.add_argument('--cum_prob_outer', default=0.9, type=float,
                      help='''Cumulative prob to cover in tiles. Default 0.9. Must be > 0.2 and < 0.95
                      and > `cum_prob_inner`''')
    parser.add_argument('--cum_prob_inner', default=0.5, type=float,
                      help='''Cumulative prob to cover in tiles. Default 0.5. Must be > 0.2 and < 0.95
                      and < `cum_prob_inner`''')

    args = parser.parse_args()
    teglon.plot_teglon(**vars(args))

    end = time.time()
    duration = (end - start)
    print("\n********* start DEBUG ***********")
    print("Teglon `plot_teglon` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")


