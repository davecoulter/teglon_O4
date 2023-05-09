import argparse
import time

from web.src.services.teglon import *


if __name__ == "__main__":
    start = time.time()

    teglon = Teglon()

    parser = argparse.ArgumentParser()
    parser.add_argument('--gw_id', default="", type=str, help='LIGO superevent name, e.g. `S190425z`')
    parser.add_argument('--healpix_dir', default="./web/events/{GWID}", type=str,
                        help='Directory for where to look for the healpix file.')
    parser.add_argument('--healpix_file', default="bayestar.fits.gz", type=str,
                        help='Healpix filename. Default=`bayestar.fits.gz`')
    parser.add_argument('--tele', default="a", type=str,
                      help='''Telescope abbreviation for the telescope to extract files for. Default: `a`.
                              Available: (s:Swope, t:Thacher, n:Nickel, t80: T80S_T80S-Cam, a:All). If `All`
                              selected,combination plot is generated (vs single plot).''')
    parser.add_argument('--band', default="r", type=str,
                      help='Only used if `tele` specified. Default: `r`. Available: (g, r, i, z, I)')
    parser.add_argument('--extinct', default=0.5, type=float,
                      help='''Extinction in mags in the specified band(s) to be less than. Default: 0.5 mag.
                      Must be > 0.0''')
    parser.add_argument('--prob_type', default="4D", type=str,
                      help='''Probability type to consider. Default `4D`. Available: (4D, 2D)''')
    parser.add_argument('--cum_prob', default=0.9, type=float,
                      help='Cumulative prob to cover. Default 0.9. Must be > 0.2 and < 0.95')
    parser.add_argument('--num_tiles', default="1000", type=int,
                      help='''Return the top `num_tiles` tiles per telescope. Must be > 1. Default 1000.''')
    parser.add_argument('--min_ra', default=-1.0, type=float,
                      help='''Optional, decimal format. If specified, `min_ra`, `max_ra`, `min_dec`,
                  and `max_dec` are all required. `min_ra` must be >= 0.0 and < `max_ra`. Within the defined box,
                  tiles will be extracted within a given telescope's declination limits.''')
    parser.add_argument('--max_ra', default=-1.0, type=float,
                      help='''Optional, decimal format. If specified, `min_ra`, `max_ra`, `min_dec`,
                      and `max_dec` are all required. `max_ra` must be <= 360.0 and > `min_ra`. Within the defined
                      box, tiles will be extracted within a given telescope's declination limits.''')
    parser.add_argument('--max_dec', default=-1.0, type=float,
                      help='''Optional, decimal format. Must be If specified, `min_ra`, `max_ra`, `min_dec`,
                      and `max_dec` are all required. `max_dec` must be <= 90.0 and > `min_dec`. Within the defined
                      box, tiles will be extracted within a given  telescope's declination limits.''')
    parser.add_argument('--min_dec', default=-1.0, type=float,
                      help='''Optional, decimal format. Must be If specified, `min_ra`, `max_ra`, `min_dec`,
                      and `max_dec` are all required. `min_dec` must be >= -90.0 and < `max_dec`. Within the defined
                      box, tiles will be extracted within a given telescope's declination limits.''')

    args = parser.parse_args()
    teglon.extract_tiles(**vars(args))

    end = time.time()
    duration = (end - start)
    print("\n********* start DEBUG ***********")
    print("Teglon `extract_tiles` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")
