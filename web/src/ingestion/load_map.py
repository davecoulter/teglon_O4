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
                      help='Directory for where to look for the healpix file (Default = ../Events/{GWID})')
    parser.add_argument('--healpix_file', default="bayestar.fits.gz", type=str, help='healpix filename')
    parser.add_argument('--analysis_mode', action="store_true", default=False,
                      help='''Upload the healpix file at the native resolution (Default = False), with all 
                        pixels. Otherwise, upload a rescaled version, and only the 90th percentile pixels.''')
    parser.add_argument('--clobber', action="store_true", default=False,
                      help='''Replace the existing (superevent, healpix_file) combination if it already exists. 
                      If `clobber` not specified and combination exists, return error message.''')
    parser.add_argument('--skip_swope', action="store_true", default=False,
                      help='''Do not register Swope tiles for this map (Default = False)''')
    parser.add_argument('--skip_thacher', action="store_true", default=False,
                      help='''Do not register Thacher tiles for this map (Default = False)''')
    parser.add_argument('--skip_t80', action="store_true", default=False,
                      help='''Do not register T80S_T80S-Cam tiles for this map (Default = False)''')
    parser.add_argument('--skip_newfirm', action="store_true", default=False,
                        help='''Do not register NEWFIRM tiles for this map (Default = False)''')

    args = parser.parse_args()
    teglon.load_map(**vars(args))

    end = time.time()
    duration = (end - start)
    print("\n********* start DEBUG ***********")
    print("Teglon `load_map` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")
