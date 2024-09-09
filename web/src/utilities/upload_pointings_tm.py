import argparse
import time

from web.src.services.teglon import *


if __name__ == "__main__":
    start = time.time()

    teglon = Teglon()

    parser = argparse.ArgumentParser()
    parser.add_argument('--tm_detector_id', type=int, help='Treasure Map detector ID')
    parser.add_argument('--teglon_detector_id', type=int, help='Teglon-specific detector ID')
    parser.add_argument('--healpix_map_id', type=int, help='''Teglon-specific healpix map ID associated 
        with the pointings to upload.''')

    args = parser.parse_args()
    teglon.upload_pointings_tm(**vars(args))

    end = time.time()
    duration = (end - start)
    print("\n********* start DEBUG ***********")
    print("Teglon `upload_pointings_tm` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")