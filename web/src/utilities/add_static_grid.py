import argparse
import time

from web.src.services.teglon import *


if __name__ == "__main__":
    start = time.time()

    teglon = Teglon()

    parser = argparse.ArgumentParser()
    parser.add_argument('--teglon_detector_id', type=int, help='Teglon Detector ID')
    parser.add_argument('--detector_prefix', type=str, help='Prefix for static grid tile numbers')


    args = parser.parse_args()
    teglon.add_static_grid(**vars(args))

    end = time.time()
    duration = (end - start)
    print("\n********* start DEBUG ***********")
    print("Teglon `add_static_grid` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")