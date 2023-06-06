import argparse
import time

from web.src.services.teglon import *


if __name__ == "__main__":
    start = time.time()

    teglon = Teglon()

    parser = argparse.ArgumentParser()
    parser.add_argument('--tm_detector_id', type=int, help='Treasure Map detector ID')
    parser.add_argument('--min_dec', type=float, default=-90.0, help='Min dec in deg. Default = -90.0')
    parser.add_argument('--max_dec', type=float, default=90.0, help='Max dec in deg. Default = +90.0')
    parser.add_argument('--detector_geometry', type=str, default="rectangle", help='''Specify shape of detector. Valid 
    options: `rectangle` or `circle`. If geometry is an arbitrary polygon, leave blank.''')
    parser.add_argument('--detector_width', type=float, help='''Only used if `detector_geometry`==`rectangle`; width in 
    deg. Ignored if `detector_geometry` is not `rectangle`.''')
    parser.add_argument('--detector_height', type=float, help='''Only used if `detector_geometry`==`rectangle`; height in 
    deg. Ignored if `detector_geometry` is not `rectangle`.''')
    parser.add_argument('--detector_radius', type=float, help='''Only used if `detector_geometry`==`circle`; radius in 
            deg. Ignored if `detector_geometry` is not `circle`.''')

    args = parser.parse_args()
    teglon.add_detector(**vars(args))

    end = time.time()
    duration = (end - start)
    print("\n********* start DEBUG ***********")
    print("Teglon `add_detector` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")