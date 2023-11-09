import argparse
import time

from web.src.services.teglon import *


if __name__ == "__main__":
    start = time.time()

    teglon = Teglon()

    parser = argparse.ArgumentParser()
    parser.add_argument('--band_name', type=str, help='Name of the band')
    parser.add_argument('--lambda_eff', type=float, default=None,
                        help='Effective wavelength')
    parser.add_argument('--f99_coef', type=float, default=None, help='F99 coeficient to convert E(B-V) to MW extinction')


    args = parser.parse_args()
    teglon.add_band(**vars(args))

    end = time.time()
    duration = (end - start)
    print("\n********* start DEBUG ***********")
    print("Teglon `add_band` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")