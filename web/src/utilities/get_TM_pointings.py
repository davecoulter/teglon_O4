# import argparse
import time
from astropy.time import Time
import urllib.request
import requests
import json
import os
from configparser import RawConfigParser
import pickle
import csv

class TeglonTemp:

    def get_skycoord_from_tm_point(self, tm_point_string):
        ra_dec_list = tm_point_string.split("(")[1].split(")")[0].split(" ")

        ra = float(ra_dec_list[0])
        dec = float(ra_dec_list[1])

        return ra, dec

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--grace_id', type="str", help='''Grace ID used in TM.''')

        return (parser)

    def main(self):

        detector_map = {
            47: "gw190426_ztf_tiles.txt",
            22: "gw190426_mmt_tiles.txt",
            12: "gw190426_swift_tiles.txt",
            71: "gw190426_goto_tiles.txt"
        }

        pointings_by_instrument_id = None
        do_TM_download = False

        if do_TM_download:
            tm_api_token = None
            tm_base_url = None

            # Sanity
            is_error = False
            try:
                configFile = './Settings.ini'
                config = RawConfigParser()
                config.read(configFile)

                tm_api_token = config.get('treasuremap', 'TM_API_TOKEN')
                tm_base_url = config.get('treasuremap', 'TM_ENDPOINT')
            except:
                print("Error! No Treasure Map keys (`TM_API_TOKEN`, `TM_ENDPOINT`) in the Settings.ini file!")
                is_error = True

            if not any(tm_api_token) or not any(tm_base_url):
                print("Can't load TM detectors with API keys or API endpoints!")
                is_error = True

            pointing_url = 'pointings'
            pointing_request = {
                "api_token": tm_api_token,
                "status": "completed",
                # Omitting bands includes all bands...
                # "bands": [
                #     "U", "B", "V", "g", "r", "i"
                # ],
                "graceid": self.options.grace_id
            }
            pointing_target_url = "{}/{}?{}".format(tm_base_url, pointing_url,
                                                    urllib.parse.urlencode(pointing_request))
            pointing_response = requests.get(url=pointing_target_url)
            pointing_list = json.loads(pointing_response.text)

            pointings_by_instrument_id = {}

            xrt_id = 13 # we're only handling UVOIR instruments
            for p in pointing_list:
                tm_detect_id = p["instrumentid"]

                if tm_detect_id == xrt_id:
                    continue

                if tm_detect_id not in pointings_by_instrument_id:
                    pointings_by_instrument_id[tm_detect_id] = []

                # source field ra dec mjd filter exp_time mag_lim p.a.
                field = p["id"]
                ra, dec = self.get_skycoord_from_tm_point(p["position"])
                obs_time = Time(p["time"], scale="utc")
                mjd = obs_time.mjd
                filter = p["band"]
                exp_time = 0.0
                mag_lim = p["depth"]
                pa = p["pos_angle"]

                pointings_by_instrument_id[tm_detect_id].append((
                    "treasure_map",
                    field,
                    ra,
                    dec,
                    mjd,
                    filter,
                    exp_time,
                    mag_lim,
                    pa
                ))

            with open('./web/src/utilities/tm_pointings.pkl', 'wb') as handle:
                pickle.dump(pointings_by_instrument_id, handle, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open('./web/src/utilities/tm_pointings.pkl', 'rb') as handle:
                pointings_by_instrument_id = pickle.load(handle)

        for key, val in pointings_by_instrument_id.items():
            file_name = detector_map[key]

            with open("./web/events/S190425z/observed_tiles/%s" % file_name, 'w') as csvfile:
                csvwriter = csv.writer(csvfile, delimiter=' ')
                csvwriter.writerow(('source', 'field', 'ra', 'dec', 'mjd', 'filter', 'exp_time', 'mag_lim', 'p.a.'))

                for row in val:
                    csvwriter.writerow(row)



if __name__ == "__main__":
    useagestring = """python get_TM_pointings.py [options]

    python get_TM_pointings.py --grace_id <Grace ID>
    """
    start = time.time()

    # teglon = Teglon()
    #
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--grace_id', type=str, help='Grace ID used in TM')
    #
    # args = parser.parse_args()
    # teglon.main(**vars(args))

    teglon = TeglonTemp()
    parser = teglon.add_options(usage=useagestring)
    options, args = parser.parse_args()
    teglon.options = options

    teglon.main()

    end = time.time()
    duration = (end - start)
    print("\n********* start DEBUG ***********")
    print("Teglon `get_TM_pointings` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")
