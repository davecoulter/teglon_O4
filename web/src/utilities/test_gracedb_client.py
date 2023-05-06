import time
import os

from ligo.gracedb.rest import GraceDb
import healpy as hp
from astropy.time import Time

import warnings
warnings.filterwarnings("ignore")

class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--superevent_id', default="", type="str",
                          help='LIGO superevent name, e.g. `S190425z`')

        parser.add_option('--healpix_dir', default='./web/events/{GWID}', type="str",
                          help='Directory for where to look for the healpix file (Default = ../Events/{GWID})')

        parser.add_option('--healpix_file', default="bayestar.fits.gz", type="str", help='healpix filename')

        return (parser)

    def main(self):

        is_error = False
        # Parameter checks
        if self.options.superevent_id == "":
            is_error = True
            print("`superevent_id` is required.")

        if is_error:
            print("Exiting...")
            return 1

        api_endpoint = "https://gracedb.ligo.org/api/"
        formatted_healpix_dir = self.options.healpix_dir
        if "{GWID}" in formatted_healpix_dir:
            formatted_healpix_dir = formatted_healpix_dir.replace("{GWID}", self.options.superevent_id)
        hpx_path = "%s/%s" % (formatted_healpix_dir, self.options.healpix_file)

        gdb_client = GraceDb(api_endpoint)
        event_json = gdb_client.superevent(self.options.superevent_id).json()
        far = event_json["far"]
        t_0 = event_json["t_0"]
        t_0_time_obj = Time(t_0, scale="utc", format="gps")
        t_0_datetime_str = t_0_time_obj.to_value(format="iso")
        url = event_json["links"]["self"]

        # build directories
        if not os.path.exists(formatted_healpix_dir):
            os.mkdir(formatted_healpix_dir)

        # Download and save map file.
        # The default self.options.healpix_file - "bayestar.fits.gz" always has the most up-to-date flat file
        #   in a search context
        file_response = gdb_client.files(self.options.superevent_id, self.options.healpix_file)
        with open(hpx_path, "wb") as f:
            f.write(file_response.data)

        # Load healpix file into memory to get header info
        (prob, distmu, distsigma, distnorm), header_gen = \
            hp.read_map(hpx_path, field=(0, 1, 2, 3), h=True)

        hdr = dict(header_gen)
        nside = hdr["NSIDE"]

        print(self.options.superevent_id, far, t_0_datetime_str, url, nside)

        test = 1


if __name__ == "__main__":
    useagestring = """python test_gracedb_client.py [options]

python test_gracedb_client.py --gw_id <gwid>

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
    print("Teglon `test_gracedb_client` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")