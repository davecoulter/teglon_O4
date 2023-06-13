import os
import time
from shutil import copyfile


class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--model_dir', default='./web/models/grb', type="str",
                          help='Working directory.')

        parser.add_option('--sub_dir', default="onaxis", type="str",
                          help='GRB Model sub directory (for batching)')

        parser.add_option('--bin_size', default="50", type="int",
                          help='Number of models to put in each bin')

        return (parser)

    def main(self):

        is_error = False

        grb_model_path = "%s/%s" % (self.options.model_dir, self.options.sub_dir)
        grb_model_files = []
        for file_index, file in enumerate(os.listdir(grb_model_path)):
            if file.endswith(".dat"):
                grb_model_files.append("%s/%s" % (grb_model_path, file))

        if len(grb_model_files) <= 0:
            is_error = True
            print("There are no models to process.")

        if is_error:
            print("Exiting...")
            return 1

        j = 0
        for i, model_file in enumerate(grb_model_files):
            current_sub_dir = "%s/%s" % (grb_model_path, j)
            if i % self.options.bin_size == 0:
                j += 1
                current_sub_dir = "%s/%s" % (grb_model_path, j)
                os.mkdir(current_sub_dir)

            base_name = os.path.basename(model_file)
            copyfile(model_file, current_sub_dir + "/%s" % base_name)
            print("`%s` moved. Deleting..." % model_file)
            os.remove(model_file)

        print("Done.")


if __name__ == "__main__":
    useagestring = """python ComputeModelDetection.py [options]

Example with healpix_dir defaulted to 'Events/<gwid>' amd model_dir to 'Events/{GWID}/Models':
python ComputeModelDetection.py --gw_id <gwid> --healpix_file <filename>

Assumes .dat files with Astropy ECSV format: 
https://docs.astropy.org/en/stable/api/astropy.io.ascii.Ecsv.html

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
    print("Teglon `GroupGRBModels` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")
