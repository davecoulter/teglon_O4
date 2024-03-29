import os
import time
from shutil import copyfile


class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--KN_dir', default="./web/models/kne", type="str",
                          help='KN Model directory (to batch)')

        parser.add_option('--bin_size', default="50", type="int",
                          help='Number of models to put in each bin')

        parser.add_option('--start_dir', default="0", type="int",
                          help='Zero-indexed start directory name of the first dir to start enumerating from.')

        return (parser)

    def main(self):

        is_error = False

        kn_model_path = self.options.KN_dir
        kn_model_files = []
        for file_index, file in enumerate(os.listdir(kn_model_path)):
            if file.endswith(".dat"):
                kn_model_files.append("%s/%s" % (kn_model_path, file))

        if len(kn_model_files) <= 0:
            is_error = True
            print("There are no models to process.")

        if is_error:
            print("Exiting...")
            return 1

        j = self.options.start_dir
        for i, model_file in enumerate(kn_model_files):
            current_sub_dir = "%s/%s" % (kn_model_path, j)
            if i % self.options.bin_size == 0:
                j += 1
                current_sub_dir = "%s/%s" % (kn_model_path, j)
                os.mkdir(current_sub_dir)

            base_name = os.path.basename(model_file)
            copyfile(model_file, current_sub_dir + "/%s" % base_name)
            print("`%s` moved. Deleting..." % model_file)
            os.remove(model_file)

        print("Done.")


if __name__ == "__main__":
    useagestring = """python group_kne_files.py --KN_dir {dir relative to this dir}
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
    print("Teglon `group_kne_files` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")