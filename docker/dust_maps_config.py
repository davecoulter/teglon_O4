from dustmaps.config import config

class Dustmaps_Config:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--data_dir', default="", type="str", help='Absolute path to data directory')

        return (parser)

    def main(self):
        config['data_dir'] = self.options.data_dir

if __name__ == "__main__":
    useagestring = """python dust_maps_config.py --data_dir { /path/to/your/data/dir }"""

    dmc = Dustmaps_Config()
    parser = dmc.add_options(usage=useagestring)
    options, args = parser.parse_args()
    dmc.options = options

    dmc.main()

    print("\n********* start DEBUG ***********")
    print("Dustmaps config set.")
    print("********* end DEBUG ***********\n")