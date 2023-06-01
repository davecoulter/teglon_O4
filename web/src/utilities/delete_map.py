from web.src.utilities.Database_Helpers import *

class Teglon:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--healpix_map_id', default="-1", type="int",
                          help='''The integer primary key for the map to remove.''')

        return (parser)

    def main(self):
        is_error = False

        if self.options.healpix_map_id < 0:
            is_error = True
            print("HealpixMap_id is required.")

        create_baks = 'CALL BackupTables(%s);' % self.options.healpix_map_id
        lock_tables = '''
            LOCK TABLE HealpixMap WRITE, HealpixPixel WRITE, HealpixPixel_Completeness WRITE, HealpixPixel_Galaxy_Weight WRITE, 
            ObservedTile WRITE, ObservedTile_HealpixPixel WRITE, StaticTile_HealpixPixel WRITE, HealpixMap_bak WRITE, HealpixPixel_bak WRITE, 
            HealpixPixel_Completeness_bak WRITE, HealpixPixel_Galaxy_Weight_bak WRITE, ObservedTile_bak WRITE, ObservedTile_HealpixPixel_bak WRITE, 
            StaticTile_HealpixPixel_bak WRITE;
        '''
        delete_map = 'CALL DeleteMap(%s);' % self.options.healpix_map_id
        unlock_tables = 'UNLOCK TABLES;'

        query_db([create_baks], commit=True)
        query_db([lock_tables], commit=True)
        query_db([delete_map], commit=True)
        query_db([unlock_tables], commit=True)


if __name__ == "__main__":
    useagestring = """python DeleteMap.py [options]
    
python DeleteMap.py --healpix_map_id <database id>
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
    print("Teglon `DeleteMap` execution time: %s" % duration)
    print("********* end DEBUG ***********\n")


