import os, stat


def chmod_outputfile(output_file_path):
    os.chmod(output_file_path, mode=stat.S_IRWXU|stat.S_IRWXG|stat.S_IROTH|stat.S_IXOTH)