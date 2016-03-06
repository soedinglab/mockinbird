import os
import subprocess


def native_wordcount(file_path):
    if not os.path.isfile(file_path):
        raise ValueError('%r is not a path to an existing file' % file_path)

    proc = subprocess.Popen(['wc', '-l', file_path], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = proc.communicate()
    line_str, *_ = stdout.split()
    return int(line_str)
