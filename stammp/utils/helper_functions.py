import os
import subprocess
import sys
import logging

logger = logging.getLogger()


def native_wordcount(file_path):
    if not os.path.isfile(file_path):
        raise ValueError('%r is not a path to an existing file' % file_path)

    proc = subprocess.Popen(['wc', '-l', file_path], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = proc.communicate()
    line_str, *_ = stdout.split()
    return int(line_str)


def prepare_output_dir(dir_path):
    if os.path.exists(dir_path) and not os.path.isdir(dir_path):
        raise ValueError('%r is not a path to a directory' % dir_path)

    if os.path.isdir(dir_path) and not os.access(dir_path, os.W_OK):
        raise ValueError('%r is not writable' % dir_path)

    if not os.path.exists(dir_path):
        try:
            os.makedirs(dir_path)
        except:
            raise ValueError('output directory %r cannot be created' % dir_path)


def execute(cmd, exit=True):
    if isinstance(cmd, list):
        cmd = ' '.join(cmd)
    logger.debug('executing command %r' % cmd)
    try:
        proc = subprocess.Popen(args=cmd, shell=True, stderr=subprocess.PIPE,
                                stdout=subprocess.PIPE, universal_newlines=True)
        retcode = proc.wait()
        stdout, stderr = proc.communicate()
        if retcode != 0:
            logger.error('\n\n%s', stderr.strip())
            raise Exception
        return stdout, stderr
    except:
        logger.error('Error while executing: %r', cmd)
        if exit:
            sys.exit(1)
