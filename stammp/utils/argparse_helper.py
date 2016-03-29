import argparse
import os


def dir_rwx(path):
    if not os.path.isdir(path):
        msg = '%r is not an existing directory' % path
        raise argparse.ArgumentTypeError(msg)
    if not os.access(path, os.R_OK | os.W_OK | os.X_OK):
        msg = 'no read-write-execute access on %r' % path
        raise argparse.ArgumentTypeError(msg)
    return path


def file_rw(path):
    if not os.path.isfile(path):
        dir_rwx(os.path.dirname(path))
    elif not os.access(path, os.R_OK | os.W_OK):
        msg = 'no read-write access on %r' % path
        raise argparse.ArgumentTypeError(msg)
    return path


def file_r(path):
    if not os.path.isfile(path):
        msg = '%r does not exist' % path
        raise argparse.ArgumentTypeError(msg)
    elif not os.access(path, os.R_OK):
        msg = 'no read access on %r' % path
        raise argparse.ArgumentTypeError(msg)
    return path


def file_rw_or_dir_rwx(path):
    if os.path.isdir(path):
        dir_rwx(path)
    else:
        file_rw(path)
    return path
