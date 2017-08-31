import argparse
import os

from mockinbird import __version__


def dir_rwx(path):
    if not os.path.isdir(path):
        msg = '%r is not an existing directory' % path
        raise argparse.ArgumentTypeError(msg)
    if not os.access(path, os.R_OK | os.W_OK | os.X_OK):
        msg = 'no read-write-execute access on %r' % path
        raise argparse.ArgumentTypeError(msg)
    return path


def dir_rwx_create(path):
    if not os.path.isdir(path):
        try:
            os.makedirs(path)
        except OSError:
            msg = '%r does not exist and cannot be created' % path
            raise argparse.ArgumentTypeError(msg)
    if not os.access(path, os.R_OK | os.W_OK | os.X_OK):
        msg = 'no read-write-execute access on %r' % path
        raise argparse.ArgumentTypeError(msg)
    return path


def dir_rx(path):
    if not os.path.isdir(path):
        msg = '%r is not an existing directory' % path
        raise argparse.ArgumentTypeError(msg)
    if not os.access(path, os.R_OK | os.X_OK):
        msg = 'no read-execute access on %r' % path
        raise argparse.ArgumentTypeError(msg)
    return path


def file_rw(path):
    path = os.path.abspath(path)
    if not os.path.isfile(path):
        dir_rwx(os.path.dirname(path))
    elif not os.access(path, os.R_OK | os.W_OK):
        msg = 'no read-write access on %r' % path
        raise argparse.ArgumentTypeError(msg)
    return path


def file_r(path):
    if os.path.isdir(path):
        msg = '%r is a directory - expected a file' % path
        raise argparse.ArgumentTypeError(msg)
    elif not os.path.isfile(path):
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


def add_version_arguments(parser):
    parser.add_argument('--version', '-v', action='version',
                        version='%(prog)s {version}'.format(version=__version__))
