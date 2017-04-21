import logging

__all__ = ['obj']

LOG_DEFAULT_FORMAT = '%(asctime)s [%(levelname)s]  %(message)s'
LOG_LEVEL_MAP = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warn': logging.WARN,
    'error': logging.ERROR,
}

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
