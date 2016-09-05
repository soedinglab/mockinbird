import logging

__version__ = '1.2.4'

__all__ = ['obj']

LOG_DEFAULT_FORMAT = '%(asctime)s [%(levelname)s]  %(message)s'
LOG_LEVEL_MAP = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warn': logging.WARN,
    'error': logging.ERROR,
}
