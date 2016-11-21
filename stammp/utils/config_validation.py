import re
from collections import namedtuple
import os
import logging
from operator import attrgetter


logger = logging.getLogger()

class Annot:
    def __init__(self, type, default=None, converter=lambda x:x,
                 warn_if_missing=True):
        self._type = type
        self._default = default
        self._converter = converter
        self._warn_if_missing = warn_if_missing

    @property
    def type(self):
        return self._type

    @property
    def default(self):
        return self._default

    @property
    def converter(self):
        return self._converter

    @property
    def warn_if_missing(self):
        return self._warn_if_missing


def rel_file_r_validator(path, cfg_path):
    if not os.path.isabs(path):
        parent_path = os.path.dirname(cfg_path)
        path = os.path.join(os.path.abspath(parent_path), path)
    return file_r_validator(path)


def rel_genome_validator(path, cfg_path):
    if not os.path.isabs(path):
        parent_path = os.path.dirname(cfg_path)
        path = os.path.join(os.path.abspath(parent_path), path)
    path = file_r_validator(path)
    file_r_validator(path + '.fai')
    return path


def dnastr_validator(dna_string):
    dna_string = dna_string.upper().replace('U', 'T')
    nuc_pat = re.compile('[ACTG]+')
    if not nuc_pat.match(dna_string):
        raise ValueError('DNA sequence %r is invalid. DNA sequences may only '
                         'contain A, C, G or T nucleotides.' % dna_string)
    return dna_string


def dnanuc_validator(dna_nuc):
    if dna_nuc not in 'ACGT':
        raise ValueError('DNA nucleotide %r is invalid. Valid nucleotides: '
                         'A, C, G, T' % dna_nuc)
    return dna_nuc.upper()


def boolean(bool_str):
    if isinstance(bool_str, bool):
        return bool_str
    if bool_str.lower() in ('no', '0', ''):
        return False
    else:
        return True


def file_r_validator(path):
    if not os.path.isfile(path):
        msg = '%r does not exists' % path
        raise ValueError(msg)
    elif not os.access(path, os.R_OK):
        msg = 'no read access on %r' % path
        raise ValueError(msg)
    return os.path.abspath(path)


def nonneg_integer(integer):
    if integer < 0:
        raise ValueError('Non-negative integer expected. Got %s.' % integer)
    return integer


def in_set_validator(item, item_set):
    if item not in item_set:
        raise ValueError('%r is not in set %r' % (item, item_set))
    return item


def is_subset_validator(item_str, item_set):
    items = []
    for item in item_str.split(','):
        in_set_validator(item, item_set)
        items.append(item)
    return items


def comma_sep_args(item_str):
    if item_str.strip() == '':
        return []
    else:
        return item_str.split(',')


def id_converter(x):
    return x


def validate_section(config, cfg_format):
    cfg_dict = {}
    for key, annot in cfg_format.items():
        if key not in config:
            if annot.default is not None:
                log_cmd = logger.warn if annot.warn_if_missing else logger.debug
                log_cmd('key %r undefined. Using default value %r.', key, annot.default)
                cfg_dict[key] = annot.default
            else:
                msg = ('mandatory configuration key %r missing' % key)
                logger.error(msg)
                raise ConfigError(msg)
        else:
            try:
                raw_data = annot.type(config[key])
                cfg_dict[key] = annot.converter(raw_data)
                logger.debug('set %r to %r', key, raw_data)
            except ValueError as ex:
                msg = 'invalid value %r for key %r' % (config[key], key)
                logger.error(msg)
                raise ConfigError(msg)
    return cfg_dict


def mand_config(config, cfg_format):
    cfg_dict = {}
    for section, fields_dict in cfg_format.items():
        if section not in config:
            logger.warn('config section %r missing.', section)
            config[section] = {}
        cfg_sec = config[section]
        sec_dict = {}
        for key, annot in fields_dict.items():
            if key not in cfg_sec:
                if annot.default is not None:
                    logger.warn('key %r missing in section %r. '
                                'Using default value %r.', key, section, annot.default)
                    sec_dict[key] = annot.default
                else:
                    msg = ('mandatory configuration key %r '
                           'missing in section %r.' % (key, section))
                    logger.error(msg)
                    raise ConfigError(msg)
            else:
                try:
                    print(cfg_sec)
                    raw_data = annot.type(cfg_sec[key])
                    sec_dict[key] = annot.converter(raw_data)
                    logger.debug('section %r: set %r to %r', section,
                                 key, raw_data)
                except ValueError as ex:
                    msg = 'invalid value for key %r in section %r' % (key, section)
                    logger.error(msg)
                    logger.error(str(ex))
                    raise ConfigError(msg)
        cfg_dict[section] = sec_dict
    return cfg_dict


def opt_config(config, cfg_format, id_field):
    cfg_list = []
    for section in config.sections():
        cfg_sec = config[section]
        if id_field not in cfg_sec:
            logger.debug('skipping section %r as it does not contain the %r key',
                         section, id_field)
            continue
        if 'skip' in cfg_sec:
            if cfg_sec.getboolean('skip'):
                logger.info('skipping section %r', section)
                continue
        type_id = cfg_sec[id_field]
        if type_id not in cfg_format:
            msg = 'unexpected value %r for key %r' % (type_id, id_field)
            logger.error(msg)
            raise ConfigError(msg)

        cfg_item = {}
        cfg_item[id_field] = type_id
        cfg_item['sec_name'] = section
        for key, annot in cfg_format[type_id].items():
            if key not in cfg_sec:
                if annot.default is not None:
                    logger.warn('key %r missing in section %r. '
                                'Using default value %r.', key, section, annot.default)
                    cfg_item[key] = annot.default
                else:
                    msg = ('mandatory configuration key %r '
                           'missing in section %r.' % (key, section))
                    logger.error(msg)
                    raise ConfigError(msg)
            else:
                try:
                    cfg_getter = type_getter[annot.type]
                    raw_data = cfg_getter(cfg_sec)(key)
                    cfg_item[key] = annot.converter(raw_data)
                    logger.debug('section %r: set %r to %r', section,
                                 key, raw_data)
                except ValueError as ex:
                    msg = 'invalid value for key %r in section %r' % (key, section)
                    logger.error(msg)
                    logger.error(str(ex))
                    raise ConfigError(msg)
        cfg_list.append(cfg_item)
    return cfg_list


class ConfigError(ValueError):
    pass
