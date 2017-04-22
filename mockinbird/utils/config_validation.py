import glob
import re
import os
import logging
from mockinbird.utils.misc import deprecated


logger = logging.getLogger()


class Annot:
    """Annot provides detailed information on a configuration option

        It stores following properties:

            type (callable): initial conversion of raw value from the config file; \
                             mostly obsolete now

            default: default value if option is not explicitly set

            converter (callable): function to convert and validate config option. \
                                  Can raise ValueError if invalid value is passed

            warn_if_missing: print a warning message if option is not provided by the user
    """
    def __init__(self, type=lambda x: x, default=None, converter=lambda x: x,
                 warn_if_missing=False):
        """
        Constructs an annotation object describing a configuration option

        Args:
            type (callable): initial conversion of raw value from the config file; \
                             mostly obsolete now
            default: default value if option is not explicitly set
            converter (callable): function to convert and validate config option. \
                                  Can raise ValueError if invalid value is passed
            warn_if_missing: print a warning message if option is not provided by the user
        """
        self._type = type
        self._default = default
        self._converter = converter
        self._warn_if_missing = warn_if_missing

    @property
    def type(self):
        """Parse function to convert the 

            Note: starting from the introduction of the yaml config files, this
            function should be obsolete.
        """
        return self._type

    @property
    def default(self):
        """Default value if configuration is not provided by the user

        The default value ``None`` makes providing the config value mandatory.
        """
        return self._default

    @property
    def converter(self):
        """Function to convert and validate the given configuration value"""
        return self._converter

    @property
    def warn_if_missing(self):
        """Print a warning if the config option is not provided and falls back to the default"""
        return self._warn_if_missing


def rel_mapindex_validator(genome_index, cfg_path):
    """Validates a genome index

    The path can either be absolute or relative to the parent folder of ``cfg_path``
    """
    if not os.path.isabs(genome_index):
        parent_path = os.path.dirname(cfg_path)
        genome_index = os.path.join(os.path.abspath(parent_path), genome_index)
    genome_index_glob = "%s*" % genome_index
    if len(glob.glob(genome_index_glob)) == 0:
        raise ValueError('genome index %r does not exist' % genome_index)
    return genome_index


def rel_file_r_validator(path, cfg_path):
    """Validates a file add assures read permissions

    The path can either be absolute or relative to the parent folder of ``cfg_path``
    """
    if not os.path.isabs(path):
        parent_path = os.path.dirname(cfg_path)
        path = os.path.join(os.path.abspath(parent_path), path)
    return file_r_validator(path)


def rel_file_rw_validator(path, cfg_path):
    """Validates a file and assures read and write permissions

    The path can either be absolute or relative to the parent folder of ``cfg_path``
    """
    if not os.path.isabs(path):
        parent_path = os.path.dirname(cfg_path)
        path = os.path.join(os.path.abspath(parent_path), path)
    return file_rw_validator(path)


def rel_genome_validator(path, cfg_path):
    """Validates a genome and assures read permissions. Asserts the presence of a fasta index

    The path can either be absolute or relative to the parent folder of ``cfg_path``. The fasta
    index can be created by ``samtools faidx </path/to/file.fasta>``. The fasta index has to have
    the same name and end with ``.fai``.
    """
    if not os.path.isabs(path):
        parent_path = os.path.dirname(cfg_path)
        path = os.path.join(os.path.abspath(parent_path), path)
    path = file_r_validator(path)
    file_r_validator(path + '.fai')
    return path


def dnastr_validator(dna_string):
    """Validates that a string contains only the bases 'A', 'C', 'G' and 'T'

    Uracil ('U') letters are converted to the DNA equivalent Thymin ('T')
    """
    dna_string = dna_string.upper().replace('U', 'T')
    nuc_pat = re.compile('[ACTG]+')
    if not nuc_pat.match(dna_string):
        raise ValueError('DNA sequence %r is invalid. DNA sequences may only '
                         'contain A, C, G or T nucleotides.' % dna_string)
    return dna_string


def dnanuc_validator(dna_nuc):
    """Validates that the character is one of the four bases 'A', 'C', 'G' and 'T'"""
    if dna_nuc not in 'ACGT':
        raise ValueError('DNA nucleotide %r is invalid. Valid nucleotides: '
                         'A, C, G, T' % dna_nuc)
    return dna_nuc.upper()


def boolean(bool_str):
    """Converts a string to bool

    Only False, 'no', '0' and '' are interpreted as False, all other inputs are converted to True
    """
    if isinstance(bool_str, bool):
        return bool_str
    if bool_str.lower() in ('no', '0', ''):
        return False
    else:
        return True


def file_r_validator(path):
    """Validates a file and assures read permissions"""
    if not os.path.isfile(path):
        msg = '%r does not exist' % path
        raise ValueError(msg)
    elif not os.access(path, os.R_OK):
        msg = 'no read access on %r' % path
        raise ValueError(msg)
    return os.path.abspath(path)


def file_rw_validator(path):
    """Validates a file and assures read and write permissions"""
    if not os.path.isfile(path):
        parent_dir = os.path.dirname(path)
        if not os.access(parent_dir, os.W_OK | os.R_OK):
            msg = 'neither the file nor the parent directoy exists.'
            raise ValueError(msg)
    elif not os.access(path, os.W_OK | os.R_OK):
        msg = 'no read/write access on %r' % path
        raise ValueError(msg)
    return os.path.abspath(path)


def nonneg_integer(integer):
    """Validates that the input is a non-negative integer"""
    if integer < 0:
        raise ValueError('Non-negative integer expected. Got %s.' % integer)
    return integer


def in_set_validator(item, item_set):
    """Validates ``item`` is a member of the set ``item_set``"""
    if item not in item_set:
        raise ValueError('%r is not in set %r' % (item, item_set))
    return item


@deprecated('is_subset_validator is a relict from pre-yaml config files and should not be used')
def is_subset_validator(item_str, item_set):
    """Validates ``item_str`` is a subset of the set ``item_set``

    item_str is a comma separated list of items.
    """
    items = []
    for item in item_str.split(','):
        in_set_validator(item, item_set)
        items.append(item)
    return items


@deprecated('comma_set_args is a relict from pre-yaml config files and should not be used')
def comma_sep_args(item_str):
    """Split the input string after the comma delimiter into a list"""
    if item_str.strip() == '':
        return []
    else:
        return item_str.split(',')


def id_converter(x):
    """Return the input. Equivalent to ``lambda x: x``"""
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
                logger.error(ex)
                raise ConfigError(msg)
    return cfg_dict


@deprecated('mand_config is a relict from pre-yaml config files and should not be used')
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


@deprecated('opt_config is a relict from pre-yaml config files and should not be used')
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
