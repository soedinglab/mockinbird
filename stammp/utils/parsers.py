from collections import namedtuple, OrderedDict
from urllib.parse import quote, unquote
import logging

logger= logging.getLogger()


GFFRecord = namedtuple('GFFRecord', [
    'seqid', 'source', 'type', 'start', 'end', 'score',
    'strand', 'phase', 'attributes',
])


def attr_type(attr_str):
    d = OrderedDict()
    try:
        for mapping in attr_str.split(';'):
            key, value = mapping.split('=')
            d[unquote(key)] = unquote(value)
    except ValueError:
        d['annot'] = attr_str
    return d


def gff_dict_mapper(attr_dict):
    items = []
    for key, value in attr_dict.items():
        items.append('%s=%s' % (quote(bytes(key, 'utf-8')), quote(bytes(value, 'utf-8'))))
    return ';'.join(items)


def str_or_none_mapper(string):
    if string is None:
        return '.'
    else:
        return quote(str(string))

GFF_FIELD_TYPES = [str, str, str, int, int, float, str, int, attr_type]
GFF_FIELD_MAPPERS = [str_or_none_mapper] * 6 + [str] + [str_or_none_mapper] + [gff_dict_mapper]


class GFF3Parser:

    def __init__(self, handle):
        self._handle = handle

    def parse(self):
        for line in self._handle:
            if line.startswith('#'):
                continue
            line = line.strip()
            tokens = line.split('\t')
            if len(tokens) != 9:
                continue
            values = []
            for value, field_type in zip(tokens, GFF_FIELD_TYPES):
                if value == '.':
                    values.append(None)
                else:
                    values.append(field_type(value))
            yield GFFRecord(*values)

    def write_record(self, record, **repl):
        items = []
        dict_repr = record._asdict()
        for key, value in repl.items():
            if key in dict_repr:
                dict_repr[key] = value

        vals = dict_repr.values()
        for value, mapper in zip(vals, GFF_FIELD_MAPPERS):
            items.append(mapper(value))
        print(*items, sep='\t', file=self._handle)


PCRecord = namedtuple('PCRecord', [
    'seqid', 'position', 'transitions', 'coverage', 'pvalue',
    'strand', 'occupancy',
])


PC_FIELD_TYPES = [str, int, int, int, float, str, float]


class PCTableParser:

    def __init__(self, handle):
        self._handle = handle
        header = handle.readline().split()
        mand_fields = [
            'seqid', 'position', 'transitions', 'coverage', 'score',
            'strand', 'occupancy'
        ]
        mand_conv = [str, int, int, int, float, str, float]

        if header[:7] != mand_fields:
            raise ValueError('data is not a PCTable')

        self._fields = header
        self._converters = dict(zip(mand_fields, mand_conv))

    def register_converter(self, field, converter):
        self._converters[field] = converter

    def parse(self):
        # we already read the header, first line is line 2
        fields = self._fields
        Record = namedtuple('PCRecord', fields)
        conv = list(map(self._converters, fields))
        for line_no, line in enumerate(self._handle, start=2):
            tokens = line.split()
            if len(tokens) != len(self._fields):
                logger.warn('PCTable: skipping line %s due to an unexpected number of fields',
                            line_no)
                continue
            values = []
            for value, field_type in zip(tokens, conv):
                if value == '.':
                    values.append(None)
                else:
                    values.append(field_type(value))
            yield Record(*values)
