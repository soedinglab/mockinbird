from collections import namedtuple, OrderedDict
from urllib.parse import quote, unquote

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

    def parse(self):
        self._handle.readline()
        for line in self._handle:
            tokens = line.split()
            if len(tokens) != 7:
                continue
            values = []
            for value, field_type in zip(tokens, PC_FIELD_TYPES):
                if value == '.':
                    values.append(None)
                else:
                    values.append(field_type(value))
            yield PCRecord(*values)
