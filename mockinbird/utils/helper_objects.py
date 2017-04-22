from collections import namedtuple, defaultdict
import bisect
import os
import copy

import numpy as np

from mockinbird.ivtree import IVTree
from mockinbird.ivtree import GenomicInterval as Interval
from mockinbird.obj.functions import makeReverseComplement as rev_complement
from mockinbird.obj.gff import GFF
from mockinbird.utils.parsers import PCTableParser, GFF3Parser
from mockinbird.utils import parsers
from mockinbird.utils.misc import deprecated


SeqInfo = namedtuple('SeqInfo', ['bp_length', 'byte_offset', 'bp_linewidth', 'byte_linewidth'])


def parse_fasta_index(fasta_index):
    seq_info_dict = {}
    with open(fasta_index) as ind:
        for line in ind:
            try:
                seq_id, length, start_byte, bp_linewidth, byte_linewidth = line.split()
                length = int(length)
                bp_linewidth = int(bp_linewidth)
                byte_linewidth = int(byte_linewidth)
                start_byte = int(start_byte)
                seq_info = SeqInfo(length, start_byte, bp_linewidth, byte_linewidth)
                seq_info_dict[seq_id] = seq_info
            except ValueError:
                raise CorruptedFastaIndexException()
    return seq_info_dict


class CorruptedFastaIndexException(ValueError):
    pass


class NotInitializedException(ValueError):
    pass


class EfficientGenome:

    def __init__(self, fasta_file):
        self._fasta_file = fasta_file

        if not os.path.exists(fasta_file):
            raise IOError('fasta file %r does not exist' % fasta_file)
        fasta_index = fasta_file + '.fai'
        if not os.path.exists(fasta_index):
            raise IOError('fasta index %r does not exist' % fasta_index)
        self._seq_info = parse_fasta_index(fasta_index)

    def __enter__(self):
        self._fa_handle = open(self._fasta_file)
        self._initialized = True
        return self

    def __exit__(self, error, err_type, traceback):
        self._initialized = False
        self._fa_handle.close()
        return False

    def get_sequence(self, seqid, start, end, strand='+'):
        if not self._initialized:
            raise NotInitializedException()
        try:
            length, byte_offset, lw_bp, lw_bytes = self._seq_info[seqid]
        except KeyError:
            raise ValueError('unknown sequence identifier: %r' % seqid)

        if start <= 0:
            raise ValueError('start position must be greater than 0')
        if end > length:
            raise ValueError('end index %s greater than the sequence length %s' % (end, length))

        start -= 1
        stride = lw_bytes - lw_bp

        start_byte = byte_offset
        line_no = start // lw_bp
        start_byte += line_no * lw_bytes
        start_byte += start % lw_bp

        seqs = []
        read_bp = 0
        self._fa_handle.seek(start_byte)
        while read_bp < end - start:
            line_bp = lw_bp - ((start + read_bp) % lw_bp)
            bp = min(end - start - read_bp, line_bp)
            seq = self._fa_handle.read(bp)
            self._fa_handle.read(stride)
            read_bp += bp
            seqs.append(seq)
        nt_seq = ''.join(seqs).upper()
        if strand == '-':
            nt_seq = rev_complement(nt_seq)
        return nt_seq

Annotation = namedtuple('Annotation', ['chrom', 'start', 'end', 'strand', 'into', 'type'])


class GFFAnnotation:

    def __init__(self, gff_filepath):
        gff_object = GFF(gff_filepath)
        self._gff_obj = gff_object

        chr_index_mapping = defaultdict(list)
        for i, chrom in enumerate(gff_object.chr):
            chr_index_mapping[chrom].append(i)
        self._chrom_ind_map = chr_index_mapping

    def __len__(self):
        return self._gff_obj.size()

    def draw_random_annotation(self, seqid=None):
        if seqid is None:
            annot_index = np.random.choice(self._gff_obj.size())
        else:
            if seqid not in self._chrom_ind_map:
                return None
            annot_index = np.random.choice(self._chrom_ind_map[seqid])
        annot = Annotation(
            self._gff_obj.chr[annot_index],
            self._gff_obj.start[annot_index],
            self._gff_obj.stop[annot_index],
            self._gff_obj.strand[annot_index],
            self._gff_obj.info[annot_index],
            self._gff_obj.typeof[annot_index]
        )
        return annot


@deprecated('use ParclipSiteContainer instead')
class PCContainer:
    def __init__(self):
        self._sites = defaultdict(IVTree)

    @classmethod
    def read_handle(cls, handle):
        container = cls()
        parser = PCTableParser(handle)
        for rec in parser.parse():
            site = Interval(rec.position, rec.position)
            site._data = rec
            container._sites[rec.seqid].insert(site)
        return container

    def occ_profile(self, seq_id, start, end, strand):
        occ_arr = np.zeros(end - start + 1)
        iv = Interval(start, end)
        for site in self._sites[seq_id].query_all_overlaps(iv):
            data = site._data
            if data.strand == strand:
                occ_arr[data.position - start] += data.occupancy
        return occ_arr


class ParclipSiteContainer:

    def __init__(self, records=[]):
        data_list = []
        data_tree = defaultdict(IVTree)
        for record in records:
            iv = Interval(record.position, record.position)
            iv.data = record
            data_list.append(iv)
            data_tree[record.seqid].insert(iv)
        self._data_list = data_list
        self._data_tree = data_tree
        self._sort_keys = None
        self._descending = False

    @classmethod
    def from_file(cls, file_path):
        with open(file_path) as handle:
            parser = PCTableParser(handle)
            records = []
            for rec in parser.parse():
                records.append(rec)
            return cls(records)

    def get_fields(self):
        if len(self._data_list) > 0:
            return self._data_list[0].data._fields
        else:
            return parsers.PC_MANDATORY_FIELDS

    def save2File(self, file_path):
        with open(file_path, 'w') as handle:
            header = self.get_fields()
            print(*header, sep='\t', file=handle)
            if self._descending:
                iter_list = self._data_list[::-1]
            else:
                iter_list = self._data_list
            for iv in iter_list:
                print(*iv.data, sep='\t', file=handle)

    def save2Fasta(self, genome, filename, width=15):

        if not self._descending:
            data_list = self._data_list
        else:
            data_list = self._data_list[::-1]

        with open(filename, 'w') as out_handle:
            for i, iv in enumerate(data_list):
                data = iv.data
                chrom = data.seqid
                start = data.position - width
                end = data.position + width
                strand = data.strand

                try:
                    seq = genome.get_sequence(chrom, start, end, strand)
                    data = [
                        '>seq', i, chrom, data.position, strand, round(data.occupancy, 2),
                        data.transitions, round(data.transitions / data.occupancy, 2)
                    ]
                    header = ':'.join(str(x) for x in data)
                    print(header, file=out_handle)
                    print(seq, file=out_handle)
                except ValueError:
                    continue

    def sort(self, by=None, ascending=True):
        if by is None:
            sort_key = lambda r: (r.data.seqid, r.data.position)
        else:
            sort_key = lambda r: r.data.__getattribute__(by)
        self._data_list.sort(key=sort_key)
        self._sort_keys = list(map(sort_key, self._data_list))
        self._descending = not ascending
        self._key_fun = sort_key

    def __getitem__(self, indexer):
        con = copy.copy(self)
        data_tree = defaultdict(IVTree)
        con._data_tree = data_tree
        if not self._descending:
            data = self._data_list
            keys = self._sort_keys
        else:
            data = self._data_list[::-1]
            keys = self._sort_keys[::-1]

        new_data = data[indexer]
        if keys is not None:
            new_keys = keys[indexer]
        else:
            new_keys = None
        for iv in new_data:
            con._data_tree[iv.data.seqid].insert(iv)

        if self._descending:
            new_data.reverse()
            if new_keys is not None:
                new_keys.reverse()
        con._data_list = new_data
        con._sort_keys = new_keys

        return con

    def get_all_sequences(self, genome, width):

        if not self._descending:
            data_list = self._data_list
        else:
            data_list = self._data_list[::-1]

        seqs = []
        for iv in data_list:
            rec = iv.data
            chrom = rec.seqid
            start = rec.position - width
            end = rec.position + width
            strand = rec.strand

            try:
                seq = genome.get_sequence(chrom, start, end, strand)
                seqs.append(seq)
            except ValueError:
                continue
        return seqs

    def remove_gff_sites(self, gff_file, extend_bp=0):
        with open(gff_file) as gff_handle:
            parser = GFF3Parser(gff_handle)
            for rec in parser.parse():
                self.remove_ovl_sites(rec.seqid, rec.start - extend_bp,
                                      rec.end + extend_bp, rec.strand)

    def remove_ovl_sites(self, seqid, start, end, strand):
        ovl_iv = Interval(start, end)
        ivs = list(self._data_tree[seqid].query_all_overlaps(ovl_iv))
        for iv in ivs:
            data = iv.data
            if data.strand != strand:
                continue
            self._data_tree[data.seqid].remove(iv)
            if self._sort_keys is None:
                for i, iv_hit in enumerate(self._data_list):
                    if iv_hit is iv:
                        del_ind = i
                        break
                del self._data_list[del_ind]
            else:
                search_key = self._key_fun(iv)
                left_i = bisect.bisect_left(self._sort_keys, search_key)
                if self._data_list[left_i] is iv:
                    del self._data_list[left_i]
                    del self._sort_keys[left_i]
                else:
                    right_i = bisect.bisect_right(self._sort_keys, search_key)
                    del_ind = None
                    for i in range(left_i, right_i):
                        if self._data_list[i] is iv:
                            del_ind = i
                            break
                    assert del_ind is not None
                    del self._data_list[del_ind]
                    del self._sort_keys[del_ind]

    def get_overlaps(self, seqid, start, end, strand):
        ovl_iv = Interval(start, end)
        ivs = list(self._data_tree[seqid].query_all_overlaps(ovl_iv))
        for iv in ivs:
            if iv.data.strand != strand:
                continue
            yield iv.data

    def get_occ_profile(self, seqid, start, end, strand):
        occ_arr = np.zeros(end - start + 1)
        iv = Interval(start, end)
        for site in self._data_tree[seqid].query_all_overlaps(iv):
            data = site.data
            if data.strand == strand:
                occ_arr[data.position - start] += data.occupancy
        return occ_arr

    # legacy method for supporting the HeatMap plots. Should be depecated and replaced by
    # get_occ_profile.
    def getValues(self, seqid, start, strand, sense, upstream, downstream):
        if strand == '+':
            extr_start = start - upstream
            extr_end = start + downstream
        else:
            extr_start = start - downstream
            extr_end = start + upstream
        if sense:
            extr_strand = strand
        else:
            extr_strand = '+' if strand == '-' else '+'

        profile = self.get_occ_profile(seqid, extr_start, extr_end, extr_strand)
        if strand == '-':
            return profile[::-1]
        else:
            return profile

    def __iter__(self):
        if not self._descending:
            self._current = -1
        else:
            self._current = len(self._data_list)
        return self

    def __next__(self):
        if self._descending:
            self._current -= 1
            if self._current >= 0:
                return self._data_list[self._current].data
            else:
                raise StopIteration
        else:
            self._current += 1
            if self._current < len(self._data_list):
                return self._data_list[self._current].data
            else:
                raise StopIteration

    def __len__(self):
        return len(self._data_list)
