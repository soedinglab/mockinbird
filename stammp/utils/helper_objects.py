from collections import namedtuple, defaultdict
import os

import numpy as np

from pyivtree import IVTree, GenomicInterval
from stammp.obj.functions import makeReverseComplement as rev_complement
from stammp.obj.gff import GFF
from stammp.obj import parclipsites
from stammp.utils.parsers import PCTableParser
# from stammp.utils.misc import deprecated


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


class PCContainer:
    def __init__(self):
        self._sites = defaultdict(IVTree)

    @classmethod
    def read_handle(cls, handle):
        container = cls()
        parser = PCTableParser(handle)
        for rec in parser.parse():
            site = GenomicInterval(rec.position, rec.position)
            site._data = rec
            container._sites[rec.seqid].insert(site)
        return container

    def occ_profile(self, seq_id, start, end, strand):
        occ_arr = np.zeros(end - start + 1)
        iv = GenomicInterval(start, end)
        for site in self._sites[seq_id].query_all_overlaps(iv):
            data = site._data
            if data.strand == strand:
                occ_arr[data.position - start] += data.occupancy
        return occ_arr


class ParclipSiteContainer:

    def __init__(self):
        self._pcs = parclipsites.ParclipSites()

    def loadFromFile(self, filename):
        self._pcs.loadFromFile(filename)

    def save2File(self, filename):
        self._pcs.save2File(filename)

    def print(self, start, stop):
        self._pcs.print(start, stop)

    def head(self):
        self._pcs.head()

    def tail(self):
        self._pc.tail()

#    @deprecated('please use the len() function instead')
    def size(self):
        return self._pcs.size()

    def sort(self, key='occ'):
        self._pcs.sort(key)

    def removeSites(self, index):
        self._pcs.removeSite(index)

    def addSite(self, chrname, position, m, r, pvalue, strand, occupancy):
        self._pcs.addSite(chrname, position, m, r, pvalue, strand, occupancy)

    def removeSitesLocatedInGFF(self, gff, width=0):
        self._pcs.removeSitesLocatedInGFF(gff, width)

    def save2Fasta(self, genome, filename, start, stop, width=15):
        with open(filename, 'w') as out_handle:
            for i, seq in self.get_sequences_it(genome, start, stop, width):
                if seq is None:
                    continue
                pc = self._pcs
                header = ':'.join(str(x) for x in ['>seq', i, pc.chrs[i], pc.pos[i], pc.strand[i],
                                  round(pc.occ[i], 2), pc.m[i], round(pc.m[i] / pc.occ[i], 2)])
                print(header, file=out_handle)
                print(seq, file=out_handle)

    def getSequences(self, genome, start, stop, width):
        return [seq for i, seq in self.get_sequences_it(genome, start, stop, width)]

    def get_sequences_it(self, genome, start, stop, width):
        pc_count = len(self._pcs)
        start_max = min(start, pc_count)
        end_max = min(stop, pc_count)

        for i in range(start_max, end_max):
            chrom = self._pcs.chrs[i]
            start = self._pcs.pos[i] - width
            end = self._pcs.pos[i] + width
            strand = self._pcs.strand[i]
            try:
                yield i, genome.get_sequence(chrom, start, end, strand)
            except ValueError:
                yield i, None

    def getChromosomePositions(self):
        self._pcs.getChromsomePositions()

    def exactSearch(self, chrname, position, strand, width=0):
        self._pcs.exactSearch(chrname, position, strand, width)

    def getValues(self, chrname, position, strand, sense, upstream, downstream):
        return self._pcs.getValues(chrname, position, strand, sense, upstream, downstream)

    def __len__(self):
        return len(self._pcs)

    def __getattr__(self, attr):
        return self._pcs.__getattribute__(attr)
