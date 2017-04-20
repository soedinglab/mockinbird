import unittest
import tempfile
import shutil
import os
import random

from stammp.utils.helper_objects import EfficientGenome
from stammp.utils.helper_objects import parse_fasta_index

fasta = """
>chr1
ACGTACGTACGT
AAAAACCCCCCA
TTTCCGGGCACT
TTATGGCGAGAG
GGGGGTTTTGGG
TTTGTATGTATA
TTTTT
>chr2
AAAAAAAAAAAA
CGTGATATGGGG
TTTAAAAAAATT
GGCCCCCGGAAG
GGGTATTTATAA
AAAA
""".strip()

index = """
chr1    77      6       12      13
chr2    64      96      12      13
""".strip()


class GenomeTests(unittest.TestCase):

    def test_foo(self):
        with EfficientGenome(GenomeTests._fasta_file) as genome:
            chrs = list(GenomeTests._seq_dict.keys())
            for i in range(10000):
                seq = random.choice(chrs)
                length, byte_offset, lw_bp, lw_bytes = GenomeTests._seq_info[seq]
                bp_len = random.randint(1, length)
                start = random.randint(0, length - bp_len)
                expected_seq = GenomeTests._seq_dict[seq][start:start + bp_len]
                observed_seq = genome.get_sequence(seq, start + 1, start + bp_len)
                self.assertEqual(expected_seq, observed_seq)

    @classmethod
    def setUpClass(cls):
        tmp_dir = tempfile.mkdtemp()

        fasta_file = os.path.join(tmp_dir, 'test.fa')
        with open(fasta_file, 'w') as handle:
            handle.write(fasta)
        cls._fasta_file = fasta_file

        index_file = os.path.join(tmp_dir, 'test.fa.fai')
        with open(index_file, 'w') as handle:
            handle.write(index)
        cls._tmp_dir = tmp_dir
        cls._seq_info = parse_fasta_index(index_file)

        seq_dict = {}
        with open(fasta_file) as f:
            header = f.readline().strip()
            header = header.lstrip('>')
            seqs = []
            for line in f:
                if line.startswith('>'):
                    seq_dict[header] = ''.join(seqs)
                    header = line.strip()
                    header = header.lstrip('>')
                    seqs = []
                else:
                    seqs.append(line.strip())
            seq_dict[header] = ''.join(seqs)
        cls._seq_dict = seq_dict

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._tmp_dir)
