import argparse
import itertools
import os
from collections import Counter
from stammp.utils.argparse_helper import dir_rwx, file_r
import pysam


CIGAR_CLIP = 4


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('bam_file', type=file_r)
    parser.add_argument('--top_n', default=10, type=int)
    parser.add_argument('output_dir', type=dir_rwx)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    hits_5p = []
    hits_3p = []
    with pysam.AlignmentFile(args.bam_file, 'rb') as infile:
        for ali in infile:
            first_op, first_bp = ali.cigartuples[0]
            if first_op == CIGAR_CLIP:
                hits_5p.append(ali.seq[:first_bp])
            last_op, last_bp = ali.cigartuples[-1]
            if last_op == CIGAR_CLIP:
                hits_3p.append(ali.seq[-last_bp:])

    len_sort = lambda x: len(x)
    for filename, hits in [('5prime_clipped.txt', hits_5p), ('3prime_clipped.txt', hits_3p)]:
        outfile = os.path.join(args.output_dir, filename)
        with open(outfile, 'w') as outhandle:
            for clip_size, seqs in itertools.groupby(sorted(hits, key=len_sort), len_sort):
                seqs = list(seqs)
                print('clipped bases: %s | clipped sequences: %s' % (clip_size, len(seqs)),
                      file=outhandle)
                print(file=outhandle)
                counter = Counter(seqs)
                for clip_seq, seq_count in counter.most_common(args.top_n):
                    print('{:30} {:10}'.format(clip_seq, seq_count), file=outhandle)
                print(file=outhandle)
            print(file=outhandle)


if __name__ == '__main__':
    main()
