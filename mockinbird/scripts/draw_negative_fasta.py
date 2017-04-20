import argparse
from collections import Counter

import numpy as np

from mockinbird.utils import argparse_helper as aph
from mockinbird.utils.helper_objects import EfficientGenome


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('parclip_sites', type=aph.file_r,
                        help='path to parclip sites fasta file')
    parser.add_argument('genome_fasta', type=aph.file_r,
                        help='path to reference genome fasta file')
    parser.add_argument('gff_annotation', type=aph.file_r,
                        help='path to gff annotation')
    parser.add_argument('flanking_bp', type=int,
                        help='size of flank in nucleotides')
    parser.add_argument('n_sequences', type=int,
                        help='number of sequences')
    parser.add_argument('fasta_out', type=aph.file_rw_or_dir_rwx,
                        help='path to output fasta file')
    parser.add_argument('--retries', default=10,
                        help='number of unsuccessful draws before giving up')
    parser.add_argument('--seed', default=42, type=int)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    np.random.seed(args.seed)

    genome = EfficientGenome(args.genome_fasta)
    gff_annot = GFFAnnotation(args.gff_annotation)

    chrom_counter = Counter()
    with open(args.parclip_sites) as pc_sites:
        for line in pc_sites:
            chrom, *_ = line.split()
            chrom_counter[chrom] += 1

    total_count = sum(chrom_counter.values())
    chrom_list = []
    chrom_prob = []
    for chrom, count in chrom_counter.items():
        chrom_list.append(chrom)
        chrom_prob.append(count / total_count)

    with open(args.output_fasta) as outfile:
        sampled_seqs = 0
        while sampled_seqs < args.n_sequences:
            rand_chrom = np.random.choice(chrom_list, p=chrom_prob)
            rand_annot = gff_annot.draw_random_annotation(seqid=rand_chrom)
            if rand_annot is None:
                print('WARNING: no annotations for chromosome %s available'
                      % rand_chrom, file=sys.stderr)
            







if __name__ == '__main__':
    main()
