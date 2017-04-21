import argparse
from collections import Counter
import sys


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('pileup_file')
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    print('chrom', 'pos', 'strand', 'k', 'n', sep='\t')
    with open(args.pileup_file) as infile:
        for line in infile:
            try:
                chrom, pos, ref_nuc, cov, cov_str, qual_str = line.split()
            except ValueError:
                continue
            if ref_nuc not in ('T', 'A'):
                continue
            cov_symb = Counter(cov_str)
            if ref_nuc == 'T':
                k = cov_symb['C']
                n = k + cov_symb['.']
                strand = '+'
            elif ref_nuc == 'A':
                k = cov_symb['g']
                n = k + cov_symb[',']
                strand = '-'
            print(chrom, pos, strand, k, n, sep='\t')


if __name__ == '__main__':
    main()
