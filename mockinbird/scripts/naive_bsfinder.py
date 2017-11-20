import argparse
import math
from collections import Counter
from mockinbird.utils import argparse_helper as aph
from mockinbird.utils.parsers import PC_MANDATORY_FIELDS


def create_parser():
    description = ('mockinbird-bsfinder detects protein RNA binding sites from PAR-CLIP '
                   'experiments in mpileup files based on transitions only.')
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('pileup_file', help='path to the inputfile *.pileup', type=aph.file_r)
    parser.add_argument('output_table', help='define output file *.table')
    parser.add_argument('--min_transitions', help='minimum number of transitions required',
                        type=int, default=2)
    parser.add_argument('--reference', '-r', help='set default reference nucleotide',
                        choices=['A', 'C', 'G', 'T'], default='T')
    parser.add_argument('--mutation', '-m', help='set default mutation nucleotide',
                        choices=['A', 'C', 'G', 'T'], default='C')
    return parser


def run():

    parser = create_parser()
    args = parser.parse_args()

    complement = {}
    for X, Y in [('A', 'T'), ('C', 'G'), ('G', 'C'), ('T', 'A')]:
        complement[X] = Y
        complement[X.lower()] = Y.lower()

    sense_reference = args.reference
    sense_transition = args.mutation
    antisense_reference = complement[sense_reference]
    antisense_transition = complement[sense_transition].lower()

    with open(args.pileup_file) as pileup, open(args.output_table, 'w') as outfile:

        # write output header
        print(*PC_MANDATORY_FIELDS, sep='\t', file=outfile)

        for line in pileup:
            try:
                chrom, pos, nucleotide, coverage, cov_str, _ = line.split()
            except ValueError:
                continue

            cov_toks = Counter(cov_str)
            if nucleotide == sense_reference:
                transitions = cov_toks[sense_transition]
                strand = '+'
            elif nucleotide == antisense_reference:
                transitions = cov_toks[antisense_transition]
                strand = '-'
            else:
                continue

            occupancy = math.nan
            if transitions >= args.min_transitions:
                print(chrom, pos, transitions, coverage, transitions, strand, occupancy,
                      sep='\t', file=outfile)


if __name__ == '__main__':
    run()
