import argparse
import json
import sys

import numpy as np
import pysam

from mockinbird.utils.parsers import GFF3Parser


def create_parser():
    parser = argparse.ArgumentParser('estimate_bam_statistics')
    parser.add_argument('bam_file')
    parser.add_argument('output_json')
    parser.add_argument('--gff_file')
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    total_reads = 0
    total_coverage = 0
    with pysam.AlignmentFile(args.bam_file) as reads:
        if args.gff_file:
            with open(args.gff_file) as gff_handle:
                parser = GFF3Parser(gff_handle)
                for rec in parser.parse():
                    total_reads += reads.fetch(rec.seqid, rec.start, rec.end)
                    # sums up coverage arrays, one for each base
                    for cov_arr in reads.count_coverage(rec.seqid, rec.start, rec.end):
                        total_coverage += np.sum(cov_arr)

        else:
            for ali in reads.fetch():
                total_reads += 1
                total_coverage += ali.reference_length

    if total_coverage == 0:
        print('Error: zero coverage counted. The bam index file might be '
              'corrupt. Please investigate.', file=sys.stderr)
        sys.exit(1)
    output = {
        'total_reads': total_reads,
        'total_coverage': total_coverage,
    }

    with open(args.output_json, 'w') as json_handle:
        json.dump(output, json_handle)


if __name__ == '__main__':
    main()
