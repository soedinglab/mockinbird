import argparse
import re

import pysam


def create_parser():
    parser = argparse.ArgumentParser('filter_edge_mutations')
    parser.add_argument('calmd_bam_file')
    parser.add_argument('output_bam_file')
    parser.add_argument('--edge_nucleotides', type=int, default=3)
    return parser


def filter_edge_mutations(input_bam, output_bam, nbases):
    mm_pat = re.compile('([0-9]+)[ACGT][0-9]+', re.IGNORECASE)
    with pysam.AlignmentFile(input_bam, 'rb') as infile:
        with pysam.AlignmentFile(output_bam, 'wb', template=infile) as outfile:
            for entry in infile:
                if entry.is_unmapped:
                    continue
                tag_dict = dict(entry.tags)
                if 'NM' in tag_dict:
                    nm = tag_dict['NM']
                elif 'nM' in tag_dict:
                    nm = tag_dict['nM']
                else:
                    raise ValueError('NM field not found.')

                if nm == 0:
                    continue

                if 'MD' not in tag_dict:
                    raise ValueError('MD field not found.')
                mm_hit = mm_pat.match(tag_dict['MD'])
                assert mm_hit is not None
                mm_pos = int(mm_hit.group(1))

                read_len = len(entry.seq)
                if nbases <= mm_pos < read_len - nbases:
                    outfile.write(entry)


if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    filter_edge_mutations(args.calmd_bam_file, args.output_bam_file, args.edge_nucleotides)
