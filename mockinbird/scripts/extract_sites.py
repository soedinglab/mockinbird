import argparse

from mockinbird.utils.parsers import GFF3Parser
from mockinbird.utils.helper_objects import EfficientGenome


def create_parser():
    parser = argparse.ArgumentParser('extract_sites')
    parser.add_argument('--transition_from', default='T', choices='ACGT')
    parser.add_argument('fasta_file')
    parser.add_argument('--gff_file')
    parser.add_argument('output_file')
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    transition_map = {
        'A': ('A', 'T'),
        'C': ('C', 'G'),
        'G': ('G', 'C'),
        'T': ('T', 'A'),
    }

    pos_nuc, neg_nuc = transition_map[args.transition_from]
    with open(args.output_file, 'w') as out:
        print('chrom', 'pos', 'base', 'strand', sep='\t', file=out)
        if args.gff_file:
            with EfficientGenome(args.fasta_file) as genome:
                with open(args.gff_file) as gff_handle:
                    parser = GFF3Parser(gff_handle)
                    for rec in parser.parse():
                        seq = genome.get_sequence(rec.seqid, rec.start, rec.end, strand='+')
                        for i, base in enumerate(seq):
                            if base == pos_nuc:
                                strand = '+'
                            elif base == neg_nuc:
                                strand = '-'
                            else:
                                continue
                            print(rec.seqid, args.start + i, base, strand, sep='\t', file=out)

        else:
            with open(args.fasta_file) as fasta_handle:
                for line in fasta_handle:
                    line = line.strip()
                    if line == '':
                        continue
                    if line.startswith('>'):
                        seqid = line.lstrip('>')
                        seqid = seqid.strip()
                        pos = 0
                        continue

                    for base in line:
                        pos += 1
                        if base == pos_nuc:
                            strand = '+'
                        elif base == neg_nuc:
                            strand = '-'
                        else:
                            continue
                        print(seqid, pos, base, strand, sep='\t', file=out)


if __name__ == '__main__':
    main()
