import argparse
from mockinbird.utils.helper_objects import ParclipSiteContainer, EfficientGenome


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('table_file')
    parser.add_argument('genome_fasta')
    parser.add_argument('output_fasta')
    parser.add_argument('--flank_bp', type=int, default=12)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    table = ParclipSiteContainer.from_file(args.table_file)
    with EfficientGenome(args.genome_fasta) as genome:
        table.save2Fasta(genome, args.output_fasta, width=args.flank_bp)


if __name__ == '__main__':
    main()
