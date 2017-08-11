import argparse
import pysam


def register_arguments(parser):
    parser.add_argument('input_bam', help='path to paired-end bam file')
    parser.add_argument('output_bam', help='path to output bam file')


def main():
    parser = argparse.ArgumentParser()
    register_arguments(parser)
    args = parser.parse_args()
    run(args)


def run(args):
    with pysam.AlignmentFile(args.input_bam, 'rb') as infile, \
        pysam.AlignmentFile(args.flipped_bam, 'wb', template=infile) as outfile:
        for ali in infile:
            if not ali.is_proper_pair:
                    continue
            if ali.is_read2:
                    ali.is_reverse = not ali.is_reverse
            else:
                    ali.mate_is_reverse = not ali.mate_is_reverse
            outfile.write(ali)


if __name__ == '__main__':
        main()
