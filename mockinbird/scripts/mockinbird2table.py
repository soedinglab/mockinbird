import argparse
from mockinbird.utils.parsers import PC_MANDATORY_FIELDS
from math import log
from itertools import chain


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('posterior_file')
    parser.add_argument('output_file')
    parser.add_argument('--post_thresh', default=0.01, type=float)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    post_thresh = args.post_thresh

    header = PC_MANDATORY_FIELDS + ['posterior']
    with open(args.posterior_file) as infile, open(args.output_file, 'w') as outfile:
        infile.readline()
        print(*header, sep='\t', file=outfile)
        for line in infile:
            toks = line.split()
            if float(toks[-1]) < post_thresh:
                continue
            score = log(float(toks[7])) + log(float(toks[8]))
            print(*chain(toks[0:4]), score, toks[6], -1, toks[-1], sep='\t', file=outfile)


if __name__ == '__main__':
    main()
