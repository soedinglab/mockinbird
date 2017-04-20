"""
Plot all possible conditional mutation probabilities based on pileup data.

**Usage:** stammp-makeNucleotideProbabilities [-h] [-c COVERAGE] [-v] [-r]
                                            inputfile outdir prefix

**Positional arguments:**
  =========      ====================
  inputfile      path to the \*.pileup
  outdir         output directory
  prefix         prefix of filenames
  =========      ====================

**Optional arguments:**
  =============  ========================================
  -h, --help     show this help message and exit
  -c INT         minimum coverage [default: 5]
  -v, --verbose  verbose output
  -r, --remove   remove temporary files. [default: false]
  =============  ========================================

.. image:: img/img_nuc_probabilities.png
    :align: center
    :height: 250px
    :alt: alternate text
"""
import argparse
import math
import os
from stammp.obj import functions
from stammp.utils import argparse_helper as aph
from stammp.utils import execute


def create_parser():
    parser = argparse.ArgumentParser(
        description=('Plots all possible conditional mutation probabilities '
                     'based on pileup data.')
    )
    parser.add_argument('inputfile', help='path to the *.pileup', type=aph.file_r)
    parser.add_argument('outdir', help='output directory', type=aph.dir_rwx)
    parser.add_argument('prefix', help='prefix of filenames')
    parser.add_argument('--coverage', '-c', help='minimum coverage',
                        default=5, type=int)
    parser.add_argument('--limit', '-l', help='y-axis limit',
                        default=0.0, type=float)
    parser.add_argument('--verbose', '-v', action='store_true', help='verbose output')
    parser.add_argument('--remove', '-r', action='store_true',
                        help='remove temporary files')
    return parser


def getCountMat(file_pileup, minCoverage, verbose):
    alphabet = ['A', 'C', 'G', 'T']
    translate = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    if verbose:
        lines = 0
        with open(file_pileup) as infile:
            for line in infile:
                lines += 1

    mat = [[0] * 4, [0] * 4, [0] * 4, [0] * 4]
    with open(file_pileup) as file_pileup:
        count = 0
        percent_old = 0
        percent_new = 0

        if verbose:
            functions.showProgress(count, lines, 'Processing Pileup')

        for line in file_pileup:
            count += 1
            split = line.split('\t')
            nuc = split[2].upper()
            if nuc != 'N':
                tmp_counts = functions.getCounts(split[4], forward=True)
                if tmp_counts[0] >= minCoverage:
                    for c in alphabet:
                        if c == nuc:
                            mat[translate[nuc]][translate[c]] += tmp_counts[0] - tmp_counts[2]
                        else:
                            mat[translate[nuc]][translate[c]] += tmp_counts[1][c]
            if verbose:
                percent_new = math.trunc(count / lines * 100)
                if percent_new > percent_old:
                    functions.showProgress(count, lines, 'Processing Pileup')
                    percent_old = percent_new

        if verbose:
            print()

        return mat


def saveMat(mat, filename):
    with open(filename, 'w') as fc:
        for row in mat:
            print(*row, sep='\t', file=fc)


def main(infile, outfile, coverage, verbose):
    m = getCountMat(infile, coverage, verbose)
    saveMat(m, outfile)


def run():
    scriptPath = os.path.dirname(os.path.realpath(__file__))
    plot_script = os.path.join(scriptPath, 'plotNucleotideProbabilities.R')

    parser = create_parser()
    args = parser.parse_args()

    functions.checkExistence(args.inputfile)
    functions.checkPath(args.outdir)
    outfile = os.path.join(args.outdir, args.prefix + '_nuc_mutations.table')
    outfile_img = os.path.join(args.outdir, args.prefix + '_nuc_mutations.pdf')
    main(args.inputfile, outfile, args.coverage, args.verbose)

    cmd = [
        'R',
        '-q',
        '--slave',
        '-f %r' % plot_script,
        '--args',
        '%r' % outfile,
        '%r' % outfile_img,
        args.limit,
    ]
    execute(cmd)

    if args.remove:
        os.remove(outfile)


if __name__ == '__main__':
    run()
