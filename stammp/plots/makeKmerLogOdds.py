"""
Plot log-odds for all kmers as a heatmap. Bins of the x-axis are sorted according to *key*. Y-axis is sorted according to the most enriched kmers in the first 3 columns. Make sure to provide the appropriate negative set for the chosen *kmer*. 

Log-odds are calculated as: log2(p(kmer|PARCLIP)/p(kmer|background))

Colorscale from red-to-white-to-green so that log-odd = 0 is white.

Colorscale limits:
  ====    =============
  kmer    log-odd limit
  2       [-2.5,2.5]
  3       [-2.5,2.5]
  4       [-3.0,3.0]
  5       [-4.0,4.0]
  >5      [-4.0,4.0]
  ====    =============

.. warning:: Make sure to build an approriate negative set first :mod:`~stammp.scripts.makeNegSets`.

.. warning:: If you have a low number of PAR-CLIP sites in your PAR-CLIP table and you use *-q* bins will overlap which can be misleading in this case!

**Usage:** stammp-makeKmerLogOdds [-h] [--gff GFF] [--kmer KMER] [--key KEY] [-q]
                              [-v]
                              parclip outdir prefix genome negset

**Positional arguments:**
  =======          ===================================
  parclip          PAR-CLIP file \*.table
  outdir           output directory
  prefix           prefix of filenames
  genome           path to genome
  kmer             kmer length [default = 4]
  negset           path to correct k-mer negative set
  =======          ===================================

**Optional arguments:**
  ===============  ===========================================================
  -h, --help       show this help message and exit
  --gff GFF        remove PAR-CLIP sites overlapping with annotations
  --key KEY        set key that is used for PAR-CLIP site ordering [default =
                   occ], options: [occ, m, r, mr, pvalue]
  -q, --quantiles  use quantiles for binarization instead of fixed bin size.
                   Note, if you have a small number of bindng sites the bins
                   based on quantiles can overlap!
  -v, --verbose    verbose output
  ===============  ===========================================================

Example::

    $ stammp-makeKmerLogOdds path/to/parclip.table outputdir/ prefix path/to/genome.fa kmer path/to/negset.table

.. image:: img/img_plotKmerLogOdds.png
   :align: center
   :width: 341px
   :alt: alternate text
"""
import os
import math
import argparse

from stammp.obj import functions, gff, parclipsites, genome
from stammp.utils import execute
from stammp.utils import argparse_helper as aph


def loadNegTable(filename):
    with open(filename) as fc:
        negset = {}
        total = 0
        for line in fc:
            s = line.split('\t')
            negset[s[0]] = float(s[1].split('\n')[0])
            total += negset[s[0]]
    for k in negset:
        negset[k] = negset[k] / total
    return negset


def getkmerLogs(s, g, neg, kmers, start, stop, width):
    total = 0
    kmer = len(kmers[0])
    kmer_freqs = {}
    for k in kmers:
        kmer_freqs[k] = 1
        total += 1
    for i in range(start, stop):
        seq = g.getSequence(s.chrs[i], (s.pos[i] - 1 - width), (s.pos[i] - 1 + width + 1))
        if seq != -1:
            if s.strand[i] == '-':
                seq = functions.makeReverseComplement(seq)
            for j in range(len(seq) - kmer):
                kmer_freqs[seq[j:(j + kmer)]] += 1
                total += 1
    for k in kmers:
        kmer_freqs[k] = math.log(((kmer_freqs[k] / total) / neg[k]), 2)
    return(kmer_freqs)


def sortAndSave(oddlist, outfile, kmers):
    sortedKmers = []
    for k in range(len(kmers)):
        if len(oddlist) > 2:
            counts = oddlist[0][kmers[k]] + oddlist[1][kmers[k]] + oddlist[2][kmers[k]]
            sortedKmers.append([kmers[k], counts])
        else:
            sortedKmers.append([kmers[k], oddlist[0][kmers[k]]])
    sortedKmers = sorted(sortedKmers, key=lambda d: d[1])[::-1]

    with open(outfile, 'w') as fc:
        for i in range(len(sortedKmers)):
            fc.write(sortedKmers[i][0] + '\t')
            for j in range(len(oddlist)):
                fc.write(str(oddlist[j][sortedKmers[i][0]]) + '\t')
            fc.write('\n')


def main(parclip, outdir, prefix, genomepath, negset, gfffile, kmer, key, useQuantiles, verbose):
    scriptPath = os.path.dirname(os.path.realpath(__file__))
    plot_script = os.path.join(scriptPath, 'plotKmerLogOdds.R')
    pc = parclipsites.ParclipSites()
    pc.loadFromFile(parclip)
    genomeseq = genome.Genome(genomepath)

    if gfffile is not None:
        anno = gff.GFF(gfffile)
        pc.removeSitesLocatedInGFF(anno)
    pc.sort(key)

    kmers = functions.makekmers(kmer, list('ACGT'))[kmer - 1]
    negfreq = loadNegTable(negset)

    allfreqs = []
    fileprefix = '%s_logodds_%smer_sort_%s' % (prefix, kmer, key)
    if useQuantiles:
        fileprefix = fileprefix+'_quantiles'
        out_file = os.path.join(outdir, fileprefix + '.warning')
        with open(out_file, 'w') as fc:
            allfreqs.append(getkmerLogs(pc, genomeseq, negfreq, kmers, 0, 1000, 15))
            quantiles = [
                0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.125, 0.15, 0.175,
                0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4,
                0.45, 0.5, 0.55, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9
            ]
            count = 1
            stop = 1000
            for q in quantiles:
                if verbose:
                    functions.showProgress(count, len(quantiles),
                                           'Getting kmer log-odds from quantiles...')
                old_stop = stop
                start = functions.getQuantileIndex(pc.size(), q) - 500
                stop = functions.getQuantileIndex(pc.size(), q) + 500
                if start < 0:
                    start = 0
                if stop > (pc.size() - 2):
                    break
                count = count + 1
                if (stop - 500) < old_stop:
                    msg_pat = 'Bin %s and %s are overlapping by %s sites!'
                    # 2x quantiles[count - 2] is probably a bug
                    msg = msg_pat % (quantiles[count - 2], quantiles[count - 2],
                                     old_stop - (stop - 500))
                    print(msg, file=fc)
                allfreqs.append(getkmerLogs(pc, genomeseq, negfreq, kmers, start, stop, 15))
    else:
        maxsize = 50000
        stepsize = 1000
        start = 0
        stop = 1000
        run = True
        while run:
            if stop > (pc.size() - 2) or stop > maxsize:
                print()
                print('STOP at: %s' % + stop)
                run = False
                break
            if verbose:
                functions.showProgress(stop, maxsize, 'Getting kmer log-odds from bins...')
            allfreqs.append(getkmerLogs(pc, genomeseq, negfreq, kmers, start, stop, 15))
            start = stop
            stop = stop + stepsize

    table_file = os.path.join(outdir, fileprefix + '.table')
    pdf_file = os.path.join(outdir, fileprefix + '.pdf')
    sortAndSave(allfreqs, table_file, kmers)

    cmd = [
        'R',
        '-q',
        '--slave',
        '-f %r' % plot_script,
        '--args',
        '%r' % table_file,
        '%r' % pdf_file,
    ]
    execute(cmd)


def run():
    parser = argparse.ArgumentParser(
        description=('Plot log-odds for all kmers as a heatmap. Make sure to '
                     'provide the appropriate negative set.')
    )
    parser.add_argument('parclip', help='PAR-CLIP file *.table', type=aph.file_r)
    parser.add_argument('outdir', help='output directory', type=aph.dir_rwx)
    parser.add_argument('prefix', help='prefix of filenames')
    parser.add_argument('genome', help='path to genome')
    parser.add_argument('--kmer', help='kmer length [default = 4]', type=int, default=4)
    parser.add_argument('negset', help='path to correct k-mer negative set')
    parser.add_argument('--gff', help='remove PAR-CLIP sites overlapping with annotations')
    key_choices = ['occ', 'm', 'r', 'mr', 'pvalue']
    key_help = ('set key that is used for PAR-CLIP site ordering [default = occ], '
                'options: [occ, m, r, mr, pvalue]')
    parser.add_argument('--key', help=key_help, choices=key_choices, default='occ')
    quantile_help = ('use quantiles for binarization instead of fixed bin size. '
                     'Note, if you have a small number of bindng sites the bins '
                     'based on quantiles can overlap!')
    parser.add_argument('-q', '--quantiles', dest='useQuantiles', action='store_true',
                        help=quantile_help)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='verbose output')
    args = parser.parse_args()

    main(args.parclip, args.outdir, args.prefix, args.genome, args.negset,
         args.gff, args.kmer, args.key, args.useQuantiles, args.verbose)


if __name__ == '__main__':
    run()
