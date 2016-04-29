"""
Plots kmer occurences per sequence position for all kmers of a given kmer-
length for selected PAR-CLIP sites.

**Usage:** stammp-makeKmerPerPosition [-h] [--kmer KMER] [--start START]
                                  [--stop STOP] [--width WIDTH] [--key KEY]
                                  [--filterGFF FILTERGFF] [--awidth AWIDTH]
                                  [-r]
                                  sites genome outdir prefix

**Positional arguments:**
  =========             ======================
  inputfile             PAR-CLIP file \*.table
  genome                path to genome
  outdir                output directory
  prefix                prefix
  =========             ======================

**Optional arguments:**
  ====================== ======================================================
  -h, --help             show this help message and exit
  --kmer KMER            kmer-length [default=3]
  --start START          start index of PAR-CLIP sites [default=0]
  --stop STOP            stop index of PAR-CLIP sites [default=1500]
  --width WIDTH          number of nt +/- the crosslink site [default=50]
  --key KEY              set key that is used for PAR-CLIP site ordering
                         [default = occ], options: [occ, m, r, mr, pvalue]
  --filterGFF FILTERGFF
                         set path to GFF if sites should be removed that
                         overlap with the GFF. Default = '' means that no sites
                         are filtered out.
  --awidth AWIDTH        number of nt that are added to the start/stop indices
                         of the GFF annotations
  -r, --remove           remove temporary text files. [default: false]
  ====================== ======================================================

Example::

    $ stammp-makeKmerPerPosition parclip.table genome.fa output/ prefix --kmer 4 --start 0 --stop 2000 --width 50 --key occ -r


.. image:: img/img_kmerPerPosition.png
   :align: center
   :width: 700px
   :alt: alternate text
"""
import argparse
import os
from itertools import chain
from stammp.obj import functions, genome, parclipsites, gff
from stammp.utils import execute
from stammp.utils import argparse_helper as aph


def getKmerOccurences(listofsequences, outfile, kmer=3, verbose=False):
    kmers = functions.makekmers(kmer + 1, list('ACGT'))[kmer]

    with open(outfile, 'w') as file_center:
        if verbose:
            print('Searching %smer      #' % (kmer + 1))
        kmercount = 0
        for km in kmers:
            if verbose:
                print('', km, '%s/%s\r' % (kmercount, len(kmers)), sep='\t')
            counts = [0] * len(listofsequences[0])
            for s in listofsequences:
                hits = functions.findAllSubstrings(s, km)
                for h in hits:
                    counts[h] = counts[h] + 1

            print(*chain([km], counts), sep='\t', file=file_center)
            kmercount += 1
    if verbose:
        print('', km, '%s/%s' % (kmercount, len(kmers)), sep='\t')


def loadFasta(filename):
    with open(filename, 'r') as fc:
        seqs = []
        for line in fc:
            if line[0] != '>':
                seqs.append(line.split('\n')[0])
    return seqs


def run():
    scriptPath = os.path.dirname(os.path.realpath(__file__))
    plot_script = os.path.join(scriptPath, 'plotKmerPerPosition.R')

    parser = argparse.ArgumentParser(
        description=(
            'Plots kmer occurences per sequence position for all kmers of a given '
            'kmer-length for selected PAR-CLIP sites.'
        )
    )
    parser.add_argument('inputfile', help='PAR-CLIP file *.table', type=aph.file_r)
    parser.add_argument('genome', help='path to genome', type=aph.file_r)
    parser.add_argument('outdir', help='output directory', type=aph.dir_rwx)
    parser.add_argument('prefix', help='prefix')
    parser.add_argument('--kmer', help='kmer-length [default=3]', type=int, default=3)
    parser.add_argument('--start', help='start index of PAR-CLIP sites [default=0]',
                        type=int, default=0)
    parser.add_argument('--stop', help='stop index of PAR-CLIP sites [default=1500]',
                        type=int, default=1500)
    parser.add_argument('--width', help='number of nt +/- the crosslink site [default=50]',
                        type=int, default=50)
    key_choices = ['occ', 'm', 'r', 'mr', 'pvalue']
    key_help = ('set key that is used for PAR-CLIP site ordering [default = occ], '
                'options: [occ, m, r, mr, pvalue]')
    parser.add_argument('--key', help=key_help, choices=key_choices, default='occ')
    gff_help = ('set path to GFF if sites should be removed that overlap with the GFF. '
                'Default = \'\' means that no sites are filtered out.')
    parser.add_argument('--filterGFF', help=gff_help, default='')
    awidth_help = 'number of nt that are added to the start/stop indices of the GFF annotations'
    parser.add_argument('--awidth', help=awidth_help, type=int, default=20)
    parser.add_argument('-r', '--remove', dest='remove', action='store_true',
                        help='remove temporary text files. [default: false]')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    yeast = genome.Genome(location=args.genome, verbose=False)
    sites = parclipsites.ParclipSites()
    sites.loadFromFile(args.inputfile)

    if args.filterGFF != '':
        anno = gff.GFF(args.filterGFF)
        sites = sites.removeSitesLocatedInGFF(anno, args.awidth)

    sites.sort(args.key)
    seqs = sites.getSequences(yeast, args.start, args.stop, args.width)

    prefix_fmt = '%s_kmerPerPosition_kmer%s_start%s_stop%s_width%s_sort_%s'

    prefix = prefix_fmt % (args.prefix, args.kmer, args.start, args.stop,
                           args.width, args.key)
    outfile_table = os.path.join(args.outdir, prefix + '.table')
    outfile_pdf = os.path.join(args.outdir, prefix + '.pdf')
    getKmerOccurences(seqs, outfile_table, kmer=(args.kmer - 1), verbose=args.verbose)

    cmd = [
        'R',
        '-q',
        '--slave',
        '-f %s' % plot_script,
        '--args',
        outfile_table,
        outfile_pdf,
        args.width,
        0,
        args.width + 1
    ]
    execute(cmd)
    if args.remove:
        os.remove(outfile_table)


if __name__ == '__main__':
    run()
