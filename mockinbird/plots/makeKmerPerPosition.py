"""
Plots kmer occurences per sequence position for all kmers of a given kmer-
length for selected PAR-CLIP sites.

**Usage:** mockinbird-makeKmerPerPosition [-h] [--kmer KMER] [--start START]
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

    $ mockinbird-makeKmerPerPosition parclip.table genome.fa output/ prefix --kmer 4 --start 0 --stop 2000 --width 50 --key occ -r


.. image:: img/img_kmerPerPosition.png
   :align: center
   :width: 700px
   :alt: alternate text
"""
import argparse
import os
from itertools import chain
from mockinbird.obj import functions
from mockinbird.utils import execute
from mockinbird.utils import argparse_helper as aph
from mockinbird.utils import ParclipSiteContainer, EfficientGenome
from mockinbird.utils.postprocess_modules import sort_keys


def create_parser():
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
    parser.add_argument('--kmer', help='kmer-length', type=int, default=3)
    parser.add_argument('--start', help='start index of PAR-CLIP sites',
                        type=int, default=0)
    parser.add_argument('--stop', help='stop index of PAR-CLIP sites',
                        type=int, default=1500)
    parser.add_argument('--width', help='number of nt +/- the crosslink site',
                        type=int, default=50)
    key_help = ('set key that is used for PAR-CLIP site ordering')
    parser.add_argument('--key', help=key_help, choices=sort_keys, default='occupancy')
    gff_help = ('set path to GFF if sites should be removed that overlap with the GFF. '
                'By default no sites are filtered out.')
    parser.add_argument('--filterGFF', help=gff_help, default='')
    awidth_help = 'number of nt that are added to the start/stop indices of the GFF annotations'
    parser.add_argument('--awidth', help=awidth_help, type=int, default=20)
    parser.add_argument('--remove', '-r', action='store_true',
                        help='remove temporary text files')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='verbose output')
    return parser


def getKmerOccurences(listofsequences, width, outfile, kmer=3, verbose=False):
    kmers = functions.makekmers(kmer + 1, list('ACGT'))[kmer]

    with open(outfile, 'w') as file_center:
        if verbose:
            print('Searching %smer      #' % (kmer + 1))
        kmercount = 0
        for km in kmers:
            if verbose:
                print('', km, '%s/%s\r' % (kmercount, len(kmers)), sep='\t')
            counts = [0] * width
            for s in listofsequences:
                assert len(s) == width
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

    parser = create_parser()
    args = parser.parse_args()

    sites = ParclipSiteContainer.from_file(args.inputfile)

    if args.filterGFF != '':
        sites.remove_gff_sites(args.filterGFF, args.awidth)

    sites.sort(by=args.key, ascending=False)

    with EfficientGenome(args.genome) as genome:
        sites = sites[args.start:args.stop]
        seqs = sites.get_all_sequences(genome, args.width)

    prefix_fmt = '%s_kmerPerPosition_kmer%s_start%s_stop%s_width%s_sort_%s'

    prefix = prefix_fmt % (args.prefix, args.kmer, args.start, args.stop,
                           args.width, args.key)
    outfile_table = os.path.join(args.outdir, prefix + '.table')
    outfile_pdf = os.path.join(args.outdir, prefix + '.pdf')
    seq_len = 2 * args.width + 1
    getKmerOccurences(seqs, seq_len, outfile_table, kmer=(args.kmer - 1), verbose=args.verbose)

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
