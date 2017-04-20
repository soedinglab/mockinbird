"""
Make random negative sets (fasta and 2-5-mer tables) from random drawings of annotations in GFF
format. Negative sets are mandatory for k-mer log odd calculations or motif finding with XXmotif.

.. warning:: Negative sets heavily influence the results of motif analysis. Make sure to select a good negative set to avoid wrong results. For the analysis of protein-RNA interactions it is a good idea to start with a negative set based on randomly choosen sequences of the transcriptome.
"""
import argparse
import random
import os
from mockinbird.obj import gff, functions
from mockinbird.utils import argparse_helper as aph
from mockinbird.utils import EfficientGenome


def getRandomSequences(anno, wg, rnumber, width):
    seq = []
    i = 0
    while i < rnumber:
        rnd_anno = random.randint(0, (anno.size() - 1))
        rnd_pos = random.randint(anno.start[rnd_anno], anno.stop[rnd_anno])
        try:
            tmp_seq = wg.get_sequence(anno.chr[rnd_anno], rnd_pos - width, rnd_pos + width,
                                      anno.strand)
            if 'N' in tmp_seq:
                continue
            seq.append(tmp_seq)
            i += 1
        except ValueError:
            continue
    return seq


def getKmerCounts(seqs, kmer=3):
    kmers = functions.makekmers(kmer + 1, list('ACGT'))[kmer]
    kmer_counts = {}
    for k in kmers:
        kmer_counts[k] = 1
    for s in seqs:
        for i in range(len(s) - kmer):
            kmer_counts[s[i:(i + kmer + 1)]] += 1
    return kmer_counts


def main(gfffile, genomepath, prefix, outdir, rnumber, width, verbose):
    anno = gff.GFF(gfffile)
    with EfficientGenome(genomepath) as g:
        rnd_seqs = getRandomSequences(anno, g, rnumber, width)

    kmer_table = getKmerCounts(rnd_seqs, kmer=3)

    basename = 'rnd_sequences_%s_%s_w%s' % (prefix, rnumber, width)
    out_file = os.path.join(outdir, basename + '.fa')
    with open(out_file, 'w') as fc:
        for i, seq in enumerate(rnd_seqs):
            print('>rnd_seq_%s' % i, file=fc)
            print(seq, file=fc)

    for i in range(1, 5):
        kmer_table = getKmerCounts(rnd_seqs, kmer=i)
        keys = list(kmer_table.keys())
        keys.sort()
        kmer_bname = 'rnd_sequences_%s_%s_w%s_%smer' % (prefix, rnumber, width, i + 1)
        kmer_fn = os.path.join(outdir, kmer_bname + '.table')
        with open(kmer_fn, 'w') as fc:
            for k in keys:
                print(k, kmer_table[k], sep='\t', file=fc)


def create_parser():
    parser = argparse.ArgumentParser(
        description=('Make random negativ sets (fasta and 2-5-mer tables) from '
                     'random sampling of annotations in GFF format.')
    )
    parser.add_argument('gff', help='GFF file', type=aph.file_r)
    parser.add_argument('genome', help='path to genome', type=aph.file_r)
    parser.add_argument('prefix', help='prefix')
    parser.add_argument('outdir', help='output directory', type=aph.dir_rwx_create)
    parser.add_argument('--number', help='set number or random drawings',
                        dest='rnumber', default=10000, type=int)
    parser.add_argument('--width', help='set number or nt +/- selected position',
                        dest='width', default=20, type=int)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='verbose output')
    return parser

def run():
    parser = create_parser()
    args = parser.parse_args()
    main(args.gff, args.genome, args.prefix, args.outdir, args.rnumber, args.width, args.verbose)


if __name__ == '__main__':
    run()
