#! /usr/bin/python3
"""
Make random negativ sets (fasta and 2-5-mer tables) from random drawings of annotations in GFF
format. Negative sets are mandatory for k-mer log odd calculations or motif finding with XXmotif.

.. warning:: Negative sets heavily influence the results of motif analysis. Make sure to select a good negative set to avoid wrong results. For the analysis of protein-RNA interactions it is a good idea to start with a negative set based on randomly choosen sequences of the transcriptome.

**Usage:** stammp-makeNegSets [-h] [--number RNUMBER] [--width WIDTH] [-v]
                          gff genome prefix outdir

**Positional arguments:**
  gff               GFF file
  genome            path to genome
  prefix            prefix
  outdir            output directory

**Optional arguments:**
  -h, --help        show this help message and exit
  --number RNUMBER  set number or random drawings [default: 10000]
  --width WIDTH     set number or nt +/- selected position [default: 20]
  -v, --verbose     verbose output [default: false]

"""
import argparse
import os
import random
from stammp.obj import *

def getRandomSequences(anno, wg, rnumber, width):
    seq = []
    i = 0
    while i < rnumber:
        rnd_anno = random.randint(0,(anno.size()-1))
        rnd_pos  = random.randint(anno.start[rnd_anno], anno.stop[rnd_anno])
        tmp_seq = wg.getSequence(anno.chr[rnd_anno], rnd_pos-width, rnd_pos+width)
        if tmp_seq != -1:
            #print(anno.chr[rnd_anno]+' '+str(anno.start[rnd_anno])+' '+str(anno.stop[rnd_anno])+' RND: '+str(rnd_pos))
            if anno.strand[rnd_anno] == '+':
                seq.append(tmp_seq)
            else:
                seq.append(functions.makeReverseComplement(tmp_seq))
            i += 1
    return seq

def getKmerCounts(seqs, kmer=3):
    kmers = functions.makekmers(kmer+1, ['A','C','G','T'])[kmer]
    kmer_counts = {}
    for k in kmers:
        kmer_counts[k] = 1
    for s in seqs:
        for i in range(len(s)-(kmer)):
            kmer_counts[s[i:(i+kmer+1)]] += 1
    return kmer_counts

def main(gfffile, genomepath, prefix, outdir, rnumber, width, verbose):
    anno = gff.GFF(gfffile)
    g    = genome.Genome(genomepath)
    
    rnd_seqs = getRandomSequences(anno, g, rnumber, width)
    kmer_table = getKmerCounts(rnd_seqs, kmer=3)
    
    fc = open(outdir+'rnd_sequences_'+prefix+'_'+str(rnumber)+'_w'+str(width)+'.fa','w')
    for i in range(len(rnd_seqs)):
        fc.write('>rnd_seq_'+str(i)+'\n')
        fc.write(rnd_seqs[i]+'\n')
    fc.close()
    
    for i in range(1,5):
        print('Getting '+str(i+1)+'mer data...')
        kmer_table = getKmerCounts(rnd_seqs, kmer=i)
        keys = list(kmer_table.keys())
        keys.sort()
        fc = open(outdir+'rnd_sequences_'+prefix+'_'+str(rnumber)+'_w'+str(width)+'_'+str(i+1)+'mer.table','w')
        
        for k in keys:
            fc.write(k+'\t'+str(kmer_table[k])+'\n')
        fc.close()

def run():
    parser = argparse.ArgumentParser(description='Make random negativ sets (fasta and 2-5-mer tables) from random drawings of annotations in GFF format.', epilog="contact: torkler@genzentrum.lmu.de")
    parser.add_argument('gff', help='GFF file')
    parser.add_argument('genome', help='path to genome')
    parser.add_argument('prefix', help='prefix')
    parser.add_argument('outdir',   help='output directory')
    parser.add_argument('--number', help='set number or random drawings [default: 10000]', dest='rnumber', default=10000, type=int)
    parser.add_argument('--width',  help='set number or nt +/- selected position [default: 20]', dest='width', default=20, type=int)
    parser.add_argument('-v','--verbose', dest='verbose', action="store_true", default=False, help='verbose output [default: false]')
    args = parser.parse_args()
    main(args.gff, args.genome, args.prefix, args.outdir, args.rnumber, args.width, args.verbose)

if __name__ == '__main__':
    run()












