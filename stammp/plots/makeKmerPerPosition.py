#! /usr/bin/python3
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
import sys
from stammp.obj import *

def getKmerOccurences(listofsequences, outfile, kmer=3):
    kmers = functions.makekmers(kmer+1, ['A','C','G','T'])[kmer]
    
    file_center = open(outfile, 'w')
    print('Searching '+str(kmer+1)+'mer      #')
    kmercount = 0
    for km in kmers:
        sys.stdout.write('\t'+km+'\t'+str(kmercount)+'/'+str(len(kmers))+'\r')
        counts = [0]*len(listofsequences[0])
        for s in listofsequences:
            hits = functions.findAllSubstrings(s, km)
            for h in hits:
                counts[h] = counts[h] + 1
                
        s = km+'\t'
        for c in counts:
            s = s+str(c)+'\t'
        file_center.write(s+'\n')
        kmercount += 1
    file_center.close()
    sys.stdout.write('\t'+km+'\t'+str(kmercount)+'/'+str(len(kmers))+'\n')

def loadFasta(filename):
    fc = open(filename, 'r')
    line = fc.readline()
    seqs = []
    while line:
        if line[0] != '>':
            seqs.append(line.split('\n')[0])
        line = fc.readline()
    fc.close()
    return(seqs)

def run():
    scriptPath = os.path.dirname(os.path.realpath(__file__))
    scriptPath = scriptPath+'/'
    parser = argparse.ArgumentParser(description='Plots kmer occurences per sequence position for all kmers of a given kmer-length for selected PAR-CLIP sites.', epilog="contact: torkler@genzentrum.lmu.de")
    parser.add_argument('inputfile',      help='PAR-CLIP file *.table')
    parser.add_argument('genome',     help='path to genome')
    parser.add_argument('outdir',     help='output directory')
    parser.add_argument('prefix',     help='prefix')
    parser.add_argument('--kmer',       help='kmer-length [default=3]', type=int, default = 3)
    parser.add_argument('--start',      help='start index of PAR-CLIP sites [default=0]', type=int, default = 0)
    parser.add_argument('--stop',       help='stop index of PAR-CLIP sites [default=1500]', type=int, default = 1500)
    parser.add_argument('--width',      help='number of nt +/- the crosslink site [default=50]', type=int, default = 50)
    parser.add_argument('--key',  help='set key that is used for PAR-CLIP site ordering [default = occ], options: [occ, m, r, mr, pvalue]', default='occ')
    parser.add_argument('--filterGFF',  help='set path to GFF if sites should be removed that overlap with the GFF. Default = \'\' means that no sites are filtered out.', default='')
    parser.add_argument('--awidth', help='number of nt that are added to the start/stop indices of the GFF annotations', type=int, default = 20)
    parser.add_argument('-r','--remove', dest='remove', action="store_true", default=False, help='remove temporary text files. [default: false]')
    args = parser.parse_args()
    
    yeast   = genome.Genome(args.genome, False)
    sites   = parclipsites.ParclipSites('')
    sites.loadFromFile(args.inputfile)
    
    if args.filterGFF != '':
        anno     = gff.GFF(args.filterGFF)
        sites   = sites.removeSitesLocatedInGFF(anno, args.awidth)
    
    sites.sort(args.key)
    seqs = sites.getSequences(yeast, args.start, args.stop, args.width)
    outfile_table = args.outdir+args.prefix+'_kmerPerPosition_kmer'+str(args.kmer)+'_start'+str(args.start)+'_stop'+str(args.stop)+'_width'+str(args.width)+'_sort'+str(args.key)+'.table'
    outfile_pdf   = args.outdir+args.prefix+'_kmerPerPosition_kmer'+str(args.kmer)+'_start'+str(args.start)+'_stop'+str(args.stop)+'_width'+str(args.width)+'_sort_'+str(args.key)+'.pdf'
    getKmerOccurences(seqs, outfile_table, kmer=(args.kmer-1))
    os.system('R -q --slave -f '+scriptPath+'plotKmerPerPosition.R --args '+outfile_table+' '+outfile_pdf+' '+str(args.width)+' 0 '+str(args.width+1))
    if args.remove:
        os.remove(outfile_table)

if __name__ == '__main__':
    run()





