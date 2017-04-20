"""
Takes PAR-CLIP sites and a genome and saves genomic sequences as fasta file
around PAR-CLIP sites according to the given parameters.

**Usage:** stammp-parclip2Fasta [-h] [-v]
                            sites genome fafile filterGFF start stop width
                            additionalFilterWidth key

**Positional arguments:**
  ===================== ======================================================
  sites                 PAR-CLIP file \*.table
  genome                path to genome
  fafile                output filename
  filterGFF             set path to GFF if sites should be removed that
                        overlap with the GFF [default = ]
  start                 start index of PAR-CLIP sites [default=0]
  stop                  stop index of PAR-CLIP sites [default=1500]
  width                 number of nt +/- the crosslink site [default=15]
  additionalFilterWidth
                        number of nt that are added to the start/stop indices
                        of the GFF annotations
  key                   set key that is used for PAR-CLIP site ordering
                        [default = 'occ'], options: ['occ', 'm', 'r', 'mr',
                        'pvalue']
  ===================== ======================================================

**Optional arguments:**
  =============         ===============================
  -h, --help            show this help message and exit
  -v, --verbose         verbose output
  =============         ===============================
"""
import argparse
import os
from stammp.obj import *
from stammp.utils import ParclipSiteContainer

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Takes PAR-CLIP sites and a genome and saves genomic sequences as fasta file around PAR-CLIP sites according to the given parameters.', epilog="contact: torkler@genzentrum.lmu.de")
    parser.add_argument('sites',      help='PAR-CLIP file *.table')
    parser.add_argument('genome',     help='path to genome')
    parser.add_argument('fafile',     help='output filename')
    parser.add_argument('filterGFF',  help='set path to GFF if sites should be removed that overlap with the GFF [default = '']', default='')
    parser.add_argument('start',      help='start index of PAR-CLIP sites [default=0]', type=int, default = 0)
    parser.add_argument('stop',       help='stop index of PAR-CLIP sites [default=1500]', type=int, default = 1500)
    parser.add_argument('width',      help='number of nt +/- the crosslink site [default=15]', type=int, default = 15)
    parser.add_argument('additionalFilterWidth', help='number of nt that are added to the start/stop indices of the GFF annotations', type=int, default = 20)
    parser.add_argument('key',  help='set key that is used for PAR-CLIP site ordering [default = \'occ\'], options: [\'occ\', \'m\', \'r\', \'mr\', \'pvalue\']', default='occ')
    parser.add_argument('-v','--verbose', dest='verbose', action="store_true", default=False, help='verbose output')
    args = parser.parse_args()
    
    yeast   = genome.Genome(args.genome, False)
    sites   = ParclipSiteContainer()
    sites.loadFromFile(args.sites)
    
    if args.verbose:
        print('#sites              : '+str(sites.size()))
    if args.filterGFF != '':
        anno     = gff.GFF(args.filterGFF)
        sites   = sites.removeSitesLocatedInGFF(anno, args.additionalFilterWidth)
        print('#sites after removal: '+str(sites.size()))
    
    sites.sort(args.key)
    sites.save2Fasta(yeast, args.fafile, args.start, args.stop, args.width)

