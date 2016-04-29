#! /usr/bin/python3
"""
Set maximal occupancy to the specified quantile.

**Usage:** stammp-convert2quantile [-h] [-q QUANTILE] inputfile outputfile

**Positional arguments:**
  ==========   ====================================================
  inputfile    PAR-CLIP file \*.table [output from stammp-normalize]
  outputfile   Converted PAR-CLIP file \*.table
  ==========   ====================================================

**Optional arguments:**
  ===========  ====================================
  -h, --help   show this help message and exit
  -q QUANTILE  Set quantile [0,1.0] [Default: 0.95]
  ===========  ====================================
"""
import argparse
import os
import sys
from stammp.obj import *

def main(inputfile, outputfile, q):
    if os.path.isfile(inputfile) == False:
        print('Inputfile: '+inputfile+' does not exist')
        sys.exit(1)
    if q < 0 or q > 1.0:
        print('q: '+str(q)+' must between [0,1]')
        sys.exit(1)
    sites = parclipsites.ParclipSites()
    sites.loadFromFile(inputfile)
    # dirty hack to avoid errors on empty files
    if len(sites.chrs) > 0:
        maxocc = functions.getQuantile(sites.occ, q)
    
    for i in range(sites.size()):
        if sites.occ[i] > maxocc:
            sites.occ[i] = maxocc
    sites.save2File(outputfile)

def run():
    parser = argparse.ArgumentParser(description='Set maximal occupancy to the specified quantile.', epilog="contact: torkler@genzentrum.lmu.de")
    parser.add_argument('inputfile', help='PAR-CLIP file *.table [output from stammp-normalize]')
    parser.add_argument('outputfile', help='Converted PAR-CLIP file *.table')
    parser.add_argument('-q', help='Set quantile [0,1.0] [Default: 0.95]', default=0.95, type=float, dest='quantile') 
    args = parser.parse_args()
    main(args.inputfile, args.outputfile, args.quantile)

if __name__ == '__main__':
    run()
