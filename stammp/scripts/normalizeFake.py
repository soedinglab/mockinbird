#! /usr/bin/python3
"""
Observed PAR-CLIP mutations are divided by the observed PAR-CLIP coverage.

**Usage:** stammp-normalizeFake [-h] [-q QUANTILE] inputfile outputfile

**Positional arguments:**
  ==========   ======================================
  inputfile    PAR-CLIP file \*.table
  outputfile   Fake normalized PAR-CLIP file \*.table
  ==========   ======================================

**Optional arguments:**
  ===========  ====================================
  -h, --help   show this help message and exit
  ===========  ====================================
"""
import argparse
import os
import sys
from stammp.obj import *

def main(inputfile, outputfile):
    if os.path.isfile(inputfile) == False:
        print('Inputfile: '+inputfile+' does not exist')
        sys.exit(-1)
    sites = parclipsites.ParclipSites('')
    sites.loadFromFile(inputfile)
    
    for i in range(sites.size()):
       sites.occ[i] = sites.m[i]/sites.r[i]
    sites.save2File(outputfile)

def run():
    parser = argparse.ArgumentParser(description='Observed PAR-CLIP mutations are divided by the observed PAR-CLIP coverage..', epilog="contact: torkler@genzentrum.lmu.de")
    parser.add_argument('inputfile', help='PAR-CLIP file *.table')
    parser.add_argument('outputfile', help='Fake normalized PAR-CLIP file *.table')
    args = parser.parse_args()
    main(args.inputfile, args.outputfile)

if __name__ == '__main__':
    run()
