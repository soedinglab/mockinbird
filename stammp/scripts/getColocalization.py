#! /usr/bin/python3
"""
Returns the colocalization value for a pair of PAR-CLIP tables A and B.

**Usage:** stammp-getColocalization [-h] [--gff GFF] [--start START] [--stop STOP]
                                [-w WIDTH] [--gffMIN GFFMIN] [--gffMAX GFFMAX]
                                [--gffDIST GFFDIST] [--annoWIDTH ANNOWIDTH]
                                [-r] [-v]
                                parclipA parclipB

**Positional arguments:**
  ========              ==================================
  parclipA              path to the PAR-CLIP table \*.table
  parclipB              path to the PAR-CLIP table \*.table
  ========              ==================================

**Optional arguments:**
  ======================== ===============================================
  -h, --help               show this help message and exit
  --gff GFF                path to GFF file
  --start START            set start index for parclipA [default: 0]
  --stop STOP              set stop index for parclipA [default: 1000]
  -w WIDTH, --width WIDTH
                           set window size in nt [default: 12]
  --gffMIN GFFMIN          set minium annotation size [default: 500]
  --gffMAX GFFMAX          set maximum annotation size [default: 4000]
  --gffDIST GFFDIST        set minium distance to next anno [default: 200]
  --annoWIDTH ANNOWIDTH
                           set additional number of nt to each start stop
                           [default: 0]
  -r, --ratio              verbose output
  -v, --verbose            verbose output
  ======================== ===============================================
"""
import argparse
import os
import sys
from stammp.obj import *
import math

def main(parclipA, parclipB, start, stop, width, anno=None, annowidth=100,
         logRatio=False, verbose=False):
    tmpA     = parclipsites.ParclipSites('')
    dataB    = parclipsites.ParclipSites('')
    tmpA.loadFromFile(parclipA)
    tmpA.sort(key='occ')
    dataB.loadFromFile(parclipB)
    if start < 0 or stop < start or stop >= tmpA.size():
        print('Bullshit start and stop indices. Come on! Concentrate!')
        sys.exit()
    dataA = parclipsites.ParclipSites('')
    total = stop - start
    count = 0
    i = start
    while count < total and i < (tmpA.size()-1):
        if verbose:
            functions.showProgress(count,total-1,'Selecting PAR-CLIP sites')
        if anno == None:
            dataA.addSite(tmpA.chrs[i], tmpA.pos[i], tmpA.m[i], tmpA.r[i],
                          tmpA.result[i], tmpA.strand[i], tmpA.occ[i])
            count +=1
        else:
            if anno.isInside(tmpA.chrs[i], tmpA.pos[i], tmpA.strand[i], 
                             annowidth, annowidth)[1]:
                dataA.addSite(tmpA.chrs[i], tmpA.pos[i], tmpA.m[i], tmpA.r[i],
                              tmpA.result[i], tmpA.strand[i], tmpA.occ[i])
                count +=1
        i += 1
    coloc = 1
    count_coloc = 1
    if verbose:
        print('\n')
    for i in range(dataA.size()):
        values = dataB.getValues(dataA.chrs[i], dataA.pos[i], dataA.strand[i], True, width, width)
        if values != None:
            count_coloc += 1
            coloc += max(values)
        if verbose:
            functions.showProgress(i, (dataA.size()-1), 'Collecting colocolization data')
    coloc = coloc / count_coloc
    if verbose:
        print('')
    if logRatio:
        return math.log( coloc/functions.getQuantile(dataB.occ,0.5) ,2)
    else:
        return coloc

def run():
    scriptPath = os.path.dirname(os.path.realpath(__file__))
    scriptPath = scriptPath+'/'
    parser = argparse.ArgumentParser(description='Returns the colocalization value for a pair of PAR-CLIP tables A and B.', epilog="contact: torkler@genzentrum.lmu.de")
    parser.add_argument('parclipA', help='path to the PAR-CLIP table *.table')
    parser.add_argument('parclipB', help='path to the PAR-CLIP table *.table')
    parser.add_argument('--gff', help='path to GFF file', dest='gff', default=None)
    parser.add_argument('--start', help='set start index for parclipA [default: 0]', dest='start', default=0, type=int)
    parser.add_argument('--stop', help='set stop index for parclipA [default: 1000]', dest='stop', default=1000, type=int)
    parser.add_argument('-w','--width', help='set window size in nt [default: 12]', dest='width', default=12, type=int)
    parser.add_argument('--gffMIN', help='set minium annotation size [default: 500]', dest='gffMIN', default=500, type=int)
    parser.add_argument('--gffMAX', help='set maximum annotation size [default: 4000]', dest='gffMAX', default=4000, type=int)
    parser.add_argument('--gffDIST', help='set minium distance to next anno [default: 200]', dest='gffDIST', default=200, type=int)
    parser.add_argument('--annoWIDTH', help='set additional number of nt to each start stop [default: 0]', dest='annoWIDTH', default=0, type=int)
    parser.add_argument('-r','--ratio', dest='ratio', action="store_true", default=False, help='verbose output')
    parser.add_argument('-v','--verbose', dest='verbose', action="store_true", default=False, help='verbose output')
    args = parser.parse_args()
    
    if args.gff == None:
        anno = None
    else:
        anno = gff.GFF(args.gff)
        anno.filterSize(args.gffMIN, args.gffMAX)
        anno.takeGenesSenseMinDist(args.gffDIST, False)
    
    coloc = main(args.parclipA, args.parclipB, args.start, args.stop,
                 args.width, anno, args.annoWIDTH, args.ratio, args.verbose)
    print(str(coloc))

if __name__ == '__main__':
    run()

