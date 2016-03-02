#! /usr/bin/python3
"""
"""
import argparse
import os
import sys
from stammp.obj import *
import math

def main(parclipA, parclipB, anno, start, stop, width, annowidth, verbose):
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
        if anno.isInside(tmpA.chrs[i],tmpA.pos[i],tmpA.strand[i], annowidth)[1]:
            dataA.addSite(tmpA.chrs[i], tmpA.pos[i], tmpA.m[i], tmpA.r[i], tmpA.result[i], tmpA.strand[i], tmpA.occ[i])
            count +=1 
        i += 1
    coloc = 1
    count_coloc = 1
    if verbose:
        print('\n')
    #print(dataA.size())
    for i in range(dataA.size()):
        values = dataB.getValues(dataA.chrs[i], dataA.pos[i], dataA.strand[i], True, width, width)
        if values != None:
            count_coloc += 1
            coloc += max(values)
#            for v in values:
#                if v > 0:
#                    coloc += v
#                    count_coloc += 1
        if verbose:
            functions.showProgress(i, (dataA.size()-1), 'Collecting colocolization data')
    coloc = coloc / count_coloc
    if verbose:
        print('')
    return math.log( coloc/functions.getQuantile(dataB.occ,0.5) ,2)
    #return coloc

def run():
    scriptPath = os.path.dirname(os.path.realpath(__file__))
    scriptPath = scriptPath+'/'
    parser = argparse.ArgumentParser(description='Returns the colocalization valiue for a pair of PAR-CLIP tables.', epilog="contact: torkler@genzentrum.lmu.de")
    parser.add_argument('parclipA', help='path to the PAR-CLIP table *.table')
    parser.add_argument('parclipB', help='path to the PAR-CLIP table *.table')
    parser.add_argument('gff', help='path to GFF file')
    parser.add_argument('--start', help='set start index for parclipA [default: 0]', dest='start', default=0, type=int)
    parser.add_argument('--stop', help='set stop index for parclipA [default: 1000]', dest='stop', default=1000, type=int)
    parser.add_argument('-w','--width', help='set window size in nt [default: 12]', dest='width', default=12, type=int)
    parser.add_argument('--gffMIN', help='set minium annotation size [default: 500]', dest='gffMIN', default=500, type=int)
    parser.add_argument('--gffMAX', help='set maximum annotation size [default: 4000]', dest='gffMAX', default=4000, type=int)
    parser.add_argument('--gffDIST', help='set minium distance to next anno [default: 200]', dest='gffDIST', default=200, type=int)
    parser.add_argument('--annoWIDTH', help='set additional number of nt to each start stop [default: 0]', dest='annoWIDTH', default=0, type=int)
    parser.add_argument('-v','--verbose', dest='verbose', action="store_true", default=False, help='verbose output')
    args = parser.parse_args()
    
    anno = gff.GFF(args.gff)
    anno.filterSize(args.gffMIN, args.gffMAX)
    anno.takeGenesSenseMinDist(args.gffDIST, False)
    
    coloc = main(args.parclipA, args.parclipB, anno, args.start, args.stop, args.width, args.annoWIDTH, args.verbose)
    print(str(coloc))

if __name__ == '__main__':
    run()

