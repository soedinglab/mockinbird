#! /usr/bin/python3
"""
Returns the processing index for a pileup and an annotation file.

**Usage:** stammp-getProcessingIndex [-h] [-s SPACE] [-w WIDTH] [-v] parclip gff

**Positional arguments:**
  =======               ===========================================
  pileup                path to the pileup \*.mpileup
  gff                   GFF file. Stop positions are used for index
                        calculation.
  =======               ===========================================

**Optional arguments:**
  ======================= =====================================================
  -h, --help              show this help message and exit
  -s SPACE, --space SPACE
                          set upstream and downstream distance to stop position
                          in nt [default: 50nt]
  -w WIDTH, --width WIDTH
                          set window size in nt [default: 25nt]
  -v, --verbose           verbose output
  ======================= =====================================================

Example::
    
    stammp-makeCenterBothEnds parclip.table annotation.gff
    
"""
import argparse
import os
import sys
from stammp.obj import *
import math

def main(pileupfile, anno, space, width, verbose):
    p    = pileup.Pileup(pileupfile)
    sum_up = 1
    sum_do = 1
    count  = 0
    for g in range(anno.size()):
        if verbose:
            functions.showProgress(g, (anno.size()-1), 'Collecting read data from annotation file')
        if anno.strand[g] == '+':
            values_up = p.getValues(anno.chr[g], (anno.stop[g]-space), anno.strand[g], width, width)
            values_do = p.getValues(anno.chr[g], (anno.stop[g]+space), anno.strand[g], width, width)
        else:
            values_up = p.getValues(anno.chr[g], (anno.start[g]+space), anno.strand[g], width, width)
            values_do = p.getValues(anno.chr[g], (anno.start[g]-space), anno.strand[g], width, width)
        if values_up != None and values_do != None:
            count  += 1
            sum_up += sum(values_up)/len(values_up)
            sum_do += sum(values_do)/len(values_do)
    sum_up = sum_up / count
    sum_do = sum_do / count
    if verbose:
        print('')
    return math.log( (sum_do/max([1,(sum_up-sum_do)])) ,2)

def run():
    scriptPath = os.path.dirname(os.path.realpath(__file__))
    scriptPath = scriptPath+'/'
    parser = argparse.ArgumentParser(description='Returns the processing index for a pileup and an annotation file.', epilog="contact: torkler@genzentrum.lmu.de")
    parser.add_argument('pileup', help='path to the pileup file *.mpileup')
    parser.add_argument('gff', help='GFF file. Stop positions are used for index calculation.')
    parser.add_argument('-s','--space', help='set upstream and downstream distance to stop position in nt [default: 50nt]', dest='space', default=50, type=int)
    parser.add_argument('-w','--width', help='set window size in nt [default: 25nt]', dest='width', default=25, type=int)
    parser.add_argument('-v','--verbose', dest='verbose', action="store_true", default=False, help='verbose output')
    args = parser.parse_args()
    anno = gff.GFF(args.gff)
    pindex = main(args.pileup, anno, args.space, args.width, args.verbose)
    print(str(pindex))

if __name__ == '__main__':
    run()
