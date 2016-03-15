#! /usr/bin/python3
'''
Plots all possible conditional mutation probabilities based on pileup data.

*Usage:* stammp-makeNucleotideProbabilities [-h] [-c COVERAGE] [-v] [-r]
                                            inputfile outdir prefix

*Positional arguments:*
  =========      ====================
  inputfile      path to the *.pileup
  outdir         output directory
  prefix         prefix of filenames
  =========      ====================

*Optional arguments:*
  =============  ========================================
  -h, --help     show this help message and exit
  -c COVERAGE    minium coverage [default: 5]
  -v, --verbose  verbose output
  -r, --remove   remove temporary files. [default: false]
  =============  ========================================
'''
import argparse
import math
import os
from stammp.obj import *

def getCountMat(file_pileup, minCoverage,verbose):
    alphabet = ['A','C','G','T']
    translate = {'A':0, 'C':1, 'G':2, 'T':3}
    
    if verbose:
        lines = 0
        for line in open(file_pileup):
            lines += 1
    
    mat = [[0]*4,[0]*4,[0]*4,[0]*4]
    file_pileup = open(file_pileup, 'r')
    line        = file_pileup.readline()
    count = percent_old = percent_new = 0
    if verbose: functions.showProgress(count,lines,'Processing Pileup')
    while(line):
        count += 1
        split = line.split('\t')
        
        tmp_counts = functions.getCounts(split[4], forward=True)
        if tmp_counts[0] >= minCoverage:
            for c in alphabet:
                if c == split[2]:
                    mat[translate[split[2]]][translate[c]] += tmp_counts[0]-tmp_counts[2]
                else:
                    mat[translate[split[2]]][translate[c]] += tmp_counts[1][c]
        if verbose:
            percent_new = math.trunc((count/lines)*100)
            if(percent_new > percent_old):
                functions.showProgress(count,lines,'Processing Pileup')
                percent_old = percent_new
        line = file_pileup.readline()
    if verbose: print('')
    return mat

def saveMat(mat, filename):
    fc = open(filename, 'w')
    for i in range(len(mat)):
        for j in range(len(mat[i])):
            fc.write(str(mat[i][j])+'\t')
        fc.write('\n')
    fc.close()

def main(infile, outfile, coverage, verbose):
    m = getCountMat(infile, coverage, verbose)
    saveMat(m, outfile)

def run():
    scriptPath = os.path.dirname(os.path.realpath(__file__))
    scriptPath = scriptPath+'/'
    parser = argparse.ArgumentParser(description='Plots all possible'
    + ' conditional mutation probabilities based on pileup data.', 
    epilog="contact: torkler@genzentrum.lmu.de")
    parser.add_argument('inputfile',  help='path to the *.pileup')
    parser.add_argument('outdir', help='output directory')
    parser.add_argument('prefix', help='prefix of filenames')
    parser.add_argument('-c', help='minium coverage [default: 5]', 
                        default=5, type=int, dest='coverage')
    parser.add_argument('-l', help='y-axis limit [default: 0]', 
                        default=0.0, type=float, dest='limit')
    parser.add_argument('-v','--verbose', dest='verbose', action="store_true",
                        default=False, help='verbose output')
    parser.add_argument('-r','--remove', dest='remove', action="store_true",
                        default=False, help='remove temporary files. [default: false]')
    args = parser.parse_args()
    
    functions.checkExistence(args.inputfile)
    functions.checkPath(args.outdir)
    outfile = args.outdir + args.prefix + '_nuc_mutations.table'
    outfile_img = args.outdir + args.prefix + '_nuc_mutations.pdf'
    main(args.inputfile, outfile, args.coverage, args.verbose)
    os.system('R -q --slave -f ' + scriptPath + 'plotNucleotideProbabilities.R'\
              + ' --args ' + outfile + ' ' + outfile_img + ' ' + str(args.limit))
    if args.remove:
        os.remove(outfile)

if __name__ == '__main__':
    run()