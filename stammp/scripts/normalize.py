#! /usr/bin/python3
"""
Command line:
-------------

**Usage:** stammp-normalize [-h] [-s SWITCH] [-v] [-r SNP]
                        inputfile outputfile rnaseqfiles chrnames

Takes all PAR-CLIP sites and traverses through given pileups to get strand
specific coverage of all given pileups and divides the PAR-CLIP mutations
counts by the sum of the coverages.

**Positional arguments:**
  ===========    =============================================================
  inputfile      PAR-CLIP file \*.table
  outputfile     Normalized PAR-CLIP file \*.table
  rnaseqfiles    Comma separated list of pileup files used for normalization
                 [no whitespaces]
  chrnames       Comma separated, ordered list of chrnames [no whitespaces]
  ===========    =============================================================

**Optional arguments:**
  =============  =============================================================
  -h, --help     show this help message and exit
  -s SWITCH      Comma sperated list of 0 or 1 indicating which files have to
                 be inverted (1) [Default: None]
  -v, --verbose  verbose output
  -r SNP         Remove positions with SNP-ratio > r [Default: 0.75]
  =============  =============================================================
"""
import argparse
import os
import sys
import math
from stammp.obj import functions

def run():
    parser = argparse.ArgumentParser(description='Takes all PAR-CLIP sites and traverses through given pileups to get strand specific coverage of all given pileups and divides the PAR-CLIP mutations counts by the sum of the coverages.', epilog="contact: torkler@genzentrum.lmu.de")
    parser.add_argument('inputfile', help='PAR-CLIP file *.table')
    parser.add_argument('outputfile', help='Normalized PAR-CLIP file *.table')
    parser.add_argument('rnaseqfiles', help='Comma separated list of pileup files used for normalization [no whitespaces]')
    parser.add_argument('chrnames',  help='Comma separated, ordered list of chrnames [no whitespaces]')
    parser.add_argument('-s', help='Comma sperated list of 0 or 1 indicating which files have to be inverted (1) [Default: None]', dest='switch', default=None)
    parser.add_argument('-v','--verbose', dest='verbose', action="store_true", default=False, help='verbose output')
    parser.add_argument('-r', help='Remove positions with SNP-ratio > r [Default: 0.75]', default=0.75, type=float, dest='snp')
    args = parser.parse_args()
    
    files_rnaseq = args.rnaseqfiles.split(',')
    switch       = [0]*len(files_rnaseq)
    if args.switch != None:
        switchparameters = args.switch.split(',')
        if len(switchparameters) != len(files_rnaseq):
            print('Unequal number of normalization files and switching information')
            sys.exit(0)
        for s in range(len(switchparameters)):
            if switchparameters[s] == '1':
                switch[s] = 1
            else:
                switch[s] = 0
    
    chrnames      = {}
    splitchrnames = args.chrnames.split(',')
    for s in range(len(splitchrnames)):
        chrnames[splitchrnames[s]] = s
    
    try:
        fc_input  = open(args.inputfile, 'r')
    except:
        print('No such input file: '+args.inputfile)
        sys.exit(0)
    
    try:
        fc_output = open(args.outputfile, 'w')
    except:
        print('Cannot open file connection to: '+args.outputfile)
        sys.exit(0)
    
    fc_rnaseq = {}
    fnames    = []
    print('RNA-seq data using for normalization:')
    for f in range(len(files_rnaseq)):
        fname = files_rnaseq[f].split('/')
        fname = fname[len(fname)-1].split('.')[0]
        try:
            fc_rnaseq[files_rnaseq[f]] = open(files_rnaseq[f], 'r')
            print('\t'+files_rnaseq[f]+'  s = '+str(switch[f]))
        except:
            print('No such RNA-seq pileup file: '+files_rnaseq[f])
            sys.exit(0)
        fnames.append(fname)
    
    
    line_input = fc_input.readline()
    line_input = fc_input.readline()
    fc_output.write('chromosome\tposition\tm\tr\tpvalue\tstrand\tocc\n')
    while(line_input):
        split_input = line_input.split('\t')
        total_rna = total_mut = count_SNP = 0
        #print(split_input[0]+'\t'+split_input[1]+'\t'+split_input[5]+'\t'+split_input[2])
        for f in range(len(fnames)):
            line_rna = fc_rnaseq[files_rnaseq[f]].readline()
            while(line_rna):
                split_rna = line_rna.split('\t')
                if split_rna[0] == split_input[0]:
                    if int(split_rna[1]) == int(split_input[1]):
                        tmp_counts_forward = functions.getCounts(split_rna[4], forward=True)
                        tmp_counts_reverse = functions.getCounts(split_rna[4], forward=False)
                        if split_input[5] == '+':
                            if switch[f] == 0:
                                #print('\t'+files_rnaseq[f]+'\t'+split_rna[0]+'\t'+str(tmp_counts_forward[0])+'\t'+split_input[1])
                                total_rna += tmp_counts_forward[0]
                                #total_mut += tmp_counts_forward[1]['C']
                                total_mut += tmp_counts_forward[2]
                            else:
                                #print('\t'+files_rnaseq[f]+'\t'+split_rna[0]+'\t'+str(tmp_counts_reverse[0])+'\t'+split_input[1])
                                total_rna += tmp_counts_reverse[0]
                                #total_mut += tmp_counts_reverse[1]['G']
                                total_mut += tmp_counts_reverse[2]
                            break
                        else:
                            if switch[f] == 0:
                                #print('\t'+files_rnaseq[f]+'\t'+split_rna[0]+'\t'+str(tmp_counts_reverse[0])+'\t'+split_input[1])
                                total_rna += tmp_counts_reverse[0]
                                #total_mut += tmp_counts_reverse[1]['G']
                                total_mut += tmp_counts_reverse[2]
                            else:
                                #print('\t'+files_rnaseq[f]+'\t'+split_rna[0]+'\t'+str(tmp_counts_forward[0])+'\t'+split_input[1])
                                total_rna += tmp_counts_forward[0]
                                #total_mut += tmp_counts_forward[1]['C']
                                total_mut += tmp_counts_forward[2]
                            break
                    if int(split_rna[1]) >= int(split_input[1]):
                        break
                else:
                    if chrnames[split_input[0]] < chrnames[split_rna[0]]:
                        break
                line_rna = fc_rnaseq[files_rnaseq[f]].readline()
        if total_rna > 1:
            occ = float(split_input[2])/total_rna
        else:
            occ = float(split_input[2])
        #print('  totalRNA '+str(total_rna))
        #print('       occ '+str(occ))
        if total_mut/(total_rna+1) < args.snp:
            fc_output.write(split_input[0]+'\t'+split_input[1]+'\t'+split_input[2]+'\t'+split_input[3]+'\t'+split_input[4]+'\t'+split_input[5]+'\t'+str(occ)+'\n')
            count_SNP += 1
        else:
            count_SNP += 1
        #fc_output.write(split_input[0]+'\t'+split_input[1]+'\t'+split_input[2]+'\t'+split_input[3]+'\t'+split_input[4]+'\t'+split_input[5]+'\t'+str(occ)+'\n')
        line_input = fc_input.readline()
    print(str(count_SNP)+' SNPs removed')
    fc_input.close()
    fc_output.close()

if __name__ == '__main__':
    run()








