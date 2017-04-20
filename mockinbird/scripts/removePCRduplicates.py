import argparse
import os
import sys
import collections
import math

from stammp.utils import native_wordcount as wccount
from stammp.utils import argparse_helper as aph


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_fastq', type=aph.file_r,
                        help='path to input fastq file')
    parser.add_argument('output_fastq', type=aph.file_rw_or_dir_rwx,
                        help='path to output fastq file.')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='verbose output')
    return parser


def showProgress(count, total, message):
    frac = math.trunc((count/total)*10)
    bar = '\t0% ['
    for i in range(1,11):
        if i <= frac:
            bar = bar+'='
        else:
            bar = bar+'.'
    frac = round((count/total)*100, 2)
    bar = bar+'] 100%\t'+str(frac)+'%\t'+message+'\r'
    sys.stdout.write(bar)

def readFastq(inputfile, verbose, nblines, outfile):
    #readcount = 0
    index = 0
    info = ""
    idxline = 0
    dic = {}
    with open(inputfile,"r") as f:
        with open(outfile,"w") as o:
            for l in f:
                idxline += 1
                if index == 4:# and l[0]=="@":
                    # We've seen the last line for this read
                    if seq in dic:
                        # duplicate
                        dic[seq] += 1
                    else:
                        # first time this read is seen
                        dic[seq] = 1
                        o.write(info)
                    # Start new read
                    #readcount += 1
                    info = l
                    index = 1
                else:
                    index += 1
                    info += l
                    if index==2:
                        seq = l.strip()
                if verbose and idxline%100000==0:
                    showProgress(idxline, nblines, "Reading FastQ")
            # Process the last read
            if index == 4:# and l[0]=="@":
                # We've seen the last line for this read
                if seq in dic:
                    # duplicate
                    dic[seq] += 1
                else:
                    # first time this read is seen
                    dic[seq] = 1
                    o.write(info)
            
    if verbose:
        showProgress(idxline, nblines, "Reading FastQ")
        print("")
    print("There are %i non-redundant reads (%.2f%%)"%(len(dic),len(dic)/nblines*4*100))
    counter = collections.defaultdict(int)
    for k in dic:
        counter[dic[k]]+=1
    # Write duplicate distribution to file
    with open(outfile+".hist","w") as o:
        for i in sorted(counter):
            o.write("%15i  %15i\n"%(i,counter[i]))

def run():
    parser = create_parser()
    args = parser.parse_args()

    print("### %s ###"%os.path.basename(args.input_fastq))
    nblines = wccount(args.input_fastq)
    print("Input file has %i lines, %i reads"%(nblines,nblines/4))
    readFastq(args.input_fastq, args.verbose, nblines, args.output_fastq)

if __name__ == '__main__':
    run()
