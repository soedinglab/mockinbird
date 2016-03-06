#! /usr/bin/python3
import argparse
import os
import sys
import collections
import math

from stammp.utils import native_wordcount as wccount


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
    parser = argparse.ArgumentParser(description='', epilog="contact: jessica.andreani@mpibpc.mpg.de")
    parser.add_argument('inputfile',     help='Input file (trimmed/quality-filtered FastQ).')
    parser.add_argument('outputfile',     help='Output file (FastQ).')
    parser.add_argument('-v','--verbose', dest='verbose', action="store_true", default=False, help='verbose output')
    args = parser.parse_args()

    if not os.path.exists(args.inputfile):
        sys.exit("Input file does not exist! %s"%(args.inputfile))
    print("### %s ###"%os.path.basename(args.inputfile))
    nblines = wccount(args.inputfile)
    print("Input file has %i lines, %i reads"%(nblines,nblines/4))
    readFastq(args.inputfile, args.verbose, nblines, args.outputfile)

if __name__ == '__main__':
    run()
