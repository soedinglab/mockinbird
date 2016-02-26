#! /usr/bin/python3
import argparse
import os
import sys
import re

#################################################################
######## Functions ##############################################
#################################################################

def keepReadStrict(read, adapter, seed, barcode):
    #print(read[barcode:(barcode+len(adapter))])
    #print(adapter)
    if read[barcode:(barcode+len(adapter))]==adapter:
        return read[(barcode+len(adapter)):(len(read)+1)]
    return read

def keepRead(read, adapter, seed, barcode):
    # clip reads that start with BARCODE nts + the last s bases of the adapter (with s >= seed)
    for s in range(len(adapter),seed-1,-1):
        #print(s)
        #print(read[barcode:(s+barcode)])
        #print(adapter[(len(adapter)-s):(len(adapter)+1)])
        if read[barcode:(s+barcode)]==adapter[(len(adapter)-s):(len(adapter)+1)]:
            return read[(s+barcode):(len(read)+1)]
    return read

def keepReadClipAnywhere(read, adapter, seed, barcode):
    # clip reads that start with any number of nts + the last s bases of the adapter (with s >= seed)
    for s in range(len(adapter),seed-1,-1):
        pattern = adapter[(len(adapter)-s):(len(adapter)+1)]
        match = re.search(pattern[::-1],read[::-1])
        if match!=None:
            #print(match.group(0))
            b=(len(read)-match.end())
            return read[(s+b):(len(read)+1)]
        # for b in range(len(read)-s,-1,-1):
        #     #print(b,s)
        #     #print(read[b:(s+b)])
        #     #print(adapter[(len(adapter)-s):(len(adapter)+1)])
        #     if read[b:(s+b)]==pattern:
        #         #print("b",b)
        #         return read[(s+b):(len(read)+1)]
    return read
        

def removeAdapter(inputfile, outputfile, adapter, strict, seed, barcode, clipanywhere):
    if strict:
        keepReadFunction = keepReadStrict
    else:
        if not clipanywhere:
            keepReadFunction = keepRead
        else:
            keepReadFunction = keepReadClipAnywhere

    with open(inputfile,"r") as f:
        with open(outputfile,"w") as o:
            for ite,l in enumerate(f):
                if ite%4==0:
                    # if ite%400000 == 0:
                    # print('Processed reads %i'%(ite/4))
                    currlines = l
                    keep = True
                elif ite%4==1:
                    # sequence line
                    # decide what to keep from this read
                    keep = keepReadFunction(l.strip(), adapter, seed, barcode)
                    currlines += keep+"\n"
                elif ite%4==2:
                    currlines += l
                elif ite%4==3:
                    if len(keep)>0: # write only non-empty clipped reads
                        qualkeep = l.strip()
                        qualkeep = qualkeep[(len(qualkeep)-len(keep)):(len(qualkeep)+1)]
                        currlines += qualkeep+"\n"
                        o.write(currlines)

def run():
    parser = argparse.ArgumentParser(description='remove 5prime adapters from fastq file', epilog="contact: jessica.andreani@mpibpc.mpg.de")
    parser.add_argument('inputfile',     help='Path to the Inputfile *.fastqsanger')
    parser.add_argument('outputfile',     help='Path to the Outputfile (fastqsanger format)')
    parser.add_argument('--adapter', help='Adapter sequence [default: %(default)s]', default="GTTCAGAGTTCTACAGTCCGACGATC", type=str)
    parser.add_argument('--strict', dest='strict', action="store_true", default=False, help='Strict removal: clip only reads starting with BARCODE nts + the full-length adapter [default: %(default)s]')
    parser.add_argument('--seed', help='Clip all reads with at least SEED bases from the adapter sequence [default: %(default)s]', default=8, type=int)
    parser.add_argument('--barcode', help='Clip reads starting with at least BARCODE nts assumed to be from the random barcode [default: %(default)s]', default=5, type=int)
    parser.add_argument('--clipanywhere',  dest='clipanywhere', action="store_true", default=False, help='Clip the read up to the longest match to the adapter sequence located anywhere in the read [default: %(default)s]')
    args = parser.parse_args()
    
    #if args.strict:
    #    print('Warning, --strict option specified, only reads starting with BARCODE nts + the full-length adapter will be clipped (the --seed and --clipanywhere options do not apply)')
    
    #if not args.strict and args.clipanywhere:
    #    print('Warning, reads will be cleaved up to a (partial) adapter sequence located anywhere (the --barcode option does not apply)')
    
    if os.path.isfile(args.inputfile) == False:
        #Check if necessary files like input files exist
        print('Inputfile: '+args.inputfile+' does not exist')
        sys.exit(0)
    
    removeAdapter(args.inputfile, args.outputfile, args.adapter, args.strict, args.seed, args.barcode, args.clipanywhere)

if __name__ == '__main__':
    run()

