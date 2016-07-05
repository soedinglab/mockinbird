import argparse
import os
import sys
from stammp.obj import *
from stammp.utils import argparse_helper as aph


def create_parser():
    description = 'Set maximum occupancy to the specified quantile.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('inputfile', help='normalized PAR-CLIP file *.table', type=aph.file_r)
    parser.add_argument('outputfile', help='converted PAR-CLIP file *.table',
                        type=aph.file_rw_or_dir_rwx)
    parser.add_argument('--quantile', '-q', help='quantile [0, 1.0]', default=0.95,
                        type=float)
    return parser


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
    parser = create_parser()
    args = parser.parse_args()
    main(args.inputfile, args.outputfile, args.quantile)


if __name__ == '__main__':
    run()
