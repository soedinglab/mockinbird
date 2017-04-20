import argparse
import os
from stammp.obj import *
from stammp.utils import argparse_helper as aph
from stammp.utils import ParclipSiteContainer


def create_parser():

    parser = argparse.ArgumentParser(
        description='Plots Jaccard-index for two given PAR-CLIP data files.'
    )
    parser.add_argument('parclipA', help='path to the PAR-CLIP *.table', type=aph.file_r)
    parser.add_argument('parclipB', help='path to the PAR-CLIP *.table', type=aph.file_r)
    parser.add_argument('outdir', help='output directory', type=aph.dir_rwx_create)
    parser.add_argument('nameA', help='factor name A')
    parser.add_argument('nameB', help='factor name B')
    parser.add_argument('--zlimit', help='set maximum jaccard index for plotting',
                        default=1.0, type=float)
    parser.add_argument('--width', help='nt +/- parclip site',
                        default=5, type=int)
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='verbose output')
    return parser


def getEntries(sites, q1, q2):
    i = 0
    lower = functions.getQuantile(sites.occ,q1)
    upper = functions.getQuantile(sites.occ,q2)
    pc = ParclipSiteContainer()
    count = 0
    size = sites.size()
    for i in range(sites.size()):
        if sites.occ[i] > lower and sites.occ[i] <= upper:
            count += 1
            pc.addSite(sites.chrs[i], sites.pos[i], sites.m[i], sites.r[i], sites.result[i], sites.strand[i], sites.occ[i])
    pc.getChromosomePositions()
    return pc

def main(parclipA, parclipB, outfile, width, verbose):
    quantiles = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    total       = (len(quantiles)-1)*(len(quantiles)-1)
    total_count = 0
    if verbose:
       functions.showProgress(total_count, total, 'Calculating Jaccard-Index')
    fc = open(outfile, 'w')
    for q1 in range(len(quantiles)-1):
        a = ParclipSiteContainer()
        a.loadFromFile(parclipA)
        aq = getEntries(a,quantiles[q1], quantiles[q1+1])
        #removeEntries(a,quantiles[q1], quantiles[q1+1])
        for q2 in range(len(quantiles)-1):
            b = ParclipSiteContainer()
            b.loadFromFile(parclipB)
            #removeEntries(b,quantiles[q2], quantiles[q2+1])
            bq = getEntries(b,quantiles[q2], quantiles[q2+1])
            intersect = 0
            for j in range(bq.size()):
                if aq.exactSearch(bq.chrs[j], bq.pos[j], bq.strand[j], width=width)[1]:
                    intersect += 1
            jaccard = intersect/(aq.size()+bq.size()-intersect)
            #print('q1: '+str(quantiles[q1])+' q2: '+str(quantiles[q2])+' '+str(round(jaccard,4)))
            fc.write(str(round(jaccard,4))+'\t')
            total_count += 1
            if verbose:
                functions.showProgress(total_count, total, 'Calculating Jaccard-Index')
        fc.write('\n')
    print('')
    fc.close()

def run():
    scriptPath = os.path.dirname(os.path.realpath(__file__))
    scriptPath = scriptPath+'/'
    parser = create_parser()
    args = parser.parse_args()
    outfile = os.path.abspath(args.outdir)
    outfile = outfile+'/jaccard_'+args.nameA+'_'+args.nameB+'.table'
    main(args.parclipA, args.parclipB, outfile, args.width, args.verbose)
    outfile_pdf = os.path.abspath(args.outdir)
    outfile_pdf = outfile_pdf+'/jaccard_'+args.nameA+'_'+args.nameB+'_zlim_'+str(args.zlimit)+'.pdf'
    os.system('R -q --slave -f '+scriptPath+'plotJaccard.R --args ' \
    + outfile + ' ' + str(args.zlimit) + ' ' +outfile_pdf +' ' \
    +args.nameA+' '+args.nameB)

if __name__ == '__main__':
    run()
