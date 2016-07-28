import argparse
import os
import sys
import math
import time
import scipy.misc
import scipy.optimize
import scipy.special
import random
from stammp.obj import functions
from stammp.utils import argparse_helper as aph


def create_parser():
    description = ('stammp-bsfinder detects protein RNA binding sites from PAR-CLIP '
                   'experiments in mpileup files.')
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('inputfile', help='path to the inputfile *.pileup', type=aph.file_r)
    parser.add_argument('outputfile', help='define output file *.table')
    parser.add_argument('--threshold', '-p', default=0.005, type=float,
                        help='set maximum p-value for site selection')
    parser.add_argument('--mincov', '-c', help='set minimum coverage', default=3, type=int)
    parser.add_argument('--reference', '-r', help='set default reference nucleotide',
                        default='T')
    parser.add_argument('--mutation', '-m', help='set default mutation nucleotide',
                        default='C')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='verbose output')
    return parser


def processPileup(file_pileup, verbose, reference='T', mutation='C'):
    pseudocount = 1
    maxr        = 1000 #upper limit for sequencing depth which is used for parameter estimation
    reference_reverse   = functions.makeReverseComplement(reference)
    mutation_reverse    = functions.makeReverseComplement(mutation)
    
    lines = 0
    for line in open(file_pileup):
        lines += 1
    
    mr_list_neg    = [[0]] #initializing empty matrices for parameter estimation
    mr_list_signal = [[0]]
    for i in range(1,maxr):
        mr_list_neg.append([pseudocount]*(i+1))
        mr_list_signal.append([pseudocount]*(i+1))
    
    file_pileup = open(file_pileup, 'r')
    line        = file_pileup.readline()
    count       = percent_old = percent_new = 0
    if verbose: functions.showProgress(count,lines,'Processing Pileup')
    while(line):
        count += 1
        split = line.split('\t')
        if split[2] == reference:
            tmp_counts = functions.getCounts(split[4], forward=True)
            counts     = [tmp_counts[0], tmp_counts[1][mutation]]
            if counts[0] < maxr:
                mr_list_signal[counts[0]][counts[1]] += 1
        if split[2] == reference_reverse:
            tmp_counts = functions.getCounts(split[4], forward=False)
            counts     = [tmp_counts[0], tmp_counts[1][mutation_reverse]]
            if counts[0] < maxr:
                mr_list_signal[counts[0]][counts[1]] += 1
        if split[2] != reference and split[2] != reference_reverse:
            tmp_counts = functions.getCounts(split[4], forward=True)
            counts     = [tmp_counts[0], max(tmp_counts[1].values())]
            if counts[0] < maxr:
                mr_list_neg[counts[0]][counts[1]] += 1
        percent_new = math.trunc((count/lines)*100)
        if(percent_new > percent_old):
            if verbose: functions.showProgress(count,lines,'Processing Pileup')
            percent_old = percent_new
        line = file_pileup.readline()
    if verbose: print('')
    file_pileup.close()
    return([mr_list_neg, mr_list_signal])

def logit(message, verbose, logfile):
    if verbose:
        sys.stdout.write(message)
    logfile.write(message)

#################################################################
######## Model fiting functions
#################################################################

def estimateParametersIter(mr_neg, min_r=5, max_r=10, iterations=5):
    est_params = []
    est_fail   = 0
    total_fits = 0
    for i in range(min_r, (max_r+1), 1):
        try:
#                                                                         rho, a0,  a1,w
            bg_params    = scipy.optimize.fmin_bfgs(f=llh_background, x0=[2.94,-4.8,-2,0], args=(mr_neg[i], -1), disp=0, full_output=1)
            total_fits += 1
        except:
            pass
            #print('Unexpected Scipy Error')
        for iter_count in range(iterations):
            try:
                bg_tmp    = scipy.optimize.fmin_bfgs(f=llh_background, x0=[random.uniform(-4,4),random.uniform(-4,4),random.uniform(-4,4),random.uniform(-4,4)], args=(mr_neg[i], -1), disp=0, full_output=1)
                total_fits += 1
                if bg_tmp[1] < bg_params[1] and bg_tmp[6] == 0:
                    bg_params = bg_tmp
            except:
                pass
                #print('Unexpected Scipy Error')
            try:
                if bg_params[6] != 0:
                    est_fail += 1
                else:
                    est_params.append(bg_params[0])
            except:
                est_fail += 1
            
    par = [0,0,0,0]
    if len(est_params) > 0:
        for i in range(len(par)):
            tmp = 0
            for j in range(len(est_params)):
                tmp += est_params[j][i]
            par[i] = tmp/len(est_params)
    else:
        par = None
        sys.exit('\tParameter-estimation failed! Too little information available.')
    #print('Total fits: '+str(total_fits)+' \t total fails: '+str(est_fail))
    return(par)

#def llh_background(params, data, sign=1.0):
#    llh    = 0
#    number = 0
#    for d in data:
#        number      += d[2]
#        llh         += d[2] * math.log( prob_bg(d[1],d[0],params) )
#    return((sign*llh)/number)

def llh_background(params, data, sign=1.0):
    llh    = 0
    number = 0
    for i in range(len(data)):
        number      += data[i]
        llh         += data[i] * math.log( prob_bg(i,(len(data)-1),params) )
    return((sign*llh)/number)

def prob_bg(k,N, par):
    rho    = 1/(1+math.exp(par[0]))
    alpha0 = math.exp(par[1])
    alpha1 = math.exp(par[2])
    r      = 1/(1+math.exp(par[3]))
    
    nck          = scipy.misc.comb(N,k,exact=1)
    betaBinomial = math.exp( math.log(scipy.special.beta(k+alpha0, N-k+alpha1)) - math.log(scipy.special.beta(alpha0, alpha1)) + math.log(nck))
    binomial     = math.exp( math.log(nck) + k * math.log(rho) + (N-k) * math.log(1-rho))
    
    return(r * binomial + (1-r) * betaBinomial)

def precalculateDistributions(minN, maxN, par, verbose):
    distributions = [0]*(maxN+1)
    for i in range(minN,(maxN+1)):
        tmp_probabilities = []
        for j in range((i+1)):
            tmp_probabilities.append(prob_bg(j, i, par))
        distributions[i] = tmp_probabilities
        if verbose: functions.showProgress(i+1, (maxN+1), 'precalculating probability distributions')
    if verbose: print('')
    return(distributions)

def getPvalue(k,N, probabilities, SNPlikely = False):
    pvalue = 0
    
    if SNPlikely == True:
        if k >= N:
            pvalue = 1
        else:
            for p in probabilities[N][k:N]:
                pvalue += p
    else:
        pvalue = sum(probabilities[N][k::])
#        for p in probabilities[N][k:(N+1)]:
#            pvalue += p
    return(pvalue)

def checkSNP(start, N, par):
    p = 0
    for i in range(start,N):
        p += prob_bg(i, N, par)
    if prob_bg(N, N, par) > p:
        return(True)
    else:
        return(False)

def findPvalueParclipInPileup(pileup, outputfile, mincov, maxcov, probabilities, verbose, reference='T', mutation='C', maxPvalue=0.001, SNPlikely = False):
    reference_reverse   = functions.makeReverseComplement(reference)
    mutation_reverse    = functions.makeReverseComplement(mutation)
    
    found_sites = 0
    lines = 0
    for line in open(pileup):
        lines += 1
    file_pileup = open(pileup, 'r')
    file_table  = open(outputfile, 'w')
    file_table.write('chromosome\tposition\tm\tr\tpvalue\tstrand\tocc\n')
    line  = file_pileup.readline()
    linecount = percent_old = percent_new = 0
    if verbose: functions.showProgress(linecount,lines,'Processing Pileup')
    while(line):
        linecount += 1
        split = line.split('\t')
        counts = [0,0]
        pvalue = 1
        if split[2] == reference:
            tmp_counts = functions.getCounts(split[4], forward=True)
            counts     = [tmp_counts[0], tmp_counts[1][mutation]]
            if counts[0] > mincov and counts[1] > 0:
                if counts[0] > 500:
                    pvalue = getPvalue(round((counts[1]/counts[0])*500), 500, probabilities, SNPlikely)
                else:
                    pvalue = getPvalue(counts[1], counts[0], probabilities, SNPlikely)
                if pvalue <= maxPvalue:
                    file_table.write(split[0]+'\t'+split[1]+'\t'+str(counts[1])+'\t'+str(counts[0])+'\t'+str(pvalue)+'\t+\t0\n')
                    found_sites += 1
        if split[2] == reference_reverse:
            tmp_counts = functions.getCounts(split[4], forward=False)
            counts     = [tmp_counts[0], tmp_counts[1][mutation_reverse]]
            if counts[0] > mincov and counts[1] > 0:
                if counts[0] > 500:
                    pvalue = getPvalue(round((counts[1]/counts[0])*500), 500, probabilities, SNPlikely)
                else:
                    pvalue = getPvalue(counts[1], counts[0], probabilities, SNPlikely)
                if pvalue <= maxPvalue:
                    file_table.write(split[0]+'\t'+split[1]+'\t'+str(counts[1])+'\t'+str(counts[0])+'\t'+str(pvalue)+'\t-\t0\n')
                    found_sites += 1
        percent_new = math.trunc((linecount/lines)*100)
        if(percent_new > percent_old):
            if verbose: functions.showProgress(linecount,lines,'Processing Pileup')
            percent_old = percent_new
        line = file_pileup.readline()
    file_table.close()
    file_pileup.close()
    print('Found %s PAR-CLIP sites.' % found_sites)

def main(inputfile, outputfile, threshold, mincov, reference, mutation, verbose):
    """
    Finds PAR-CLIP binding sites in a given pileup file.
    
    Args:
        inputfile (str): Pileup file
        outputfile (str): Output file
        threshold (float): p-value cutoff. Only sites with a p-value < threshold are reported
        mincov (int): Only sites with a minium number of reads are considered
        reference (char): Nucleotide which is supposed to be present without a crosslink [typically a 'T']
        mutation (char): Nucleotide which is the result of the crosslink induced mutation [typically a 'C']
        verbose (bool): Set to False to avoid any messages to stdout.
    
    Can be accessed directly::
        
        $ stammp-bsfinder /path/to/input.mpileup /path/to/output.table
    """
    if os.path.isfile(inputfile) == False:
        print('Inputfile: '+inputfile+' does not exist')
        sys.exit(0)
    
    CONST_OUTDIR = outputfile.rpartition('/')[0]
    CONST_OUTDIR = CONST_OUTDIR+'/'
    
    if os.path.exists(CONST_OUTDIR) == False:
        os.makedirs(CONST_OUTDIR)
        #sys.exit(0)
    
    LOG_3P      = open(outputfile.rpartition('.')[0]+'_logfile.txt', 'w')
    logit('###### '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())+'\n', verbose, LOG_3P)
    logit('Input                 : '+inputfile+'\n', verbose, LOG_3P)
    logit('Output                : '+outputfile+'\n', verbose, LOG_3P)
    logit('max pvalue            : '+str(threshold)+'\n', verbose, LOG_3P)
    logit('min coverage          : '+str(mincov)+'\n', verbose, LOG_3P)
    
    pileupdata    = processPileup(inputfile, verbose, reference=reference, mutation=mutation)
    par           = estimateParametersIter(pileupdata[0], min_r=5, max_r=10, iterations=100)
    logit('Parameter             :', verbose, LOG_3P)
    for p in par:
        logit(str(round(p,6))+', ', verbose, LOG_3P)
    logit('\n', verbose, LOG_3P)
    distributions = precalculateDistributions(mincov, 500, par, verbose)
    SNPlikely     = checkSNP(5, 11, par)
    logit('Searching and saving PAR-CLIP sites...\n', verbose, LOG_3P)
    findPvalueParclipInPileup(inputfile, outputfile, mincov, 10000, distributions, verbose, reference=reference, mutation=mutation, maxPvalue=threshold, SNPlikely = SNPlikely)
    logit('Searching and saving PAR-CLIP sites...\tDone!\n', verbose, LOG_3P)
    LOG_3P.close()

def run():
    parser = create_parser()
    args = parser.parse_args()
    main(args.inputfile, args.outputfile, args.threshold, args.mincov, args.reference, args.mutation, args.verbose)

if __name__ == '__main__':
    run()
