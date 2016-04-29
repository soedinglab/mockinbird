"""
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

Example::

    $ stammp-normalize /path/parclipsites.table /path/parclipsites_normalized.table RNAseq1.pileup,...,RNAseqN.pileup chr1,...,chrN -s 0,0

"""
import argparse
from collections import defaultdict, Counter

from stammp.utils.argparse_helper import file_r, file_rw_or_dir_rwx


def run():
    parser = argparse.ArgumentParser(
        description=(
            'Takes all PAR-CLIP sites and traverses through given pileups to get strand '
            'specific coverage of all given pileups and divides the PAR-CLIP mutations '
            'counts by the sum of the coverages.'
        )
    )
    parser.add_argument('input_file', help='PAR-CLIP file *.table', type=file_r)
    parser.add_argument('output_file', help='Normalized PAR-CLIP file *.table',
                        type=file_rw_or_dir_rwx)
    parser.add_argument('RNAseq_pileup', type=file_r)
    rm_snp_help = 'Remove positions with SNP-ratio > r [Default: 0.75]'
    parser.add_argument('--mut_snp_ratio', '-r', help=rm_snp_help, default=0.75, type=float)
    args = parser.parse_args()

    # first pass: mark pc_sites
    pc_sites = defaultdict(dict)
    with open(args.input_file) as site_table:
        site_table.readline()
        for line in site_table:
            chrom, pos, _, _, _, strand, _ = line.split()
            pc_sites[chrom][pos] = strand

    # count coverages
    with open(args.RNAseq_pileup) as mpileup:
        for line in mpileup:
            try:
                chrom, pos, _, _, cov_str, _ = line.split()
            except ValueError:
                continue

            if pos in pc_sites[chrom]:
                cov_toks = Counter(cov_str)
                strand = pc_sites[chrom][pos]
                if strand == '+':
                    cov = cov_toks['.']
                    mut = cov_toks['A'] + cov_toks['C'] + cov_toks['G'] + cov_toks['T']
                else:
                    cov = cov_toks[',']
                    mut = cov_toks['a'] + cov_toks['c'] + cov_toks['g'] + cov_toks['t']
                pc_sites[chrom][pos] = (cov, mut)

    # second pass
    with open(args.input_file) as site_table, open(args.output_file, 'w') as out_file:
        snp_count = 0
        header = site_table.readline()
        print(*header.split(), sep='\t', file=out_file)
        for line in site_table:
            toks = line.split()
            chrom, pos, m, r, pval, strand, occ = toks
            cov_mut = pc_sites[chrom][pos]
            if isinstance(cov_mut, str):
                mut = 0
                total_cov = 0
                occ = m
            else:
                cov, mut = cov_mut
                total_cov = cov + mut
                if total_cov == 0:
                    occ = m
                else:
                    occ = float(m) / (total_cov)

            # update the occupancy
            toks[-1] = occ

            # not clear how to interpret the +1
            if mut / (total_cov + 1) < args.mut_snp_ratio:
                print(*toks, sep='\t', file=out_file)
            else:
                snp_count += 1

    print('Removed %s snps' % snp_count)


if __name__ == '__main__':
    run()
