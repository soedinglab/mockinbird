import argparse
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd

from stammp.obj import gff, parclipsites
from stammp.utils import execute
from stammp.utils.argparse_helper import file_r, dir_rwx


def create_parser():
    parser = argparse.ArgumentParser(
        description=(
            'Plot PAR-CLIP data in sense and anti-sense direction around start '
            'and stop positions given in the GFF file. The output filename is constructed '
            'from parameters: outputdir/prefix+u+g+d+min+max+labelCenterA+labelBody+'
            'labelCenterB+plotSmooth.pdf\n\n'
            'Blue color represents PAR-CLIP signal on the sense strand, green color '
            'represents PAR-CLIP signal on the anti-sense strand.'
        )
    )
    parser.add_argument('pc_table', help='path to the PAR-CLIP *.table', type=file_r)
    parser.add_argument('outputdir', help='output directory', type=dir_rwx)
    parser.add_argument('prefix', help='prefix of filenames')
    parser.add_argument('gff_file', help='GFF file used for plotting', type=file_r)
    parser.add_argument('--downstream_bp', '-d', help='set downstream range',
                        default=1000, type=int)
    parser.add_argument('--upstream_bp', '-u', help='set upstream range',
                        default=1000, type=int)
    parser.add_argument('--gene_bp', '-g', help='set gene range',
                        default=750, type=int)
    parser.add_argument('--min_ts_len', help='minimum transcript size',
                        default=0, type=int)
    parser.add_argument('--max_ts_len', help='maximum transcript size',
                        default=5000, type=int)
    smooth_help = 'window size used for running mean smoothing'
    parser.add_argument('--smooth_window', help=smooth_help, default=41, type=int)
    label_cenA_help = 'plot label for the first center position'
    parser.add_argument('--labelCenterA', help=label_cenA_help, default='TSS')
    parser.add_argument('--labelBody', help='for body (between A and B)',
                        default='gene')
    label_cenB_help = 'plot label for the second center position'
    parser.add_argument('--labelCenterB', help=label_cenB_help, default='pA')
    parser.add_argument('--title', help='plot title')
    parser.add_argument('--cleanup', action='store_true',
                        help='remove temporary files')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--n_bs_iterations', default=500, type=int)
    parser.add_argument('--n_processes', default=4, type=int)
    return parser

global data_dict

def calc_profile(indices, bp_size):
    avg_vec = np.zeros(bp_size)
    for ind in indices:
        avg_vec += data_dict[ind]
    avg_vec /= len(indices)
    return avg_vec


def main():
    parser = create_parser()
    args = parser.parse_args()

    np.random.seed(args.seed)

    anno = gff.GFF(args.gff_file)
    anno.filterSize(args.min_ts_len, args.max_ts_len)
    pc = parclipsites.ParclipSites()
    pc.loadFromFile(args.pc_table)

    cut_len = args.downstream_bp + args.upstream_bp + 2 * args.gene_bp + 2

    def aggregate_data(g, sense=True):
        if anno.strand[g] == '+':
            values_upstream = pc.getValues(anno.chr[g], anno.start[g], anno.strand[g],
                                           sense, args.upstream_bp, args.gene_bp)
            values_dostream = pc.getValues(anno.chr[g], anno.stop[g], anno.strand[g],
                                           sense, args.gene_bp, args.downstream_bp)
        else:
            values_upstream = pc.getValues(anno.chr[g], anno.stop[g], anno.strand[g],
                                           sense, args.upstream_bp, args.gene_bp)
            values_dostream = pc.getValues(anno.chr[g], anno.start[g], anno.strand[g],
                                           sense, args.gene_bp, args.downstream_bp)
        if values_upstream is not None and values_dostream is not None:
            data = np.hstack((values_upstream, values_dostream))
            smooth_data = pd.Series(data).rolling(window=args.smooth_window, center=True,
                                                  min_periods=0).mean()
            return smooth_data

    def write_out_data(sense, out_file, bs_file):
        global data_dict
        data_dict = {}
        for g in range(anno.size()):
            data = aggregate_data(g, sense)
            data_dict[g] = data
        ids = np.array(list(data_dict.keys()))

        # actual smoothed curve
        with open(out_file, 'w') as out:
            avg_vec = np.zeros(cut_len)
            for ts_ind in ids:
                avg_vec += data_dict[ts_ind]
            avg_vec /= len(ids)
            print(*avg_vec, sep='\t', file=out)


        with open(bs_file, 'w') as out:
            with Pool(args.n_processes) as pool:
                jobs = []
                for bs in range(args.n_bs_iterations):
                    ts_ind = np.random.choice(ids, size=len(ids))
                    job = pool.apply_async(calc_profile, args=(ts_ind, cut_len))
                    jobs.append(job)

                for job in jobs:
                    res_vec = job.get()
                    print(*res_vec, sep='\t', file=out)

    prefix_fmt = '%s_centerBoth_up%s_gene%s_do%s_min%s_max%s'
    fmt_args = (
        args.prefix,
        args.upstream_bp,
        args.gene_bp,
        args.downstream_bp,
        args.min_ts_len,
        args.max_ts_len,
    )
    fn_prefix = prefix_fmt % fmt_args

    sense_table = os.path.join(args.outputdir, fn_prefix + '_sense.table')
    sense_bs_table = os.path.join(args.outputdir, fn_prefix + '_sense_bs.table')

    asense_table = os.path.join(args.outputdir, fn_prefix + '_asense.table')
    asense_bs_table = os.path.join(args.outputdir, fn_prefix + '_asense_bs.table')

    write_out_data(True, sense_table, sense_bs_table)
    write_out_data(False, asense_table, asense_bs_table)

    scriptPath = os.path.dirname(os.path.realpath(__file__))
    plot_script = os.path.join(scriptPath, 'plotCenterBothEnds_bs.R')

    if not args.title:
        title = args.prefix
    else:
        title = args.title

    outfile_pdf = os.path.join(args.outputdir, fn_prefix + '_sm%s.pdf' % args.smooth_window)

    plot_cmd = [
        'R',
        '-q',
        '--slave',
        '-f %r' % plot_script,
        '--args',
        '%r' % sense_table,
        '%r' % sense_bs_table,
        '%r' % asense_table,
        '%r' % asense_bs_table,
        '%r' % outfile_pdf,
        '%r' % title,
        args.upstream_bp,
        args.downstream_bp,
        args.gene_bp,
        args.smooth_window,
        '%r' % args.labelCenterA,
        '%r' % args.labelBody,
        '%r' % args.labelCenterB,
    ]
    execute(plot_cmd)
    if args.cleanup:
        os.remove(sense_table)
        os.remove(sense_bs_table)
        os.remove(asense_table)
        os.remove(asense_bs_table)


if __name__ == '__main__':
    main()
