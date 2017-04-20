import argparse
import os
from itertools import chain
import datetime

from mockinbird.obj import functions, gff
from mockinbird.utils import ParclipSiteContainer
from mockinbird.utils import execute
from mockinbird.utils.argparse_helper import file_r, dir_rwx


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
    parser.add_argument('parclip', help='path to the PAR-CLIP *.table', type=file_r)
    parser.add_argument('outputdir', help='output directory', type=dir_rwx)
    parser.add_argument('prefix', help='prefix of filenames')
    parser.add_argument('gff', help='GFF file used for plotting', type=file_r)
    parser.add_argument('--downstream', '-d', help='set downstream range',
                        default=1000, type=int)
    parser.add_argument('--upstream', '-u', help='set upstream range',
                        default=1000, type=int)
    parser.add_argument('--gene', '-g', help='set gene range',
                        default=750, type=int)
    parser.add_argument('--min', help='minimum transcript size',
                        default=0, type=int)
    parser.add_argument('--max', help='maximum transcript size',
                        default=5000, type=int)
    smooth_help = 'half of the window size used for the running mean'
    parser.add_argument('--plotSmooth', help=smooth_help, default=20, type=int)
    label_cenA_help = 'plot label for the first center position'
    parser.add_argument('--labelCenterA', help=label_cenA_help, default='TSS')
    parser.add_argument('--labelBody', help='for body (between A and B)',
                        default='gene')
    label_cenB_help = 'plot label for the second center position'
    parser.add_argument('--labelCenterB', help=label_cenB_help, default='pA')
    parser.add_argument('--title', help='plot title')
    parser.add_argument('--remove', '-r', action='store_true',
                        help='remove temporary files')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='verbose output')
    return parser


def main(parclipfile, outputfile, gfffile, downstream, upstream, gene, sense, minSize,
         maxSize, verbose, vstring=''):
    anno = gff.GFF(gfffile)
    anno.filterSize(minSize, maxSize)
    pc = ParclipSiteContainer()
    pc.loadFromFile(parclipfile)
    with open(outputfile, 'w') as fc_out:
        for g in range(anno.size()):
            if verbose:
                functions.showProgress(g, (anno.size() - 1), vstring)
            if anno.strand[g] == '+':
                values_upstream = pc.getValues(anno.chr[g], anno.start[g], anno.strand[g],
                                               sense, upstream, gene)
                values_dostream = pc.getValues(anno.chr[g], anno.stop[g], anno.strand[g],
                                               sense, gene, downstream)
            else:
                values_upstream = pc.getValues(anno.chr[g], anno.stop[g], anno.strand[g],
                                               sense, upstream, gene)
                values_dostream = pc.getValues(anno.chr[g], anno.start[g], anno.strand[g],
                                               sense, gene, downstream)
            if values_upstream is not None and values_dostream is not None:
                print(*chain(values_upstream, values_dostream), sep='\t', file=fc_out)
        if verbose:
            print()


def run():
    scriptPath = os.path.dirname(os.path.realpath(__file__))
    plot_script = os.path.join(scriptPath, 'plotCenterBothEnds.R')
    parser = create_parser()
    args = parser.parse_args()

    if not args.title:
        title = args.prefix
    else:
        title = args.title

    prefix_fmt = '%s_centerBoth_up%s_gene%s_do%s_min%s_max%s'
    fmt_args = (
        args.prefix,
        args.upstream,
        args.gene,
        args.downstream,
        args.min,
        args.max,
    )
    fn_prefix = prefix_fmt % fmt_args
    outfile_sense = os.path.join(args.outputdir, fn_prefix + '_sense.table')
    outfile_asense = os.path.join(args.outputdir, fn_prefix + '_asense.table')
    outfile_pdf = os.path.join(args.outputdir, fn_prefix + '_sm%s.pdf' % args.plotSmooth)

    start_time = datetime.datetime.now()

    main(args.parclip, outfile_sense, args.gff, args.downstream, args.upstream,
         args.gene, True, args.min, args.max, args.verbose,
         'Collecting PAR-CLIP sense data')
    main(args.parclip, outfile_asense, args.gff, args.downstream, args.upstream,
         args.gene, False, args.min, args.max, args.verbose,
         'Collecting PAR-CLIP anti-sense data')

    if args.verbose:
        end_time = datetime.datetime.now()
        run_time = end_time - start_time
        print()
        print('time: %s seconds' % run_time.seconds)

    plot_cmd = [
        'R',
        '-q',
        '--slave',
        '-f %r' % plot_script,
        '--args',
        '%r' % outfile_sense,
        '%r' % outfile_asense,
        '%r' % outfile_pdf,
        '%r' % title,
        args.upstream,
        args.downstream,
        args.gene,
        args.plotSmooth,
        '%r' % args.labelCenterA,
        '%r' % args.labelBody,
        '%r' % args.labelCenterB,
    ]
    execute(plot_cmd)
    if args.remove:
        os.remove(outfile_sense)
        os.remove(outfile_asense)


if __name__ == '__main__':
    run()
