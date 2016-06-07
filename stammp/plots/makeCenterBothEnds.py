"""
Plot PAR-CLIP data in sense and anti-sense direction around start and stop
positions given in the GFF file. Outputfilename is constructed from
parameters: outputdir/prefix+u+g+d+min+max+labelCenterA+labelBody+labelCenterB
+plotSmooth.pdf

Blue color represents PAR-CLIP signal on the sense strand, green color represents PAR-CLIP signal on the anti-sense strand.

**Usage:** stammp-makeCenterBothEnds [-h] [-d DOWNSTREAM] [-u UPSTREAM] [-g GENE]
                                 [--min MIN] [--max MAX]
                                 [--plotSmooth PLOTSMOOTH]
                                 [--labelCenterA LABELCENTERA]
                                 [--labelBody LABELBODY]
                                 [--labelCenterB LABELCENTERB] [-r]
                                 parclip outputdir prefix gff

**Positional arguments:**
  =========             ==============================
  parclip               path to the PAR-CLIP \*.pileup
  outputdir             output directory
  prefix                prefix of filenames
  gff                   GFF file used for plotting
  =========             ==============================

**Optional arguments:**
  ==================    ======================================================
  -h, --help            show this help message and exit
  -d 1000               set downstream range [default: 1000nt]
  -u 1000               set upstream range [default: 1000nt]
  -g 750                set gene range [default: 750nt]
  --min 0               minium transcript size [default: 0nt]
  --max 5000            maximum transcript size [default: 5000nt]
  --plotSmooth 20       half of the window size used for the running mean
                        [default: 20nt]
  --labelCenterA TSS    plot label for the first center position [default:
                        TSS]
  --labelBody gene      plot label for body (between A and B) [default: gene]
  --labelCenterB pA     plot label for the second center position [default:
                        pA]
  -r, --remove          remove temporary text files. [default: false]
  ==================    ======================================================

Example::

    stammp-makeCenterBothEnds parclip.table outputdirectory/ annotation.gff -d 1000 -u 1000 -g 750 --min 1500 --max 4000 --plotSmooth 20 --labelCenterA TSS --labelBody Gene --labelCenterB pA


.. image:: img/img_plotCenterBoth.png
    :align: center
    :height: 689px
    :alt: alternate text
"""
import argparse
import os
from itertools import chain
import datetime

from stammp.obj import functions, gff, parclipsites
from stammp.utils import execute
from stammp.utils.argparse_helper import file_r, dir_rwx


def main(parclipfile, outputfile, gfffile, downstream, upstream, gene, sense, minSize,
         maxSize, verbose, vstring=''):
    anno = gff.GFF(gfffile)
    anno.filterSize(minSize, maxSize)
    pc = parclipsites.ParclipSites('')
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
    parser = argparse.ArgumentParser(
        description=(
            'Plot PAR-CLIP data in sense and anti-sense direction around start '
            'and stop positions given in the GFF file. Outputfilename is constructed '
            'from parameters: outputdir/prefix+u+g+d+min+max+labelCenterA+labelBody+'
            'labelCenterB+plotSmooth.pdf'
        )
    )
    parser.add_argument('parclip', help='path to the PAR-CLIP *.table', type=file_r)
    parser.add_argument('outputdir', help='output directory', type=dir_rwx)
    parser.add_argument('prefix', help='prefix of filenames')
    parser.add_argument('gff', help='GFF file used for plotting', type=file_r)
    parser.add_argument('-d', help='set downstream range [default: 1000nt]',
                        dest='downstream', default=1000, type=int)
    parser.add_argument('-u', help='set upstream range [default: 1000nt]',
                        dest='upstream', default=1000, type=int)
    parser.add_argument('-g', help='set gene range [default: 750nt]', dest='gene',
                        default=750, type=int)
    parser.add_argument('--min', help='minimum transcript size [default: 0nt]',
                        default=0, type=int)
    parser.add_argument('--max', help='maximum transcript size [default: 5000nt]',
                        default=5000, type=int)
    smooth_help = 'half of the window size used for the running mean [default: 20nt]'
    parser.add_argument('--plotSmooth', help=smooth_help, default=20, type=int)
    label_cenA_help = 'plot label for the first center position [default: TSS]'
    parser.add_argument('--labelCenterA', help=label_cenA_help, default='TSS')
    parser.add_argument('--labelBody', help='for body (between A and B) [default: gene]',
                        default='gene')
    label_cenB_help = 'plot label for the second center position [default: pA]'
    parser.add_argument('--labelCenterB', help=label_cenB_help, default='pA')
    parser.add_argument('--title')
    parser.add_argument('-r', '--remove', action='store_true',
                        help='remove temporary files. [default: false]')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose output')
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
        '%r' % args.title,
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
