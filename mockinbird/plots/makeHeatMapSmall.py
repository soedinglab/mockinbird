import argparse
import os
from mockinbird.utils import argparse_helper as aph
from mockinbird.utils import execute
from mockinbird.obj import gff, functions
from mockinbird.utils import ParclipSiteContainer


def create_parser():
    parser = argparse.ArgumentParser(
        description=(
            'Plot PAR-CLIP data in sense and anti-sense direction as heat maps '
            'based on annotations of a GFF file. The plot is centered at the '
            'start coordinate given in the GFF. The data '
            '[(start-UPSTREAM),(stop+downstream)] is plotted. Note, '
            'that no binning in y-direction is performed if the value of --ybins '
            'is smaller compared to the number of entries in the GFF.'
        )
    )
    parser.add_argument('parclip', help='path to the PAR-CLIP *.table', type=aph.file_r)
    parser.add_argument('outputdir', help='output directory', type=aph.dir_rwx)
    parser.add_argument('prefix', help='prefix of filenames')
    parser.add_argument('gff', help='GFF file used for plotting', type=aph.file_r)
    parser.add_argument('--downstream', '-d', help='set downstream range',
                        default=500, type=int)
    parser.add_argument('--upstream', '-u', help='set upstream range',
                        default=1000, type=int)
    parser.add_argument('--min', help='minium transcript size',
                        default=0, type=int)
    parser.add_argument('--max', help='maximum transcript size',
                        default=5000, type=int)
    parser.add_argument('--xbins', help='number of bins in x direction',
                        default=500, type=int)
    parser.add_argument('--ybins', help='number of bins in y direction',
                        default=500, type=int)
    parser.add_argument('--xpx', help='width of final plot in px',
                        default=500, type=int)
    parser.add_argument('--ypx', help='height of final plot in px',
                        default=500, type=int)
    parser.add_argument('--remove', '-r', action='store_true',
                        help='remove temporary text files')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='verbose output')
    return parser


def main(parclipfile, gfffile, upstream, downstream, sense, minSize,
         maxSize, verbose, xbins, ybins, vstring=''):
    anno = gff.GFF(gfffile)
    anno.filterSize(minSize, maxSize)
    totalsize = upstream + maxSize + 1 + downstream
    anno.sort2size()
    pc = ParclipSiteContainer.from_file(parclipfile)
    mat = []
    annosize = []
    for g in range(anno.size()):
        tmp = [-1] * totalsize
        if verbose:
            functions.showProgress(g, (anno.size() - 1), vstring)
        if anno.strand[g] == '+':
            values = pc.getValues(anno.chr[g], anno.start[g], anno.strand[g],
                                  sense, upstream,
                                  (anno.stop[g] - anno.start[g]) + downstream)
        else:
            values = pc.getValues(anno.chr[g], anno.stop[g], anno.strand[g],
                                  sense, upstream,
                                  (anno.stop[g] - anno.start[g]) + downstream)
        if values is not None:
            tmp[0:(len(values) - 1)] = values
        mat.append(functions.shrinkValues(tmp, xbins))
        annosize.append(anno.stop[g] - anno.start[g])
    smat = []
    sannosize = []
    if ybins >= anno.size():
        print('Warning: --ybins >= entries in ' + gfffile)
        ybins = anno.size()
    ystep = round(anno.size() / ybins)
    ystart = 0
    ystop = ystep
    while ystop < anno.size():
        tmp = [0] * xbins
        for i in range(xbins):
            count = 0
            tmpanno = 0
            for j in range(ystart, ystop):
                tmp[i] += mat[j][i]  # [row][col]
                tmpanno += annosize[j]
                count += 1
            tmp[i] = tmp[i] / count
            tmpanno = tmpanno / count
        smat.append(tmp)
        sannosize.append(tmpanno)
        ystart = ystop
        ystop += ystep
    return smat, sannosize
    if verbose:
        print()


def saveMat(outfile, mat, upstream, downstream, annosize, xbins, totalsize):
    with open(outfile, 'w') as fc:
        for i in range(len(mat)):
            fc.write(str(((upstream + annosize[i]) / totalsize) * xbins) + '\t')
            for j in range(len(mat[i])):
                fc.write(str(mat[i][j]) + '\t')
            fc.write('\n')


def run():
    scriptPath = os.path.dirname(os.path.realpath(__file__))
    plot_script = os.path.join(scriptPath, 'plotHeatMapSmall.R')
    parser = create_parser()
    args = parser.parse_args()

    functions.checkExistence(args.parclip)
    functions.checkExistence(args.gff)

    prefix_pattern = '%s_mat_sm_up%s_do%s_min%s_max%s_xbins%s_ybins%s'
    out_prefix = prefix_pattern % (args.prefix, args.upstream, args.downstream,
                                   args.min, args.max, args.xbins, args.ybins)

    outfile_mat_sense = os.path.join(args.outputdir, out_prefix + '_sense.table')
    outfile_mat_asense = os.path.join(args.outputdir, out_prefix + '_asense.table')
    outfile_img_sense = os.path.join(args.outputdir, out_prefix + '_sense.png')
    outfile_img_asense = os.path.join(args.outputdir, out_prefix + '_asense.png')

    sense = main(args.parclip, args.gff, args.upstream, args.downstream, True,
                 args.min, args.max, args.verbose, args.xbins, args.ybins,
                 'Collecting data from sense strand')

    asense = main(args.parclip, args.gff, args.upstream, args.downstream, False,
                  args.min, args.max, args.verbose, args.xbins, args.ybins,
                  'Collecting data from anti-sense strand')

    total = args.upstream + args.max + 1 + args.downstream
    saveMat(outfile_mat_sense, sense[0], args.upstream, args.downstream,
            sense[1], args.xbins, total)

    saveMat(outfile_mat_asense, asense[0], args.upstream, args.downstream,
            asense[1], args.xbins, total)

    start = args.upstream / total * args.xbins

    cmd = [
        'R',
        '-q',
        '--slave',
        '-f %r' % plot_script,
        '--args',
        '%r' % outfile_mat_sense,
        '%r' % outfile_mat_asense,
        '%r' % outfile_img_sense,
        '%r' % outfile_img_asense,
        0.98,  # hard-coded qvalue
        start,
        args.ypx,
        args.xpx,
    ]
    execute(cmd)

    if args.remove:
        os.remove(outfile_mat_sense)
        os.remove(outfile_mat_asense)


if __name__ == '__main__':
    run()
