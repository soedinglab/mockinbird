#! /usr/bin/python3
"""
Plot PAR-CLIP data in senseand anti-sense as heat map from given annotation
from a GFF file

*Usage:* stammp-makeHeatMap [-h] [-d DOWNSTREAM] [-u UPSTREAM] [--min MIN]
                            [--max MAX] [--xsize XSIZE] [--ysize YSIZE] [-r]
                            [-v]
                            parclip outputdir prefix gff

*Positional arguments:*
  =========      ============================
  parclip        path to the PAR-CLIP *.table
  outputdir      output directory
  prefix         prefix of filenames
  gff            GFF file used for plotting
  =========      ============================

*Optional arguments:*
  =============  ==============================================
  -h, --help     show this help message and exit
  -d DOWNSTREAM  set downstream range [default: 1000nt]
  -u UPSTREAM    set upstream range [default: 4000nt]
  --min MIN      minium transcript size [default: 0nt]
  --max MAX      maximum transcript size [default: 5000nt]
  --xsize XSIZE  number of bins in x direction [default: 500]
  --ysize YSIZE  number of bins in y direction [default: 500]
  -r, --remove   remove temporary text files. [default: false]
  -v, --verbose  verbose output
  =============  ==============================================

.. image:: tut_pub1_heatmap_sense.png
    :align: center
    :height: 500px
    :alt: alternate text
    
.. image:: tut_pub1_heatmap_asense.png
    :align: center
    :height: 500px
    :alt: alternate text
"""
import argparse
import os
from stammp.obj import *

def main(parclipfile, gfffile, upstream, downstream, sense, minSize, 
         maxSize, verbose, xbins, ybins, vstring=''):
    anno = gff.GFF(gfffile)
    anno.filterSize(minSize, maxSize)
    anno.sort2size()
    pc = parclipsites.ParclipSites('')
    pc.loadFromFile(parclipfile)
    mat = []
    annosize = []
    for g in range(anno.size()):
        if verbose:
            functions.showProgress(g, (anno.size()-1), vstring)
        if anno.strand[g] == '+':
            values = pc.getValues(anno.chr[g], anno.start[g], anno.strand[g], 
                                  sense, upstream, downstream)
        else:
            values = pc.getValues(anno.chr[g], anno.stop[g], anno.strand[g], 
                                  sense, upstream, downstream)
        if values != None:
            mat.append(functions.shrinkValues(values, xbins))
            annosize.append(anno.stop[g]-anno.start[g])
    smat = []
    sannosize = []
    ystep = round(anno.size()/ybins)
    ystart = 0
    ystop = ystep
    while ystop < anno.size():
        tmp = [0]*xbins
        for i in range(xbins):
            count = 0
            tmpanno = 0
            for j in range(ystart,ystop):
                tmp[i] += mat[j][i] #[row][col]
                tmpanno += annosize[j]
                count += 1
            tmp[i] = tmp[i]/count
            tmpanno = tmpanno/count
        smat.append(tmp)
        sannosize.append(tmpanno)
        ystart = ystop + 1
        ystop += ystep
    return(smat, sannosize)
    if verbose:
        print('')

def saveMat(outfile, mat, upstream, downstream, annosize, xbins):
    fc = open(outfile, 'w')
    for i in range(len(mat)):
        fc.write(str(((upstream+annosize[i]) / (upstream+downstream)) * xbins)+'\t')
        for j in range(len(mat[i])):
            fc.write(str(mat[i][j])+'\t')
        fc.write('\n')
    fc.close()

def run():
    scriptPath = os.path.dirname(os.path.realpath(__file__))
    scriptPath = scriptPath+'/'
    parser = argparse.ArgumentParser(description='Plot PAR-CLIP data in sense'
    + ' and anti-sense as heat map from given annotation from a GFF file', 
    epilog="contact: torkler@genzentrum.lmu.de")
    parser.add_argument('parclip', help='path to the PAR-CLIP *.table')
    parser.add_argument('outputdir', help='output directory')
    parser.add_argument('prefix', help='prefix of filenames')
    parser.add_argument('gff', help='GFF file used for plotting')
    parser.add_argument('-d', help='set downstream range [default: 1000nt]', 
                        dest='downstream', default=1000, type=int)
    parser.add_argument('-u', help='set upstream range [default: 4000nt]', 
                        dest='upstream', default=4000, type=int)
    parser.add_argument('--min', help='minium transcript size [default: 0nt]',  
                        default=0, type=int)
    parser.add_argument('--max', help='maximum transcript size [default: 5000nt]', 
                        default=5000, type=int)
    parser.add_argument('--xsize', help='number of bins in x direction [default: 500]',  
                        default=500, type=int)
    parser.add_argument('--ysize', help='number of bins in y direction [default: 500]',  
                        default=500, type=int)
    parser.add_argument('-r','--remove', dest='remove', action="store_true", 
                        default=False, help='remove temporary text files. [default: false]')
    parser.add_argument('-v','--verbose', dest='verbose', action="store_true", 
                        default=False, help='verbose output')
    args = parser.parse_args()    

    functions.checkExistence(args.parclip)
    functions.checkExistence(args.gff)
#    functions.checkPath(args.outputdir)
    
    outfile_mat_sense = args.outputdir + args.prefix+'_mat_up' \
    + str(args.upstream) + '_do' + str(args.downstream) + '_min' \
    + str(args.min) + '_max' + str(args.max) + '_sense.table'

    outfile_mat_asense = args.outputdir + args.prefix + '_mat_up' \
    + str(args.upstream) + '_do' + str(args.downstream) + '_min' \
    + str(args.min) + '_max' + str(args.max) + '_asense.table'

    outfile_img_sense = args.outputdir + args.prefix + '_mat_up' \
    + str(args.upstream) + '_do' + str(args.downstream) + '_min' \
    + str(args.min) + '_max' + str(args.max) + '_sense.png'

    outfile_img_asense = args.outputdir + args.prefix + '_mat_up' \
    + str(args.upstream) + '_do' + str(args.downstream) + '_min' \
    + str(args.min) + '_max' + str(args.max) + '_asense.png'

    sense = main(args.parclip, args.gff, args.upstream, args.downstream, True, 
                 args.min, args.max, args.verbose, args.xsize, args.ysize, 
                 'Collecting data from sense strand')
                 
    asense = main(args.parclip, args.gff, args.upstream, args.downstream, False,
                  args.min, args.max, args.verbose, args.xsize, args.ysize, 
                  'Collecting data from anti-sense strand')

    saveMat(outfile_mat_sense, sense[0], args.upstream, args.downstream, 
            sense[1], args.xsize)

    saveMat(outfile_mat_asense, asense[0], args.upstream, args.downstream, 
            asense[1], args.xsize)

    os.system('R -q --slave -f ' + scriptPath + 'plotHeatMap.R --args ' \
    + outfile_mat_sense + ' ' + outfile_mat_asense + ' ' + outfile_img_sense \
    + ' ' + outfile_img_asense + ' 0.98 ' \
    + str((args.upstream/(args.upstream+args.downstream)*args.xsize)) + ' ' \
    + str(args.xsize) + ' ' + str(args.ysize))
    
    if args.remove:
        os.remove(outfile_sense)
        os.remove(outfile_asense)

if __name__ == '__main__':
    run()