#! /usr/bin/python3
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
import sys
from stammp.obj import *
import datetime

def main(parclipfile, outputfile, gfffile, downstream, upstream, gene, sense, minSize, maxSize, verbose, vstring=''):
    anno = gff.GFF(gfffile)
    anno.filterSize(minSize, maxSize)
    pc = parclipsites.ParclipSites('')
    pc.loadFromFile(parclipfile)
    fc_out = open(outputfile, 'w')
    for g in range(anno.size()):
        if verbose:
            functions.showProgress(g, (anno.size()-1), vstring)
        if anno.strand[g] == '+':
            values_upstream = pc.getValues(anno.chr[g], anno.start[g], anno.strand[g], sense, upstream, gene)
            values_dostream = pc.getValues(anno.chr[g], anno.stop[g], anno.strand[g], sense, gene, downstream)
        else:
            values_upstream = pc.getValues(anno.chr[g], anno.stop[g], anno.strand[g], sense, upstream, gene)
            values_dostream = pc.getValues(anno.chr[g], anno.start[g], anno.strand[g], sense, gene, downstream)
        if values_upstream != None and values_dostream != None:
            for v in values_upstream:
                fc_out.write(str(v)+'\t')
            for v in values_dostream:
                fc_out.write(str(v)+'\t')
            fc_out.write('\n')
    fc_out.close()
    if verbose:
        print('')

def run():
    scriptPath = os.path.dirname(os.path.realpath(__file__))
    scriptPath = scriptPath+'/'
    parser = argparse.ArgumentParser(description='Plot PAR-CLIP data in sense and anti-sense direction around start and stop positions given in the GFF file. Outputfilename is constructed from parameters: outputdir/prefix+u+g+d+min+max+labelCenterA+labelBody+labelCenterB+plotSmooth.pdf', epilog="contact: torkler@genzentrum.lmu.de")
    parser.add_argument('parclip', help='path to the PAR-CLIP *.table')
    parser.add_argument('outputdir', help='output directory')
    parser.add_argument('prefix', help='prefix of filenames')
    parser.add_argument('gff', help='GFF file used for plotting')
    parser.add_argument('-d', help='set downstream range [default: 1000nt]', dest='downstream', default=1000, type=int)
    parser.add_argument('-u', help='set upstream range [default: 1000nt]', dest='upstream', default=1000, type=int)
    parser.add_argument('-g', help='set gene range [default: 750nt]', dest='gene', default=750, type=int)
    parser.add_argument('--min', help='minium transcript size [default: 0nt]',  default=0, type=int)
    parser.add_argument('--max', help='maximum transcript size [default: 5000nt]', default=5000, type=int)
    parser.add_argument('--plotSmooth', help='half of the window size used for the running mean [default: 20nt]', default=20, type=int)
    parser.add_argument('--labelCenterA', help='plot label for the first center position [default: TSS]', default='TSS')
    parser.add_argument('--labelBody', help='plot label for body (between A and B) [default: gene]', default='gene')
    parser.add_argument('--labelCenterB', help='plot label for the second center position [default: pA]', default='pA')
    parser.add_argument('-r','--remove', dest='remove', action="store_true", default=False, help='remove temporary text files. [default: false]')
    parser.add_argument('-v','--verbose', dest='verbose', action="store_true", default=False, help='verbose output')
    args = parser.parse_args()
    outfile_sense = args.outputdir+args.prefix+'_centerBoth_up'+str(args.upstream)+'_gene'+str(args.gene)+'_do'+str(args.downstream)+'_min'+str(args.min)+'_max'+str(args.max)+'_'+str(args.labelCenterA)+'_'+str(args.labelBody)+'_'+str(args.labelCenterB)+'_sense.table'
    outfile_asense = args.outputdir+args.prefix+'_centerBoth_up'+str(args.upstream)+'_gene'+str(args.gene)+'_do'+str(args.downstream)+'_min'+str(args.min)+'_max'+str(args.max)+'_'+str(args.labelCenterA)+'_'+str(args.labelBody)+'_'+str(args.labelCenterB)+'_asense.table'
    outfile_pdf = args.outputdir+args.prefix+'_centerBoth_up'+str(args.upstream)+'_gene'+str(args.gene)+'_do'+str(args.downstream)+'_min'+str(args.min)+'_max'+str(args.max)+'_'+str(args.labelCenterA)+'_'+str(args.labelBody)+'_'+str(args.labelCenterB)+'_sm'+str(args.plotSmooth)+'.pdf'
    
    if args.verbose:
        a = datetime.datetime.now()
    main(args.parclip, outfile_sense, args.gff, args.downstream, args.upstream, args.gene, True, args.min, args.max, args.verbose, 'Collecting PAR-CLIP sense data')
    main(args.parclip, outfile_asense, args.gff, args.downstream, args.upstream, args.gene, False, args.min, args.max, args.verbose, 'Collecting PAR-CLIP anti-sense data')
    if args.verbose:
        b = datetime.datetime.now()
        c = b - a
        print('')
        print('time: '+str(c.seconds)+' seconds')
    
    os.system('R -q --slave -f '+scriptPath+'plotCenterBothEnds.R --args '+outfile_sense+' '+outfile_asense+' '+outfile_pdf+' '+args.prefix+' '+str(args.upstream)+' '+str(args.downstream)+' '+str(args.gene)+' '+str(args.plotSmooth)+' '+str(args.labelCenterA)+' '+str(args.labelBody)+' '+str(args.labelCenterB))
    if args.remove:
        os.remove(outfile_sense)
        os.remove(outfile_asense)

if __name__ == '__main__':
    run()





