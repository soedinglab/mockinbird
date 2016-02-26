#! /usr/bin/python3
"""
Selects sequences from PAR-CLIP sites and pass them for motif search to
XXmotif.

.. note:: It is recommended to use negative sets for motif analysis of PAR-CLIP data sets as provided by :mod:`~stammp.scripts.makeNegSets`.

**Usage:** stammp-xxmotif [-h] [--negSet NEGSET] [--plotPWM PLOTPWM]
                      [--start START] [--stop STOP] [--width WIDTH]
                      [--key KEY] [--filterGFF FILTERGFF] [--awidth AWIDTH]
                      inputfile genome outdir prefix

**Positional arguments:**
  =========             =====================
  inputfile             PAR-CLIP file \*.table
  genome                path to genome
  outdir                output directory
  prefix                prefix
  =========             =====================

**Optional arguments:**
  ===================== ======================================================
  -h, --help            show this help message and exit
  --negSet NEGSET       set path to negative set if available. [default =
                        None]
  --plotPWM PLOTPWM     plot top plotPWM PWMs as pdf sequence logos. [default
                        = 0]
  --start START         start index of PAR-CLIP sites [default=0]
  --stop STOP           stop index of PAR-CLIP sites [default=1500]
  --width WIDTH         number of nt +/- the crosslink site [default=50]
  --key KEY             set key that is used for PAR-CLIP site ordering
                        [default = occ], options: [occ, m, r, mr, pvalue]
  --filterGFF FILTERGFF
                        set path to GFF if sites should be removed that
                        overlap with the GFF. Default = '' means that no sites
                        are filtered out.
  --awidth AWIDTH       number of nt that are added to the start/stop indices
                        of the GFF annotations
  ===================== ======================================================

Example::
    
    $ stammp-xxmotif parclip.table genome.fa outdir/ prefix --start 0 --stop 1000 --plotPWM 3
    

.. image:: tut_pub1_xxmotif_start0_stop1000_width12_sort_occ_0.png
   :align: center
   :width: 700px
   :alt: alternate text

"""
import argparse
import os
import sys
from stammp.obj import *

def run():
    scriptPath = os.path.dirname(os.path.realpath(__file__))
    scriptPath = scriptPath+'/'
    parser = argparse.ArgumentParser(description='Selects sequences from PAR-CLIP sites and pass them for motif search to XXmotif.', epilog="contact: torkler@genzentrum.lmu.de")
    parser.add_argument('inputfile',help='PAR-CLIP file *.table')
    parser.add_argument('genome',   help='path to genome')
    parser.add_argument('outdir',   help='output directory')
    parser.add_argument('prefix',   help='prefix')
    parser.add_argument('--negSet', help='set path to negative set if available. [default = None]', default=None)
    parser.add_argument('--plotPWM',help='plot top plotPWM PWMs as pdf sequence logos. [default = 0]', type=int, default = 0)
    parser.add_argument('--start',  help='start index of PAR-CLIP sites [default=0]', type=int, default = 0)
    parser.add_argument('--stop',   help='stop index of PAR-CLIP sites [default=1500]', type=int, default = 1500)
    parser.add_argument('--width',  help='number of nt +/- the crosslink site [default=12]', type=int, default = 12)
    parser.add_argument('--key',    help='set key that is used for PAR-CLIP site ordering [default = occ], options: [occ, m, r, mr, pvalue]', default='occ')
    parser.add_argument('--filterGFF',  help='set path to GFF if sites should be removed that overlap with the GFF. Default = \'\' means that no sites are filtered out.', default='')
    parser.add_argument('--awidth', help='number of nt that are added to the start/stop indices of the GFF annotations', type=int, default = 20)
    args = parser.parse_args()
    
    fileprefix = args.outdir+args.prefix+'_xxmotif_start'+str(args.start)+'_stop'+str(args.stop)+'_width'+str(args.width)+'_sort_'+args.key
    
    if os.path.exists(args.outdir) == False:
        os.makedirs(args.outdir)
    
    yeast   = genome.Genome(args.genome, False)
    sites   = parclipsites.ParclipSites('')
    sites.loadFromFile(args.inputfile)
    
    if args.filterGFF != '':
        anno     = gff.GFF(args.filterGFF)
        sites   = sites.removeSitesLocatedInGFF(anno, args.awidth)
    
    sites.sort(args.key)
    sites.save2Fasta(yeast, fileprefix+'.fa', args.start, args.stop, width=args.width)
    cmd = 'XXmotif '+args.outdir+' '+fileprefix+'.fa'+' --zoops --merge-motif-threshold LOW --max-match-positions 10'
    if args.negSet != None:
        cmd = cmd+' --negSet '+args.negSet
    os.system(cmd)
    
    if args.plotPWM > 0:
        os.system('R -q --slave -f '+scriptPath[0:(len(scriptPath)-8)]+'plots/'+'weblogo.R --args '+os.path.abspath(fileprefix+'.pwm')+' '+os.path.abspath(args.outdir)+'/ '+args.prefix+' '+str(args.plotPWM))

if __name__ == '__main__':
    run()
