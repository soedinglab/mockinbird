"""
Select all PAR-CLIP sites that are located within the given GFF

**Usage:** mockinbird-selectSitesInside [-h] [--min MIN] [--max MAX]
                                [--upstream UPSTREAM]
                                [--downstream DOWNSTREAM]
                                parclip outputfile gff
**Positional arguments:**
  ==========            ======================
  parclip               PAR-CLIP file \*.table
  outputfile            PAR-CLIP file \*.table
  gff                   GFF file
  ==========            ======================

**Optional arguments:**
  ======================= =========================================
  -h, --help              show this help message and exit
  --min INT               minium transcript size [default: 0nt]
  --max INT               maximum transcript size [default: 5000nt]
  --upstream INT          additional nt upstream [default: 0nt]
  --downstream INT        additional nt downstream [default: 0nt]
  -v, --verbose           verbose output
  ======================= =========================================
"""
import argparse
from mockinbird.obj import *
from mockinbird.utils import ParclipSiteContainer

def main(inputfile, outputfile, gfffile, gffmin, gffmax, upstream,
         downstream, verbose):
    sites = ParclipSiteContainer()
    sites.loadFromFile(inputfile)
    anno = gff.GFF(gfffile)
    anno.filterSize(gffmin, gffmax)
    anno.getChromosomePositions()
    if anno.size() < 10:
        print('Warning: Low number of annotation enries! '+str(anno.size()))
    fsites = ParclipSiteContainer()
    percent_old = 0
    percent_new = 0
    for i in range(sites.size()):
        if anno.isInside(sites.chrs[i], sites.pos[i], sites.strand[i], 
                         upstream, downstream)[1]:
            fsites.addSite(sites.chrs[i], sites.pos[i], sites.m[i], sites.r[i], 
                           sites.result[i], sites.strand[i], sites.occ[i])
        percent_new = round(i/sites.size()*100)
        if percent_new > percent_old:
            if verbose: functions.showProgress(i, anno.size(), 'selecting sites')
            percent_old = percent_new
    fsites.save2File(outputfile)

def run():
    parser = argparse.ArgumentParser(description='Select all PAR-CLIP sites'
    + ' that are located within the given GFF', epilog="contact: torkler@genzentrum.lmu.de")
    parser.add_argument('parclip', help='PAR-CLIP file *.table')
    parser.add_argument('outputfile', help='PAR-CLIP file *.table')
    parser.add_argument('gff', help='GFF file')
    parser.add_argument('--min', help='minium transcript size [default: 0nt]',
                        default=0, type=int)
    parser.add_argument('--max', help='maximum transcript size [default: 5000nt]',
                        default=5000, type=int)
    parser.add_argument('--upstream', help='additional nt upstream [default: 0nt]',
                        default=0, type=int)
    parser.add_argument('--downstream', help='additional nt downstream [default: 0nt]',
                        default=0, type=int)
    parser.add_argument('-v','--verbose', dest='verbose', 
                        action="store_true", default=False, help='verbose output')
    args = parser.parse_args()
    functions.checkExistence(args.parclip)
    functions.checkExistence(args.gff)
    main(args.parclip, args.outputfile, args.gff, args.min, args.max, 
         args.upstream, args.downstream, args.verbose)

if __name__ == '__main__':
    run()
