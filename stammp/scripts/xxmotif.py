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


.. image:: img/img_pub1_xxmotif_start0_stop1000_width12_sort_occ_0.png
   :align: center
   :width: 700px
   :alt: alternate text

"""
import argparse
import os
from stammp.obj import parclipsites, genome, gff
from stammp.utils import argparse_helper as aph
from stammp.utils import execute


def run():
    scriptPath = os.path.dirname(os.path.realpath(__file__))
    parser = argparse.ArgumentParser(
        description=('Selects sequences from PAR-CLIP sites and pass them for '
                     'motif search to XXmotif.')
    )
    parser.add_argument('inputfile', help='PAR-CLIP file *.table')
    parser.add_argument('genome', help='path to genome', type=aph.file_r)
    parser.add_argument('outdir', help='output directory', type=aph.dir_rwx)
    parser.add_argument('prefix', help='prefix')
    parser.add_argument('--negSet', help='set path to negative set if available.')
    plot_pwm_help = 'plot top plotPWM PWMs as pdf sequence logos. [default=0]',
    parser.add_argument('--plotPWM', help=plot_pwm_help, type=int, default=0)
    parser.add_argument('--start', help='start index of PAR-CLIP sites [default=0]',
                        type=int, default=0)
    parser.add_argument('--stop', help='stop index of PAR-CLIP sites [default=1500]',
                        type=int, default=1500)
    parser.add_argument('--width', help='number of nt +/- the crosslink site [default=12]',
                        type=int, default=12)
    sort_key_help = ('set key that is used for PAR-CLIP site ordering [default = occ], '
                     'options: [occ, m, r, mr, pvalue]')
    sort_keys = ['occ', 'm', 'r', 'mr', 'pvalue']
    parser.add_argument('--key', help=sort_key_help, choices=sort_keys, default='occ')
    filter_gff_help = ('set path to GFF if sites should be removed that overlap '
                       'with the GFF. Default = \'\' means that no sites are filtered out.')
    parser.add_argument('--filterGFF', help=filter_gff_help, default='')
    awidth_help = 'number of nt that are added to the start/stop indices of the GFF annotations'
    parser.add_argument('--awidth', help=awidth_help, type=int, default=20)
    args = parser.parse_args()

    prefix_pat = '%s_xxmotif_start%s_stop%s_width%s_sort_%s'
    file_prefix = prefix_pat % (args.prefix, args.start, args.stop, args.width, args.key)

    gen = genome.Genome(args.genome, verbose=False)
    sites = parclipsites.ParclipSites()
    sites.loadFromFile(args.inputfile)

    if args.filterGFF != '':
        anno = gff.GFF(args.filterGFF)
        sites = sites.removeSitesLocatedInGFF(anno, args.awidth)

    sites.sort(args.key)
    gen_file = os.path.join(args.outdir, file_prefix + '.fa')
    sites.save2Fasta(gen, gen_file, args.start, args.stop, width=args.width)

    cmd = [
        'XXmotif',
        args.outdir,
        gen_file,
        '--zoops',
        '--merge-motif-threshold LOW',
        '--max-match-positions 10',
    ]
    if args.negSet:
        cmd.append('--negSet %s' % args.negSet)
    execute(cmd)

    plot_script = os.path.join(scriptPath, '..', 'plots', 'weblogo.R')
    pwm_file = os.path.join(args.outdir, file_prefix + '.pwm')
    plot_cmd = [
        'R',
        '-q',
        '--slave',
        '-f %s' % plot_script,
        '--args',
        pwm_file,
        args.outdir,
        file_prefix,
        args.plotPWM,
    ]
    if args.plotPWM > 0:
        execute(plot_cmd)


if __name__ == '__main__':
    run()
