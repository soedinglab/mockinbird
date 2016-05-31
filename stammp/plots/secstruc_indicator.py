import argparse
import os
import random
import re
import shutil

import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from stammp.obj import genome, parclipsites, gff, functions
from stammp.utils import argparse_helper as aph
from stammp.utils import execute


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('pc_sites', type=aph.file_r)
    parser.add_argument('genome_fasta', type=aph.file_r)
    parser.add_argument('sample_gff', type=aph.file_r)
    parser.add_argument('--prefix')
    parser.add_argument('output_dir', type=aph.dir_rwx_create)
    parser.add_argument('--start', help='start index of PAR-CLIP sites',
                        type=int, default=0)
    parser.add_argument('--stop', help='stop index of PAR-CLIP sites',
                        type=int, default=1500)
    parser.add_argument('--width', help='number of nt +/- the crosslink site',
                        type=int, default=15)
    sort_key_help = ('set key that is used for PAR-CLIP site ordering [default = occ], '
                     'options: [occ, m, r, mr, pvalue]')
    sort_keys = ['occ', 'm', 'r', 'mr', 'pvalue']
    parser.add_argument('--key', help=sort_key_help, choices=sort_keys, default='occ')
    parser.add_argument('--keep-tmp-files', action='store_true')
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    if not args.prefix:
        bname = os.path.basename(args.pc_sites)
        prefix, ext = os.path.splitext(bname)
    else:
        prefix = args.prefix

    tmp_dir = os.path.join(args.output_dir, 'tmp_ss_energy')
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)

    full_genome = genome.Genome(args.genome_fasta, verbose=False)
    sites = parclipsites.ParclipSites()
    sites.loadFromFile(args.pc_sites)

    sites.sort(args.key)
    xsite_fasta = os.path.join(tmp_dir, prefix + '.fa')
    sites.save2Fasta(full_genome, xsite_fasta, args.start, args.stop, width=args.width)

    n_sample = args.stop - args.start
    sample_anno = gff.GFF(args.sample_gff)
    rnd_seqs = generate_rnd_seqs(sample_anno, full_genome, n_sample, args.width)
    rnd_fasta = os.path.join(tmp_dir, 'random_seqs.fa')

    with open(rnd_fasta, 'w') as rnd_file:
        for i, seq in enumerate(rnd_seqs):
            print('>rnd_seq_%s' % i, file=rnd_file)
            print(seq, file=rnd_file)

    xsite_pred = os.path.join(tmp_dir, prefix + '.pred')
    xsite_cmd = [
        'RNAfold',
        '-i %r' % xsite_fasta,
        '--noPS'
        '> %r' % xsite_pred
    ]
    execute(xsite_cmd)

    rnd_pred = os.path.join(tmp_dir, 'random_seqs.pred')
    rnd_cmd = [
        'RNAfold',
        '-i %r' % rnd_fasta,
        '--noPS'
        '> %r' % rnd_pred
    ]
    execute(rnd_cmd)

    xsite_energy = read_energy(xsite_pred)
    rnd_energy = read_energy(rnd_pred)

    energies = np.hstack([xsite_energy, rnd_energy])
    label = ['Crosslink'] * len(xsite_energy) + ['Random'] * len(rnd_energy)

    df = pd.DataFrame({
        'Sites': pd.Series(label, dtype='category'),
        'Energy': energies,
    })

    png_file = os.path.join(args.output_dir, '%s_ss_energy.png' % prefix)
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(data=df, y='Sites', x='Energy', ax=ax)
    fig.savefig(png_file)
    plt.close()

    if not args.keep_tmp_files:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def generate_rnd_seqs(annot, genome, n_sample, width):
    seqs = []
    cur_n = 0
    while cur_n < n_sample:
        rnd_anno = random.randint(0, (annot.size() - 1))
        rnd_pos = random.randint(annot.start[rnd_anno], annot.stop[rnd_anno])
        tmp_seq = genome.getSequence(annot.chr[rnd_anno], rnd_pos - width, rnd_pos + width)
        if tmp_seq != -1:
            if annot.strand[rnd_anno] == '+':
                seqs.append(tmp_seq)
            else:
                seqs.append(functions.makeReverseComplement(tmp_seq))
            cur_n += 1
    return seqs


def read_energy(predfile):
    regex = re.compile('\((\s*-?[0-9]+\.[0-9]+\s*)\)')
    energies = []
    with open(predfile) as f:
        for line in f:
            if line.startswith('>'):
                f.readline()
                pred_line = f.readline()
                hit = regex.search(pred_line)
                assert hit is not None
                energy = float(hit.group(1))
                energies.append(energy)
    return -np.array(energies)



if __name__ == '__main__':
    main()
