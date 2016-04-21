import argparse
from collections import defaultdict
import re
import os
import itertools
import functools

import pysam
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('calmd_bam_file')
    parser.add_argument('output_dir')
    return parser


def create_transition_plots(calmd_bam_file, output_dir):
    data_dic = defaultdict(list)
    lengths = []
    nuc_compl = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
    }
    mm_pat = re.compile('([0-9]+)([ACGTacgt])([0-9]+)')
    with pysam.AlignmentFile(calmd_bam_file, 'rb') as bam_file:
        for entry in bam_file.fetch():
            tag_dict = dict(entry.tags)
            if entry.is_unmapped:
                continue
            lengths.append(entry.inferred_length)

            if 'NM' in tag_dict:
                nm = tag_dict['NM']
            elif 'nM' in tag_dict:
                nm = tag_dict['nM']
            else:
                raise ValueError('NM field not found')

            if nm != 1:
                continue

            if 'MD' not in tag_dict:
                raise ValueError('MD field not found')
            hit = mm_pat.match(tag_dict['MD'])
            mm_pos = int(hit.group(1)) + 1
            ref_nuc = hit.group(2).upper()
            mm_nuc = entry.seq[mm_pos - 1].upper()
            if entry.is_reverse:
                transition = nuc_compl[ref_nuc] + nuc_compl[mm_nuc]
            else:
                transition = ref_nuc + mm_nuc

            data_dic['length'].append(entry.inferred_length)
            data_dic['transition'].append(transition)
            data_dic['mm_pos'].append(mm_pos)

    df = pd.DataFrame(data_dic)

    max_len = np.max(lengths)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(lengths, bins=np.arange(max_len) + 1)
    fig.suptitle('Read size distribution of mapped reads', fontsize=20)
    ax.set_xlabel('Read length in bp', fontsize=16)

    maplen_plot = os.path.join(output_dir, 'mapped_lengths.png')
    plt.savefig(maplen_plot)

    colors = [
        '#a6cee3',
        '#1f78b4',
        '#b2df8a',
        '#33a02c',
        '#fb9a99',
        '#e31a1c',
        '#fdbf6f',
        '#ff7f00',
        '#cab2d6',
        '#6a3d9a',
        '#b15928',
        '#ffff99',
    ]

    def count_nuc(ser, nuc):
        return len(ser[ser == nuc])

    dinucs = [''.join(d) for d in itertools.permutations('ACGT', 2)]
    tr_map = {}
    for dinuc in dinucs:
        tr_map[dinuc] = functools.partial(count_nuc, nuc=dinuc)

    tr_plot_dir = os.path.join(output_dir, 'transition_plots')
    if not os.path.exists(tr_plot_dir):
        os.makedirs(tr_plot_dir)

    for length, mm_df in df.groupby('length'):
        # too noisy for a useful plot
        if len(mm_df) < length * 10:
            continue
        dinuc_tr_df = mm_df.groupby('mm_pos').agg({
            'transition': tr_map
        })['transition']

        fig, ax = plt.subplots(figsize=(10, 8))
        for i, dinuc in enumerate(dinucs):
            ax.plot(dinuc_tr_df[dinuc], color=colors[i])
        lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, fontsize=16)
        title = fig.suptitle('Observed transitions in %s bp reads' % length, fontsize=20)
        ax.set_xlabel('Read Position', fontsize=16)
        ax.set_xlim(0, length + 1)

        tr_plot = os.path.join(tr_plot_dir, 'transition_%sbp_plot.png' % length)
        plt.savefig(tr_plot, bbox_extra_artists=(lgd, title), bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    create_transition_plots(args.calmd_bam_file, args.output_dir)
