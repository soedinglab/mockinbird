import argparse
import functools
import itertools
from collections import defaultdict
import os
import io

import pysam
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


CIGAR_MATCH = 0
CIGAR_SOFTCLIP = 4

DINUCS = [''.join(d) for d in itertools.permutations('ACGT', 2)]
DINUC_COLORS = {
    'AC': '#a6cee3',
    'AG': '#1f78b4',
    'AT': '#b2df8a',
    'CA': '#33a02c',
    'CG': '#fb9a99',
    'CT': '#e31a1c',
    'GA': '#fdbf6f',
    'GC': '#ff7f00',
    'GT': '#cab2d6',
    'TA': '#6a3d9a',
    'TC': '#b15928',
    'TG': '#ffff99',
}
NUC_COMPL = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
}


def parse_md(md_str):
    char_pos = 0
    read_pos = 0
    for i, c in enumerate(md_str):
        if c.isalpha():
            read_pos += int(md_str[char_pos:i])
            yield c.upper(), read_pos
            read_pos += 1
            char_pos = i + 1


def get_map_pos(read):
    start = 0
    end = read.infer_query_length()
    assert read.infer_query_length() == len(read.seq)
    op, bp = read.cigartuples[0]
    if op == CIGAR_SOFTCLIP:
        start = bp
    op, bp = read.cigartuples[-1]
    if op == CIGAR_SOFTCLIP:
        end -= bp
    return start, end, end - start


def create_parser():
    parser = argparse.ArgumentParser('bam_postprocess')
    parser.add_argument('input_bam_file')
    parser.add_argument('output_bam_file')
    parser.add_argument('output_directory')
    parser.add_argument('--plot_transition_profiles', action='store_true')
    parser.add_argument('--min-length', default=0, type=int)
    parser.add_argument('--mut_edge_bp', default=0, type=int)
    parser.add_argument('--max_transitions', default=1, type=int)
    return parser


def run(args):

    pre_modules = []
    if args.plot_transition_profiles:
        pre_tr_dir = os.path.join(args.output_directory, 'pre_fil_data')
        if not os.path.exists(pre_tr_dir):
            os.makedirs(pre_tr_dir)
        pre_modules.append(MutationModule(
            pre_tr_dir)
        )

    post_modules = []
    if args.plot_transition_profiles:
        post_tr_dir = os.path.join(args.output_directory, 'post_fil_data')
        if not os.path.exists(post_tr_dir):
            os.makedirs(post_tr_dir)
        post_modules.append(MutationModule(
            post_tr_dir)
        )
    filters = []
    filters.append(SizeFilter(
        args.min_length)
    )
    filters.append(EdgeClippingFilter(
        args.mut_edge_bp,
        args.max_transitions)
    )

    with pysam.AlignmentFile(args.input_bam_file, 'rb') as infile:
        with pysam.AlignmentFile(args.output_bam_file, 'wb', template=infile) as outfile:
            for read in infile:
                # prefilter modules
                for mod in pre_modules:
                    mod.update(read)

                # filtering
                cur_read = read
                for fil in filters:
                    cur_read = fil.filter(cur_read)
                    if cur_read is None:
                        break
                if cur_read is not None:
                    outfile.write(cur_read)

                    # postfilter modules
                    for mod in post_modules:
                        mod.update(cur_read)

    for mod in pre_modules:
        mod.aggregate()

    for mod in post_modules:
        mod.aggregate()


class SizeFilter:
    def __init__(self, min_size):
        self._min_size = min_size

    def filter(self, read):
        if read.is_unmapped:
            return None
        start, end, real_length = get_map_pos(read)
        if real_length < self._min_size:
            return None
        return read


class EdgeClippingFilter:
    def __init__(self, edge_bp, allowed_mm):
        self._edge_bp = edge_bp
        self._allowed_mm = allowed_mm

    def filter(self, read):
        tag_dict = dict(read.tags)
        assert 'NM' in tag_dict
        nm = tag_dict['NM']
        assert 'MD' in tag_dict
        md = tag_dict['MD']

        if nm == 0:
            return read

        start, end, real_length = get_map_pos(read)
        mm = []
        left_clip_bp = 0
        right_clip_bp = 0
        for ref_nuc, mm_pos in parse_md(md):
            if mm_pos < self._edge_bp:
                left_clip_bp = max(left_clip_bp, mm_pos + 1)
                nm -= 1
            elif mm_pos >= real_length - self._edge_bp:
                right_clip_bp = max(right_clip_bp, real_length - mm_pos)
                nm -= 1
            else:
                mm.append((ref_nuc, mm_pos))

        if len(mm) > self._allowed_mm:
            return None

        # no additional clipping needed
        if left_clip_bp == 0 and right_clip_bp == 0:
            return read

        # update md
        md_io = io.StringIO()
        cur_bps = 0
        for ref_nuc, mm_pos in mm:
            add_bp = mm_pos - left_clip_bp - cur_bps
            cur_bps += add_bp + 1
            md_io.write(str(add_bp))
            md_io.write(ref_nuc)
        total_bp = real_length - left_clip_bp - right_clip_bp
        md_io.write(str(total_bp - cur_bps))
        md_io.seek(0)
        md_str = md_io.read()

        tags = list(read.tags)
        for i, (key, value) in enumerate(tags):
            if key == 'NM':
                nm_ind = i
            elif key == 'MD':
                md_ind = i
        tags[nm_ind] = ('NM', nm)
        tags[md_ind] = ('MD', md_str)
        read.tags = tags

        # update cigar
        cigar_tup = read.cigartuples
        if left_clip_bp > 0:
            first_op, first_bp = cigar_tup[0]
            if first_op == CIGAR_SOFTCLIP:
                second_op, second_bp = cigar_tup[1]
                assert second_op == CIGAR_MATCH
                if second_bp > left_clip_bp:
                    cigar_tup[0] = (CIGAR_SOFTCLIP, first_bp + left_clip_bp)
                    cigar_tup[1] = (CIGAR_MATCH, second_bp - left_clip_bp)
                else:
                    return None
            else:
                assert first_op == CIGAR_MATCH
                if first_bp > left_clip_bp:
                    cigar_tup[0] = (CIGAR_MATCH, first_bp - left_clip_bp)
                    cigar_tup.insert(0, (CIGAR_SOFTCLIP, left_clip_bp))
                else:
                    return None
        if right_clip_bp > 0:
            last_op, last_bp = cigar_tup[-1]
            if last_op == CIGAR_SOFTCLIP:
                penu_op, penu_bp = cigar_tup[-2]
                assert penu_op == CIGAR_MATCH
                if penu_bp > right_clip_bp:
                    cigar_tup[-1] = (CIGAR_SOFTCLIP, last_bp + right_clip_bp)
                    cigar_tup[-2] = (CIGAR_MATCH, penu_bp - right_clip_bp)
                else:
                    return None
            else:
                assert last_op == CIGAR_MATCH
                if last_bp > right_clip_bp:
                    cigar_tup[-1] = (CIGAR_MATCH, last_bp - right_clip_bp)
                    cigar_tup.append((CIGAR_SOFTCLIP, right_clip_bp))
                else:
                    return None

        read.cigartuples = cigar_tup
        read.reference_start += left_clip_bp
        return read


class MutationModule:

    def __init__(self, output_dir):
        self._output_dir = output_dir
        self._real_lengths = []
        self._tr_data = defaultdict(list)

    def update(self, read):

        if read.is_unmapped:
            return

        tag_dict = dict(read.tags)
        assert 'NM' in tag_dict
        nm = tag_dict['NM']
        assert 'MD' in tag_dict
        md = tag_dict['MD']

        start, end, real_length = get_map_pos(read)
        self._real_lengths.append(real_length)

        if nm == 0:
            return

        # extract transition
        for ref_nuc, mm_pos in parse_md(md):
            mm_nuc = read.seq[start + mm_pos].upper()

            if read.is_reverse:
                transition = NUC_COMPL[ref_nuc] + NUC_COMPL[mm_nuc]
            else:
                transition = ref_nuc + mm_nuc

            self._tr_data['length'].append(real_length)
            self._tr_data['mm_pos'].append(mm_pos + 1)
            self._tr_data['transition'].append(transition)

    def aggregate(self):
        # plot read size distribution
        max_len = max(self._real_lengths)
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(self._real_lengths, bins=np.arange(max_len) + 1)
        fig.suptitle('Mapped read size distribution', fontsize=20)
        ax.set_xlabel('Mapped read length in bp', fontsize=16)
        maplen_plot = os.path.join(self._output_dir, 'mapped_lengths.png')
        plt.savefig(maplen_plot)
        plt.close()

        mm_plot_dir = os.path.join(self._output_dir, 'mismatch_profiles')
        if not os.path.exists(mm_plot_dir):
            os.makedirs(mm_plot_dir)

        tr_plot_dir = os.path.join(self._output_dir, 'transition_profiles')
        if not os.path.exists(tr_plot_dir):
            os.makedirs(tr_plot_dir)

        def count_nuc(ser, nuc):
            return len(ser[ser == nuc])

        tr_map = {}
        for dinuc in DINUCS:
            tr_map[dinuc] = functools.partial(count_nuc, nuc=dinuc)

        # position mismatch histograms
        mm_df = pd.DataFrame(self._tr_data)
        for length, df in mm_df.groupby('length'):

            # the profiles will be too noisy
            if len(df) < 20 * length:
                continue

            # mutation profile
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.hist(df.mm_pos, bins=np.arange(length) + 1)
            fig.suptitle('Mismatch frequency distribution of %s bp reads' % length,
                         fontsize=20)
            ax.set_xlabel('Read position', fontsize=16)

            mismatch_freq_plot = os.path.join(mm_plot_dir,
                                              'mismatch_freq_%sbp.png' % length)
            plt.savefig(mismatch_freq_plot)
            plt.close()

            # transition profile
            dinuc_tr_df = df.groupby('mm_pos').agg({
                'transition': tr_map
            })['transition']

            fig, ax = plt.subplots(figsize=(10, 8))
            for dinuc in DINUCS:
                ax.plot(dinuc_tr_df[dinuc], color=DINUC_COLORS[dinuc])
            lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, fontsize=16)
            title = fig.suptitle('Observed transitions in %s bp reads' % length,
                                 fontsize=20)
            ax.set_xlabel('Read Position', fontsize=16)
            ax.set_xlim(0, length + 1)

            tr_plot = os.path.join(tr_plot_dir, 'transition_%sbp_plot.png' % length)
            plt.savefig(tr_plot, bbox_extra_artists=(lgd, title), bbox_inches='tight')
            plt.close()


def main():
    parser = create_parser()
    args = parser.parse_args()
    run(args)

if __name__ == '__main__':
    main()
