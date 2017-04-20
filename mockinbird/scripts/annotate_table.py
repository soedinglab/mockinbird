import argparse
import os
from collections import Counter, defaultdict
from itertools import chain

import numpy as np

from pyivtree import IVTree, GenomicInterval as Interval

from stammp.utils.argparse_helper import file_r, dir_rwx_create
from stammp.utils.parsers import GFF3Parser, PCTableParser


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('parclip_table', type=file_r, help='path to parclip table')
    parser.add_argument('output_dir', type=dir_rwx_create, help='output directory')
    parser.add_argument('--max_n', type=int, default=np.inf)
    parser.add_argument('gff3_annot', nargs='+', type=file_r,
                        help='paths to one or more gff3 annotation files')
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    gff_tree = defaultdict(IVTree)
    for gff_file in args.gff3_annot:
        with open(gff_file) as gff:
            parser = GFF3Parser(gff)
            for record in parser.parse():
                iv = Interval(record.start, record.end)
                iv.rec = record
                gff_tree[record.seqid].insert(iv)

    table = []
    occ_vector = []
    with open(args.parclip_table) as in_table:
        pc_parser = PCTableParser(in_table)
        for rec in pc_parser.parse():
            table.append(rec)
            occ_vector.append(-rec.occupancy)
    sort_vec = np.argsort(occ_vector)
    ranks = np.empty(len(occ_vector), int)
    ranks[sort_vec] = np.arange(len(occ_vector))

    table_name = os.path.basename(args.parclip_table)
    annot_table = os.path.join(args.output_dir, table_name + '_annot')
    with open(annot_table, 'w') as out_table:
        print(*chain(pc_parser._fields, ['annotation']), sep='\t', file=out_table)
        hit_counter = Counter()
        ambig_hits = 0
        unannotated = 0
        for ind, rec in zip(ranks, table):
            if ind >= args.max_n:
                continue
            iv = Interval(rec.position, rec.position)
            hits = set()
            for gff_rec in gff_tree[rec.seqid].query_all_overlaps(iv):
                if rec.strand == gff_rec.rec.strand:
                    hits.add(gff_rec.rec.type)

            if len(hits) > 1:
                ambig_hits += 1
            elif len(hits) == 0:
                unannotated += 1
            else:
                hit, = hits
                hit_counter[hit] += 1

            if len(hits) == 0:
                hits.add('NA')
            print(*chain(list(rec), ['|'.join(hits)]), sep='\t', file=out_table)

    summary_file = os.path.join(args.output_dir, 'summary.tab')
    with open(summary_file, 'w') as sum_file:
        for annot, count in hit_counter.items():
            print(annot, count, sep='\t', file=sum_file)
        print('unannotated', unannotated, sep='\t', file=sum_file)
        print('ambiguous', ambig_hits, sep='\t', file=sum_file)


if __name__ == '__main__':
    main()
