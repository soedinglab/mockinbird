import argparse
from collections import defaultdict
from os import path
import random

import stammp.utils.argparse_helper as aph
from stammp.utils.parsers import GFF3Parser
from stammp.utils.helper_objects import EfficientGenome, ParclipSiteContainer
from pyivtree import GenomicInterval as Interval
from pyivtree import IVTree


def create_parser():
    parser = argparse.ArgumentParser('table2graphprot')
    parser.add_argument('table_file', type=aph.file_r)
    parser.add_argument('gff_file', type=aph.file_r)
    parser.add_argument('genome_fasta', type=aph.file_r)
    parser.add_argument('--viewpoint_bp', type=int, default=11)
    parser.add_argument('--flanking_bp', type=int, default=150)
    parser.add_argument('output_directory', type=aph.dir_rwx)
    parser.add_argument('--top_n_sites', type=int, default=2000)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    table_file_name = path.basename(args.table_file)
    base_name, ext = path.splitext(table_file_name)
    neg_fpath = path.join(args.output_directory, base_name + '_neg.fa')
    pos_fpath = path.join(args.output_directory, base_name + '_pos.fa')

    mRNA_trees = defaultdict(IVTree)

    with open(args.gff_file) as gff:
        parser = GFF3Parser(gff)
        for rec in parser.parse():
            iv = Interval(rec.start, rec.end)
            iv.strand = rec.strand
            mRNA_trees[rec.seqid].insert(iv)

    vp_bp = args.viewpoint_bp
    fl_bp = args.flanking_bp
    pc_sites = ParclipSiteContainer.from_file(args.table_file)
    pc_sites.sort(by='occupancy', ascending=False)
    pc_sites = pc_sites[:args.top_n_sites]

    seq_counter = 0
    with EfficientGenome(args.genome_fasta) as genome, \
            open(pos_fpath, 'w') as pos, open(neg_fpath, 'w') as neg:
        for rec in pc_sites:
            iv = Interval(rec.position, rec.position)
            for mRNA in mRNA_trees[rec.seqid].query_all_overlaps(iv):

                if rec.strand != mRNA.strand:
                    continue

                xl_site = rec.position
                assert genome.get_sequence(rec.seqid, xl_site, xl_site, mRNA.strand) == 'T'
                found = False
                while not found:
                    ran_site = random.randint(mRNA.start + vp_bp, mRNA.end - vp_bp + 1)
                    if genome.get_sequence(rec.seqid, ran_site, ran_site, mRNA.strand) == 'T':
                        found = True
                seq_counter += 1

                skip = False
                for site, handle in [(xl_site, pos), (ran_site, neg)]:
                    if skip:
                        continue
                    if site - vp_bp < mRNA.start or site + vp_bp > mRNA.end:
                        skip = True
                        continue
                    seq_start = max(site - fl_bp, mRNA.start)
                    seq_end = min(site + fl_bp, mRNA.end)
                    seq = genome.get_sequence(rec.seqid, seq_start, seq_end, mRNA.strand)

                    site_in_seq = site - seq_start
                    if rec.strand == '-':
                        site_in_seq = len(seq) - 1 - site_in_seq
                    assert seq[site_in_seq] == 'T'

                    vp_start = site_in_seq - vp_bp
                    vp_end = site_in_seq + vp_bp

                    final_seq = (
                        seq[:vp_start].lower() +
                        seq[vp_start:vp_end + 1] +
                        seq[vp_end + 1:].lower()
                    )
                    print('>seq %i' % seq_counter, file=handle)
                    print(final_seq, file=handle)


if __name__ == '__main__':
    main()
