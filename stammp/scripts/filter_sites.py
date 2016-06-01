import argparse
from collections import Counter
from pyivtree import IVTree
from pyivtree import GenomicInterval as Interval
from stammp.utils import argparse_helper as aph


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('parclip_sites', type=aph.file_r)
    parser.add_argument('filtered_file', type=aph.file_rw_or_dir_rwx)
    parser.add_argument('gff_file', type=aph.file_r)
    parser.add_argument('--padding_bp', type=int, default=10)
    parser.add_argument('--filter_features', nargs='+', default=[])
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    filter_features = set(args.filter_features)

    gff_tree = IVTree()
    with open(args.gff_file) as gff:
        for line in gff:
            if line.startswith('#'):
                continue
            toks = line.split()
            try:
                chrom, _, annot, start_str, end_str, _, strand, *_ = toks
            except ValueError:
                continue
            iv = Interval(int(start_str) - args.padding_bp, int(end_str) + args.padding_bp)
            iv.strand = strand
            iv.annot = annot
            gff_tree.insert(iv)

    with open(args.parclip_sites) as pc, open(args.filtered_file, 'w') as out:
        header = pc.readline().split()
        print(*header, sep='\t', file=out)
        fil_annots = Counter()
        for line in pc:
            toks = line.split()
            chrom, pos_str, _, _, _, strand, _ = toks
            site_iv = Interval(int(pos_str), int(pos_str))
            has_overlap = False
            for ovl in gff_tree.query_all_overlaps(site_iv):
                if ovl.strand == strand:
                    if ovl.annot in filter_features:
                        has_overlap = True
                        break
            if not has_overlap:
                print(*toks, sep='\t', file=out)
            else:
                fil_annots[ovl.annot] += 1

    for annot, count in fil_annots.items():
        print('Removed %s PARCLIP sites annotated %r' % (count, annot))


if __name__ == '__main__':
    main()
