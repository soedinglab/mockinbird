import argparse
import sys
from stammp.obj import functions
from stammp.utils import ParclipSiteContainer
from stammp.utils import argparse_helper as aph


def create_parser():
    description = 'Set maximum occupancy to the specified quantile.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('inputfile', help='normalized PAR-CLIP file *.table', type=aph.file_r)
    parser.add_argument('outputfile', help='converted PAR-CLIP file *.table',
                        type=aph.file_rw_or_dir_rwx)
    parser.add_argument('--quantile', '-q', help='quantile [0, 1.0]', default=0.95,
                        type=float)
    return parser


def main(input_file, output_file, q):
    if not 0 <= q < 1:
        print('q must lie between 0 and 1 - got %s' % q)
        sys.exit(1)
    sites = ParclipSiteContainer.from_file(input_file)

    # dirty hack to avoid errors on empty files
    occ_vals = []
    for rec in sites:
        occ_vals.append(rec.occupancy)

    if len(occ_vals) > 0:
        max_occ = functions.getQuantile(occ_vals, q)

    records = []
    for rec in sites:
        if rec.occupancy > max_occ:
            rec = rec._replace(occupancy=max_occ)
        records.append(rec)

    new_sites = ParclipSiteContainer(records)
    new_sites.save2File(output_file)


def run():
    parser = create_parser()
    args = parser.parse_args()
    main(args.inputfile, args.outputfile, args.quantile)


if __name__ == '__main__':
    run()
