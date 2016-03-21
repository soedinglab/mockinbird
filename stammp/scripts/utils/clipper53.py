import argparse
import re


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('outfile')
    parser.add_argument('prime5_adapter')
    parser.add_argument('prime3_adapter')
    parser.add_argument('--aggressive', action='store_true')
    parser.add_argument('--clip_len', type=int, default=12)
    parser.add_argument('--min_len', type=int, default=0)
    parser.add_argument('--nt_barcode_5prime', type=int, default=0)
    parser.add_argument('--verbose', action='store_true')
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    prim5_pat = re.compile(args.prime5_adapter[-args.clip_len:])
    prim3_pat = re.compile(args.prime3_adapter[:args.clip_len])
    aggressive = args.aggressive
    clip_len = args.clip_len
    min_len = args.min_len
    bc_len = args.nt_barcode_5prime
    if aggressive:
        # all possible suffixes of the search sequences
        prim3_pset = set()
        for i in range(1, clip_len):
            prim3_pset.add(args.prime3_adapter[:i])

    read_buffer = [''] * 4
    total_reads = 0
    discarded_reads = 0
    with open(args.infile) as infile, open(args.outfile, 'w') as outfile:
        for line_no, line in enumerate(infile):
            read_ind = line_no % 4
            read_buffer[read_ind] = line.rstrip()

            # just read a complete read
            if read_ind == 3:
                total_reads += 1
                nt_seq = read_buffer[1]
                read_start = 0
                read_end = len(nt_seq)

                prim5_matches = list(prim5_pat.finditer(nt_seq))
                if len(prim5_matches) > 0:
                    read_start = prim5_matches[-1].end()
                    read_start += bc_len

                prim3_matches = list(prim3_pat.finditer(nt_seq))
                if len(prim3_matches) > 0:
                    read_end = prim3_matches[0].start()
                else:
                    if aggressive:
                        prime3_trim = False
                        for i in range(clip_len - 1, 0, -1):
                            if nt_seq[-i:] in prim3_pset:
                                prime3_trim = True
                                break
                        if prime3_trim:
                            read_end -= i

                if read_end - read_start >= min_len:
                    print(read_buffer[0], file=outfile)
                    print(read_buffer[1][read_start:read_end], file=outfile)
                    print(read_buffer[2], file=outfile)
                    print(read_buffer[3][read_start:read_end], file=outfile)
                else:
                    discarded_reads += 1

    if args.verbose:
        print('total reads:     %s' % total_reads)
        print('discarded reads: %s' % discarded_reads)
        print('surviving reads: %.2f%%' % (100 - discarded_reads / total_reads * 100))


if __name__ == '__main__':
    main()
