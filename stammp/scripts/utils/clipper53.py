import argparse


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('outfile')
    parser.add_argument('prime5_adapter')
    parser.add_argument('prime3_adapter')
    parser.add_argument('--clip_len', type=int, default=12)
    parser.add_argument('--min_len', type=int, default=1)
    parser.add_argument('--nt_barcode_5prime', type=int, default=0)
    parser.add_argument('--nt_barcode_3prime', type=int, default=0)
    parser.add_argument('--verbose', action='store_true')
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    clip_len = args.clip_len
    min_len = args.min_len
    bc5_len = args.nt_barcode_5prime
    bc3_len = args.nt_barcode_3prime

    prim5_string = args.prime5_adapter[-clip_len:]
    prim3_string = args.prime3_adapter[:clip_len]

    read_buffer = [''] * 4
    total_reads = 0
    discarded_reads = 0
    total_5prime_clipped = 0
    total_3prime_clipped = 0

    with open(args.infile) as infile, open(args.outfile, 'w') as outfile:
        for line_no, line in enumerate(infile):
            read_ind = line_no % 4
            read_buffer[read_ind] = line.rstrip()

            # just read a complete read
            if read_ind == 3:
                total_reads += 1
                nt_seq = read_buffer[1]

                left_hit = nt_seq.rfind(prim5_string)
                right_hit = nt_seq.find(prim3_string)

                if left_hit < 0:
                    read_start = bc5_len
                else:
                    read_start = left_hit + len(prim5_string) + bc5_len
                    total_5prime_clipped += 1

                if right_hit < 0:
                    read_end = len(nt_seq)
                else:
                    read_end = right_hit - bc3_len
                    total_3prime_clipped += 1

                if read_end - read_start >= min_len:
                    print(read_buffer[0], file=outfile)
                    print(read_buffer[1][read_start:read_end], file=outfile)
                    print(read_buffer[2], file=outfile)
                    print(read_buffer[3][read_start:read_end], file=outfile)
                else:
                    discarded_reads += 1

    if args.verbose:
        print('total reads:     %s' % total_reads)
        print('5prime clipped:  %s' % total_5prime_clipped)
        print('3prime clipped:  %s' % total_3prime_clipped)
        print('too short reads: %s' % discarded_reads)
        print('surviving reads: %.2f%%' % (100 - discarded_reads / total_reads * 100))


if __name__ == '__main__':
    main()
