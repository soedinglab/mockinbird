import argparse
import sys
from mockinbird.utils.parsers import PC_MANDATORY_FIELDS


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('table_file')
    parser.add_argument('output_file')
    args = parser.parse_args()

    with open(args.table_file) as table, open(args.output_file, 'w') as outfile:
        header = table.readline().split()
        if header[0] == 'chromosome':
            new_header = PC_MANDATORY_FIELDS + ['p_value']
            print(*new_header, sep='\t', file=outfile)
            for line in table:
                toks = line.split()
                pval = float(toks[4])
                toks[4] = 1 - pval
                toks.append(pval)
                print(*toks, sep='\t', file=outfile)
        else:
            print('Unknown table format', file=sys.stderr)
            sys.exit(1)


if __name__ == '__main__':
    main()
