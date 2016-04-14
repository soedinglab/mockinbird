import argparse
import configparser
import glob
import os
import sys

from stammp.utils import file_rw_or_dir_rwx, file_r


def main():
    description = 'Tool to generate a new preprocessing config file'
    parser = argparse.ArgumentParser('stammp-generateConfig', description=description)
    genome_help = 'path to the genome fasta file'
    parser.add_argument('genome_fasta', help=genome_help, type=file_r)
    output_help = 'folder in which the config file should be created'
    parser.add_argument('output', help=output_help, type=file_rw_or_dir_rwx)
    index_help = 'prefix of the genome index'
    parser.add_argument('--genome-index', help=index_help)
    args = parser.parse_args()

    cur_dir = os.path.dirname(os.path.realpath(__file__))
    config_path = os.path.join(cur_dir, 'preprocess.config')

    if os.path.isdir(args.output):
        out_path = os.path.join(args.output, 'stammp-preprocess.cfg')
    else:
        out_path = args.output
    cfg = configparser.ConfigParser(
        allow_no_value=True
    )
    cfg.read(config_path)

    genome_path_abs = os.path.abspath(args.genome_fasta)
    cfg['general']['genomefasta'] = genome_path_abs

    print()
    fasta_ind = genome_path_abs + '.fai'
    if not os.path.isfile(fasta_ind):
        print('[WARNING]: %r does not exist' % (genome_path_abs + '.fai'),
              file=sys.stderr)

    if args.genome_index:
        genome_index_abs = os.path.abspath(args.genome_index)
        genome_index_glob = "%s*" % genome_index_abs
        if len(glob.glob(genome_index_glob)) == 0:
            print('[WARNING]: %r is not a prefix to an existing genome index'
                  % genome_index_abs, file=sys.stderr)
        cfg['general']['genomeindex'] = genome_index_abs
    else:
        msg = (
            '[WARNING]: no genome index specified. '
            'You have to edit the config by hand and set this path.'
        )
        print(msg, file=sys.stderr)

    with open(out_path, 'w') as out_handle:
        # the config parser automatically ignores comments, so we have to copy
        # the header manually
        with open(config_path, 'r') as cfg_handle:
            for line in cfg_handle:
                line = line.strip()
                if line.startswith('#') or line == '':
                    print(line, file=out_handle)
                else:
                    break
        # write the config to the file
        cfg.write(out_handle)

    print('[SUCCESS] Successfully created %r' % out_path)
    print()
    print('You now may want to check if the default config matches your needs!')


if __name__ == '__main__':
    main()
