import argparse
import configparser
import glob
import os
import sys

from stammp.utils import file_rw_or_dir_rwx, file_r
from stammp.utils import execute


def main():
    description = 'Tool to generate a new preprocessing config file'
    parser = argparse.ArgumentParser('stammp-generateConfig', description=description)
    output_help = 'folder in which the config file should be created'
    parser.add_argument('output', help=output_help, type=file_rw_or_dir_rwx)
    genome_help = 'path to the genome fasta file'
    parser.add_argument('genome_fasta', help=genome_help, type=file_r)
    bowtie_help = 'prefix of the bowtie index'
    parser.add_argument('--bowtie-index', help=bowtie_help)
    parser.add_argument('--create', action='store_true')
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
    # does a fasta index exist?
    fasta_ind = genome_path_abs + '.fai'
    if not os.path.isfile(fasta_ind):
        if args.create:
            cmd_args = [
                'samtools',
                'faidx',
                genome_path_abs,
            ]
            execute(cmd_args)
            print('[INFO] created %r' % fasta_ind)
        else:
            print('[WARNING]: %r does not exist' % (genome_path_abs + '.fai'),
                  file=sys.stderr)

    if args.bowtie_index:
        bowtie_index_abs = os.path.abspath(args.bowtie_index)
        bt_index_glob = "%s*" % bowtie_index_abs
        if len(glob.glob(bt_index_glob)) == 0:
            if args.create:
                index_dir = os.path.dirname(bt_index_glob)
                os.makedirs(index_dir)
                cmd_args = [
                    'bowtie-build',
                    args.genome_fasta,
                    args.bowtie_index,
                ]
                execute(cmd_args)
                print('[INFO] created bowtie index %r' % args.bowtie_index)
            else:
                print('[WARNING]: %r is not a prefix to an existing bowtie index'
                      % bowtie_index_abs, file=sys.stderr)
        cfg['general']['bowtieindex'] = bowtie_index_abs
    else:
        msg = (
            '[WARNING]: no bowtie index specified. '
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
