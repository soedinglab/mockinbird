import argparse
import os
import shutil

from stammp.utils import argparse_helper as aph


def create_parser():
    description = 'Tool to generate example config files'
    parser = argparse.ArgumentParser('stammp-generateConfig', description=description)
    output_help = 'output folder'
    parser.add_argument('output_directory', help=output_help, type=aph.dir_rwx_create)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    cur_dir = os.path.dirname(os.path.realpath(__file__))

    def generate_config(template, destination):
        if not os.path.exists(destination):
            shutil.copy(template, destination)
            print('[INFO] Successfully created %r' % destination)
        else:
            print('[WARNING] File %r already exists. I will not overwrite it.' % destination)

    preprocess_src = os.path.join(cur_dir, '../data/preprocess_template.cfg')
    preprocess_dest = os.path.join(args.output_directory, 'preprocess.cfg')
    generate_config(preprocess_src, preprocess_dest)

    postprocess_src = os.path.join(cur_dir, '../data/postprocess_template.cfg')
    postprocess_dest = os.path.join(args.output_directory, 'postprocess.cfg')
    generate_config(postprocess_src, postprocess_dest)

    print()
    print('All done. Please don\'t forget to adapt the config files to your needs!')


if __name__ == '__main__':
    main()
