import argparse
import os
import sys
import logging
import glob

from mockinbird.utils import argparse_helper as aph
from mockinbird.utils import pipeline as pl
from mockinbird.utils import module_utils as mu
from mockinbird import __version__
from mockinbird import LOG_DEFAULT_FORMAT, LOG_LEVEL_MAP


def register_arguments(parser):
    parser.add_argument('preprocess_dir', help='folder to files created by the preprocessing',
                        type=aph.dir_rx)
    prefix_help = ('preprocessing filename prefix - only required if there are multiple table '
                   'files in the specified preprocess directory')
    parser.add_argument('--prefix', help=prefix_help)
    no_pileup_help = 'do not require a pileup file and skip all tasks that depend on the pileup.'
    output_help = 'output directory - will be created if it does not exist'
    parser.add_argument('output_dir', help=output_help, type=aph.dir_rwx_create)
    config_help = 'path to the postprocessing config file'
    parser.add_argument('config_file', help=config_help, type=aph.file_r)
    log_level_help = 'verbosity level of the logger'
    parser.add_argument('--log_level', help=log_level_help, choices=LOG_LEVEL_MAP.keys(),
                        default='info')
    aph.add_version_arguments(parser)


def create_parser():
    description = 'run the PAR-CLIP postprocessing pipeline'
    parser = argparse.ArgumentParser(prog='mockinbird-postprocess', description=description)
    register_arguments(parser)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    run(args)


def run(args):
    # activate logging
    logging_file = os.path.join(args.output_dir, 'postprocess.log')
    logger = logging.getLogger()
    logger.setLevel(LOG_LEVEL_MAP[args.log_level])
    formatter = logging.Formatter(LOG_DEFAULT_FORMAT)

    console_handler = logging.StreamHandler(stream=sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    file_handler = logging.FileHandler(logging_file, mode='w')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    logger.info('mockinbird version: %s', __version__)
    logger.info('working directory: %s', os.getcwd())
    logger.info('started postprocessing via %r', ' '.join(sys.argv))

    required_ext = ['.table']

    scan_pat = os.path.join(args.preprocess_dir, '*' + required_ext[0])
    avail_prefixes = []
    for scfile in glob.glob(scan_pat):
        base_name = os.path.basename(scfile)
        prefix, ext = os.path.splitext(base_name)
        complete = True
        for req_ext in required_ext:
            pref_file = os.path.join(args.preprocess_dir, prefix + req_ext)
            if not os.path.exists(pref_file):
                complete = False
                break
        if complete:
            avail_prefixes.append(prefix)

    if len(avail_prefixes) == 0:
        logger.error('no %s files found in directory %s',
                     '|'.join(required_ext), args.preprocess_dir)
        sys.exit(1)

    if args.prefix is None:
        if len(avail_prefixes) != 1:
            logger.error('multiple table files found in directory %s',
                         args.preprocess_dir)
            logger.error('please use --prefix to select one of %s',
                         '|'.join(avail_prefixes))
            sys.exit(1)
        prefix, = avail_prefixes
    else:
        if args.prefix not in avail_prefixes:
            logger.error('no files for prefix %r found in directory %s',
                         args.prefix, args.preprocess_dir)
            logger.error('please use --prefix to select one of %s',
                         '|'.join(avail_prefixes))
            sys.exit(1)
        prefix = args.prefix

    qnormed_table = os.path.join(args.preprocess_dir, prefix + '.table')

    config = mu.parse_yaml(args.config_file)

    initial_files = {}
    if os.path.exists(qnormed_table):
        initial_files['table'] = qnormed_table

    general_cfg = {
        'general': {
            'output_dir': args.output_dir,
            'prefix': prefix,
        }
    }
    pipeline = pl.Pipeline(initial_files=initial_files, general_cfg=general_cfg,
                           cfg_path=args.config_file)

    mu.queue_pipeline(config, pipeline, def_lookup_path='mockinbird.utils.postprocess_modules')
    mu.run_pipeline(pipeline)

    for job in pipeline:
        job.cleanup()
    logger.info('All done. Bye!')


if __name__ == '__main__':
    main()
