import argparse
import os
import configparser
import sys
import logging
import functools
import glob

from stammp.utils.pipeline import CmdPipelineModule
from stammp.utils import config_validation as cv
from stammp.utils import argparse_helper as aph
from stammp.utils import pipeline as pl
from stammp import LOG_DEFAULT_FORMAT, LOG_LEVEL_MAP


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('preprocess_dir', type=aph.dir_rx)
    parser.add_argument('--prefix')
    parser.add_argument('output_dir', type=aph.dir_rwx)
    parser.add_argument('config_file', type=aph.file_r)
    parser.add_argument('--log_level', choices=LOG_LEVEL_MAP.keys(), default='info')
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

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

    logger.info('started postprocessing via %r', ' '.join(sys.argv))

    required_ext = ['.mpileup', '.qtable']
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
        logger.error('no complete set of %s files found in directory %s',
                     '|'.join(required_ext), args.preprocess_dir)
        sys.exit(1)

    if args.prefix is None:
        if len(avail_prefixes) != 1:
            logger.error('multiple sets of files found in directory %s',
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

    qnormed_table = os.path.join(args.preprocess_dir, prefix + '.qtable')
    pileup_file = os.path.join(args.preprocess_dir, prefix + '.mpileup')

    # config parsing
    config = configparser.ConfigParser(
        inline_comment_prefixes=';',
        interpolation=configparser.ExtendedInterpolation(),
    )
    config.read(args.config_file)

    def rel_file_r_validator(path):
        if not os.path.isabs(path):
            parent_path = os.path.dirname(args.config_file)
            path = os.path.join(os.path.abspath(parent_path), path)
        return cv.file_r_validator(path)

    sort_keys = ['occ', 'm', 'r', 'mr', 'pvalue']
    sort_key_validator = functools.partial(cv.in_set_validator, item_set=sort_keys)

    def opt_file_validator(path):
        if path.strip() == '':
            return None
        return rel_file_r_validator(path)

    cfg_format = {
        'center_plot': {
            'gff_file': cv.Annot(str, None, rel_file_r_validator),
            'output_prefix': cv.Annot(str, None, cv.id_converter),
            'downstream_bp': cv.Annot(int, 1000, cv.nonneg_integer),
            'upstream_bp': cv.Annot(int, 1000, cv.nonneg_integer),
            'gene_bp': cv.Annot(int, 750, cv.nonneg_integer),
            'min_trscr_size_bp': cv.Annot(int, 0, cv.nonneg_integer),
            'max_trscr_size_bp': cv.Annot(int, 5000, cv.nonneg_integer),
            'smoothing_window': cv.Annot(int, 20, cv.nonneg_integer),
            'labelCenterA': cv.Annot(str, None, cv.id_converter),
            'labelCenterB': cv.Annot(str, None, cv.id_converter),
            'labelBody': cv.Annot(str, None, cv.id_converter),
            'remove_tmp_files': cv.Annot(bool, True, cv.id_converter),
        },
        'kmer_count_plot': {
            'genome_fasta': cv.Annot(str, None, rel_file_r_validator),
            'output_prefix': cv.Annot(str, None, cv.id_converter),
            'kmer_k': cv.Annot(int, 3, cv.nonneg_integer),
            'first_index': cv.Annot(int, 0, cv.nonneg_integer),
            'last_index': cv.Annot(int, 1500, cv.nonneg_integer),
            'width': cv.Annot(int, 50, cv.nonneg_integer),
            'sort_key': cv.Annot(str, 'occ', sort_key_validator),
            'gff_exclude_path': cv.Annot(str, None, opt_file_validator),
            'gff_padding': cv.Annot(int, 20, cv.nonneg_integer),
            'remove_tmp_files': cv.Annot(bool, True, cv.id_converter),
        },
        'kmer_logodd_plot': {
            'genome_fasta': cv.Annot(str, None, rel_file_r_validator),
            'output_prefix': cv.Annot(str, None, cv.id_converter),
            'kmer_k': cv.Annot(int, 3, cv.nonneg_integer),
            'sort_key': cv.Annot(str, 'occ', sort_key_validator),
            'gff_exclude_path': cv.Annot(str, None, opt_file_validator),
            'use_quantiles': cv.Annot(bool, True, cv.id_converter),
            'negative_set': cv.Annot(str, None, rel_file_r_validator),
        },
        'xxmotif': {
            'genome_fasta': cv.Annot(str, None, rel_file_r_validator),
            'output_prefix': cv.Annot(str, None, cv.id_converter),
            'negative_set': cv.Annot(str, '', opt_file_validator),
            'plot_top_n_pwm': cv.Annot(int, 3, cv.nonneg_integer),
            'first_index': cv.Annot(int, 0, cv.nonneg_integer),
            'last_index': cv.Annot(int, 1500, cv.nonneg_integer),
            'width': cv.Annot(int, 12, cv.nonneg_integer),
            'sort_key': cv.Annot(str, 'occ', sort_key_validator),
            'gff_exclude_path': cv.Annot(str, None, opt_file_validator),
            'gff_padding': cv.Annot(int, 20, cv.nonneg_integer),
        },
        'tr_freq_plot': {
            'output_prefix': cv.Annot(str, None, cv.id_converter),
            'min_cov': cv.Annot(int, 5, cv.nonneg_integer),
            'y_axis_limit': cv.Annot(float, 0, cv.id_converter),
            'remove_tmp_files': cv.Annot(bool, True, cv.id_converter),
        },
        'heatmap_plot': {
            'gff_file': cv.Annot(str, None, rel_file_r_validator),
            'output_prefix': cv.Annot(str, None, cv.id_converter),
            'downstream_bp': cv.Annot(int, 4000, cv.nonneg_integer),
            'upstream_bp': cv.Annot(int, 1000, cv.nonneg_integer),
            'min_trscr_size_bp': cv.Annot(int, 0, cv.nonneg_integer),
            'max_trscr_size_bp': cv.Annot(int, 5000, cv.nonneg_integer),
            'xbins': cv.Annot(int, 500, cv.nonneg_integer),
            'ybins': cv.Annot(int, 500, cv.nonneg_integer),
            'x_pixels': cv.Annot(int, 500, cv.nonneg_integer),
            'y_pixels': cv.Annot(int, 500, cv.nonneg_integer),
            'remove_tmp_files': cv.Annot(bool, True, cv.id_converter),
        },
        'heatmap_small_plot': {
            'gff_file': cv.Annot(str, None, rel_file_r_validator),
            'output_prefix': cv.Annot(str, None, cv.id_converter),
            'downstream_bp': cv.Annot(int, 500, cv.nonneg_integer),
            'upstream_bp': cv.Annot(int, 1000, cv.nonneg_integer),
            'min_trscr_size_bp': cv.Annot(int, 0, cv.nonneg_integer),
            'max_trscr_size_bp': cv.Annot(int, 5000, cv.nonneg_integer),
            'xbins': cv.Annot(int, 500, cv.nonneg_integer),
            'ybins': cv.Annot(int, 500, cv.nonneg_integer),
            'x_pixels': cv.Annot(int, 500, cv.nonneg_integer),
            'y_pixels': cv.Annot(int, 500, cv.nonneg_integer),
            'remove_tmp_files': cv.Annot(bool, True, cv.id_converter),
        },
    }

    modules = {
        'center_plot': (CenterPlotModule, qnormed_table),
        'kmer_count_plot': (KmerPerPositionModule, qnormed_table),
        'kmer_logodd_plot': (KmerLogoddModule, qnormed_table),
        'xxmotif': (XXmotifModule, qnormed_table),
        'tr_freq_plot': (TransitionFrequencyModule, pileup_file),
        'heatmap_plot': (HeatmapPlotModule, qnormed_table),
        'heatmap_small_plot': (HeatmapSmallPlotModule, qnormed_table),
    }

    try:
        cfg_dicts = cv.opt_config(config, cfg_format, 'type')
    except cv.ConfigError:
        sys.exit(1)

    pipeline = pl.Pipeline(qnormed_table)
    for cfg_dict in cfg_dicts:
        module_class, infile = modules[cfg_dict['type']]
        module = module_class()
        module.prepare(infile, args.output_dir, prefix, cfg_dict)
        module.msg = 'Executing work package %r' % cfg_dict['sec_name']
        pipeline.schedule(module)

    for job in pipeline:
        if hasattr(job, 'msg'):
            logger.info(job.msg)
        res = job.execute()
        if res:
            for stdout, stderr in res:
                if stdout.strip() != '':
                    logger.info('\n\n' + stdout)
    logger.info('All done. Bye!')


class CenterPlotModule(CmdPipelineModule):

    def prepare(self, norm_table_file, outdir, prefix, cfg):
        cmd = [
            'stammp-makeCenterBothEnds',
            '%r' % norm_table_file,
            '%r' % outdir,
            cfg['output_prefix'],
            cfg['gff_file'],
            '-d %s' % cfg['downstream_bp'],
            '-u %s' % cfg['upstream_bp'],
            '--min %s' % cfg['min_trscr_size_bp'],
            '--max %s' % cfg['max_trscr_size_bp'],
            '-g %s' % cfg['gene_bp'],
            '--plotSmooth %s' % cfg['smoothing_window'],
            '--labelCenterA %r' % cfg['labelCenterA'],
            '--labelCenterB %r' % cfg['labelCenterB'],
            '--labelBody %r' % cfg['labelBody'],
        ]
        if cfg['remove_tmp_files']:
            cmd.append('-r')
        self._cmds.append(cmd)


class KmerPerPositionModule(CmdPipelineModule):

    def prepare(self, norm_table_file, outdir, prefix, cfg):
        cmd = [
            'stammp-makeKmerPerPosition',
            norm_table_file,
            cfg['genome_fasta'],
            outdir,
            cfg['output_prefix'],
            '--kmer %s' % cfg['kmer_k'],
            '--start %s' % cfg['first_index'],
            '--stop %s' % cfg['last_index'],
            '--width %s' % cfg['width'],
            '--key %s' % cfg['sort_key'],
            '--awidth %s' % cfg['gff_padding'],
        ]
        if cfg['gff_exclude_path']:
            cmd.append('--filterGFF %s' % cfg['gff_exclude_path'])
        if cfg['remove_tmp_files']:
            cmd.append('-r')
        self._cmds.append(cmd)


class KmerLogoddModule(CmdPipelineModule):

    def prepare(self, norm_table_file, outdir, prefix, cfg):
        cmd = [
            'stammp-makeKmerLogOdds',
            '%r' % norm_table_file,
            '%r' % outdir,
            '%r' % cfg['output_prefix'],
            '%r' % cfg['genome_fasta'],
            '%r' % cfg['negative_set'],
            '--kmer %s' % cfg['kmer_k'],
            '--key %s' % cfg['sort_key'],
        ]
        if cfg['gff_exclude_path']:
            cmd.append('--filterGFF %s' % cfg['gff_exclude_path'])
        if cfg['use_quantiles']:
            cmd.append('-q')
        self._cmds.append(cmd)


class XXmotifModule(CmdPipelineModule):

    def prepare(self, norm_table_file, outdir, prefix, cfg):

        cmd = [
            'stammp-xxmotif',
            norm_table_file,
            cfg['genome_fasta'],
            outdir,
            cfg['output_prefix'],
            '--plotPWM %s' % cfg['plot_top_n_pwm'],
            '--start %s' % cfg['first_index'],
            '--stop %s' % cfg['last_index'],
            '--width %s' % cfg['width'],
            '--key %s' % cfg['sort_key'],
            '--awidth %s' % cfg['gff_padding'],
        ]
        if cfg['gff_exclude_path']:
            cmd.append('--filterGFF %s' % cfg['gff_exclude_path'])
        if cfg['negative_set']:
            cmd.append('--negSet %s' % cfg['negative_set'])
        self._cmds.append(cmd)


class TransitionFrequencyModule(CmdPipelineModule):

    def prepare(self, pileup_file, outdir, prefix, cfg):

        cmd = [
            'stammp-makeNucleotideProbabilities',
            pileup_file,
            outdir,
            cfg['output_prefix'],
            '-c %s' % cfg['min_cov'],
            '-l %s' % cfg['y_axis_limit'],
        ]
        if cfg['remove_tmp_files']:
            cmd.append('-r')
        self._cmds.append(cmd)


class HeatmapPlotModule(CmdPipelineModule):

    def prepare(self, norm_table_file, outdir, prefix, cfg):
        cmd = [
            'stammp-makeHeatMap',
            norm_table_file,
            outdir,
            cfg['output_prefix'],
            cfg['gff_file'],
            '-d %s' % cfg['downstream_bp'],
            '-u %s' % cfg['upstream_bp'],
            '--min %s' % cfg['min_trscr_size_bp'],
            '--max %s' % cfg['max_trscr_size_bp'],
            '--xbins %s' % cfg['xbins'],
            '--ybins %s' % cfg['ybins'],
            '--xpx %s' % cfg['x_pixels'],
            '--ypx %s' % cfg['y_pixels'],
        ]
        if cfg['remove_tmp_files']:
            cmd.append('-r')
        self._cmds.append(cmd)


class HeatmapSmallPlotModule(CmdPipelineModule):

    def prepare(self, norm_table_file, outdir, prefix, cfg):
        cmd = [
            'stammp-makeHeatMapSmall',
            norm_table_file,
            outdir,
            cfg['output_prefix'],
            cfg['gff_file'],
            '-d %s' % cfg['downstream_bp'],
            '-u %s' % cfg['upstream_bp'],
            '--min %s' % cfg['min_trscr_size_bp'],
            '--max %s' % cfg['max_trscr_size_bp'],
            '--xbins %s' % cfg['xbins'],
            '--ybins %s' % cfg['ybins'],
            '--xpx %s' % cfg['x_pixels'],
            '--ypx %s' % cfg['y_pixels'],
        ]
        if cfg['remove_tmp_files']:
            cmd.append('-r')
        self._cmds.append(cmd)
