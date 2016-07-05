import argparse
import os
import configparser
import sys
import logging
import functools
import glob
from collections import OrderedDict
import operator

from stammp.utils.pipeline import CmdPipelineModule
from stammp.utils import config_validation as cv
from stammp.utils import argparse_helper as aph
from stammp.utils import pipeline as pl
from stammp import LOG_DEFAULT_FORMAT, LOG_LEVEL_MAP


def create_parser():
    description = 'run the PAR-CLIP postprocessing pipeline'
    parser = argparse.ArgumentParser(prog='stammp-postprocess', description=description)
    parser.add_argument('preprocess_dir', help='folder to files created by the preprocessing',
                        type=aph.dir_rx)
    prefix_help = ('preprocessing filename prefix - only required if there are multiple prefixes '
                   'in the specified preprocess directory')
    parser.add_argument('--prefix', help=prefix_help)
    no_pileup_help = 'do not require a pileup file and skip all tasks that depend on the pileup.'
    parser.add_argument('--no-pileup', help=no_pileup_help, action='store_true')
    output_help = 'output directory - will be created if it does not exist'
    parser.add_argument('output_dir', help=output_help, type=aph.dir_rwx_create)
    config_help = 'path to the postprocessing config file'
    parser.add_argument('config_file', help=config_help, type=aph.file_r)
    log_level_help = 'verbosity level of the logger'
    parser.add_argument('--log_level', help=log_level_help, choices=LOG_LEVEL_MAP.keys(),
                        default='info')
    aph.add_version_arguments(parser)
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

    required_ext = ['.qtable']
    if not args.no_pileup:
        required_ext.append('.mpileup')

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
        strict=False,
        inline_comment_prefixes=';',
        interpolation=configparser.ExtendedInterpolation(),
    )
    try:
        config.read(args.config_file)
    except configparser.Error as e:
        logger.error(e)
        sys.exit(1)

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
        'center_plot': OrderedDict([
            ('gff_file', cv.Annot(str, None, rel_file_r_validator)),
            ('output_prefix', cv.Annot(str, None, cv.id_converter)),
            ('downstream_bp', cv.Annot(int, 1000, cv.nonneg_integer)),
            ('upstream_bp', cv.Annot(int, 1000, cv.nonneg_integer)),
            ('gene_bp', cv.Annot(int, 750, cv.nonneg_integer)),
            ('min_trscr_size_bp', cv.Annot(int, 0, cv.nonneg_integer)),
            ('max_trscr_size_bp', cv.Annot(int, 5000, cv.nonneg_integer)),
            ('smoothing_window', cv.Annot(int, 20, cv.nonneg_integer)),
            ('labelCenterA', cv.Annot(str, None, cv.id_converter)),
            ('labelCenterB', cv.Annot(str, None, cv.id_converter)),
            ('labelBody', cv.Annot(str, None, cv.id_converter)),
            ('remove_tmp_files', cv.Annot(bool, True, cv.id_converter)),
        ]),
        'kmer_count_plot': OrderedDict([
            ('genome_fasta', cv.Annot(str, None, rel_file_r_validator)),
            ('output_prefix', cv.Annot(str, None, cv.id_converter)),
            ('kmer_k', cv.Annot(int, 3, cv.nonneg_integer)),
            ('first_index', cv.Annot(int, 0, cv.nonneg_integer)),
            ('last_index', cv.Annot(int, 1500, cv.nonneg_integer)),
            ('width', cv.Annot(int, 50, cv.nonneg_integer)),
            ('sort_key', cv.Annot(str, 'occ', sort_key_validator)),
            ('gff_exclude_path', cv.Annot(str, None, opt_file_validator)),
            ('gff_padding', cv.Annot(int, 20, cv.nonneg_integer)),
            ('remove_tmp_files', cv.Annot(bool, True, cv.id_converter)),
        ]),
        'kmer_logodd_plot': OrderedDict([
            ('genome_fasta', cv.Annot(str, None, rel_file_r_validator)),
            ('output_prefix', cv.Annot(str, None, cv.id_converter)),
            ('kmer_k', cv.Annot(int, 3, cv.nonneg_integer)),
            ('sort_key', cv.Annot(str, 'occ', sort_key_validator)),
            ('gff_exclude_path', cv.Annot(str, None, opt_file_validator)),
            ('use_quantiles', cv.Annot(bool, True, cv.id_converter)),
            ('negative_set_gff', cv.Annot(str, None, rel_file_r_validator)),
            ('n_negative_seqs', cv.Annot(int, 20000, cv.nonneg_integer)),
        ]),
        'xxmotif': OrderedDict([
            ('genome_fasta', cv.Annot(str, None, rel_file_r_validator)),
            ('output_prefix', cv.Annot(str, None, cv.id_converter)),
            ('negative_set_gff', cv.Annot(str, None, rel_file_r_validator)),
            ('n_negative_seqs', cv.Annot(int, 20000, cv.nonneg_integer)),
            ('plot_top_n_pwm', cv.Annot(int, 3, cv.nonneg_integer)),
            ('first_index', cv.Annot(int, 0, cv.nonneg_integer)),
            ('last_index', cv.Annot(int, 1500, cv.nonneg_integer)),
            ('width', cv.Annot(int, 12, cv.nonneg_integer)),
            ('sort_key', cv.Annot(str, 'occ', sort_key_validator)),
            ('gff_exclude_path', cv.Annot(str, None, opt_file_validator)),
            ('gff_padding', cv.Annot(int, 20, cv.nonneg_integer)),
            ('remove_tmp_files', cv.Annot(bool, True, cv.id_converter)),
        ]),
        'tr_freq_plot': OrderedDict([
            ('output_prefix', cv.Annot(str, None, cv.id_converter)),
            ('min_cov', cv.Annot(int, 5, cv.nonneg_integer)),
            ('y_axis_limit', cv.Annot(float, 0, cv.id_converter)),
            ('remove_tmp_files', cv.Annot(bool, True, cv.id_converter)),
        ]),
        'heatmap_plot': OrderedDict([
            ('gff_file', cv.Annot(str, None, rel_file_r_validator)),
            ('output_prefix', cv.Annot(str, None, cv.id_converter)),
            ('downstream_bp', cv.Annot(int, 4000, cv.nonneg_integer)),
            ('upstream_bp', cv.Annot(int, 1000, cv.nonneg_integer)),
            ('min_trscr_size_bp', cv.Annot(int, 0, cv.nonneg_integer)),
            ('max_trscr_size_bp', cv.Annot(int, 5000, cv.nonneg_integer)),
            ('xbins', cv.Annot(int, 500, cv.nonneg_integer)),
            ('ybins', cv.Annot(int, 500, cv.nonneg_integer)),
            ('x_pixels', cv.Annot(int, 500, cv.nonneg_integer)),
            ('y_pixels', cv.Annot(int, 500, cv.nonneg_integer)),
            ('remove_tmp_files', cv.Annot(bool, True, cv.id_converter)),
        ]),
        'heatmap_small_plot': OrderedDict([
            ('gff_file', cv.Annot(str, None, rel_file_r_validator)),
            ('output_prefix', cv.Annot(str, None, cv.id_converter)),
            ('downstream_bp', cv.Annot(int, 500, cv.nonneg_integer)),
            ('upstream_bp', cv.Annot(int, 1000, cv.nonneg_integer)),
            ('min_trscr_size_bp', cv.Annot(int, 0, cv.nonneg_integer)),
            ('max_trscr_size_bp', cv.Annot(int, 5000, cv.nonneg_integer)),
            ('xbins', cv.Annot(int, 500, cv.nonneg_integer)),
            ('ybins', cv.Annot(int, 500, cv.nonneg_integer)),
            ('x_pixels', cv.Annot(int, 500, cv.nonneg_integer)),
            ('y_pixels', cv.Annot(int, 500, cv.nonneg_integer)),
            ('remove_tmp_files', cv.Annot(bool, True, cv.id_converter)),
        ]),
        'gff_filter': OrderedDict([
            ('file_postfix', cv.Annot(str, 'fil', cv.id_converter)),
            ('padding_bp', cv.Annot(int, 10, cv.nonneg_integer)),
            ('features', cv.Annot(str, [], cv.comma_sep_args)),
            ('filter_gff', cv.Annot(str, None, rel_file_r_validator)),
        ]),
        'ss_indicator': OrderedDict([
            ('genome_fasta', cv.Annot(str, None, rel_file_r_validator)),
            ('output_prefix', cv.Annot(str, None, cv.id_converter)),
            ('sample_gff', cv.Annot(str, None, rel_file_r_validator)),
            ('first_index', cv.Annot(int, 0, cv.nonneg_integer)),
            ('last_index', cv.Annot(int, 1500, cv.nonneg_integer)),
            ('width', cv.Annot(int, 12, cv.nonneg_integer)),
            ('sort_key', cv.Annot(str, 'occ', sort_key_validator)),
            ('remove_tmp_files', cv.Annot(bool, True, cv.id_converter)),
        ]),
    }

    output_getter = operator.attrgetter('cur_output')
    modules = {
        'center_plot': (CenterPlotModule, output_getter),
        'kmer_count_plot': (KmerPerPositionModule, output_getter),
        'kmer_logodd_plot': (KmerLogoddModule, output_getter),
        'xxmotif': (XXmotifModule, output_getter),
        'tr_freq_plot': (TransitionFrequencyModule, lambda x: pileup_file),
        'heatmap_plot': (HeatmapPlotModule, output_getter),
        'heatmap_small_plot': (HeatmapSmallPlotModule, output_getter),
        'gff_filter': (GffFilterModule, output_getter),
        'ss_indicator': (SSIndicatorModule, output_getter),
    }

    try:
        cfg_dicts = cv.opt_config(config, cfg_format, 'type')
    except cv.ConfigError:
        sys.exit(1)

    fil_cfg_dicts = []
    if args.no_pileup:
        for cfg_dict in cfg_dicts:
            if cfg_dict['type'] == 'tr_freq_plot':
                logger.warn('Scheduled transition frequency plot conflicts with '
                            'the option %r', '--no-pileup')
            else:
                fil_cfg_dicts.append(cfg_dict)
        cfg_dicts = fil_cfg_dicts

    pipeline = pl.Pipeline(qnormed_table)
    for cfg_dict in cfg_dicts:
        try:
            module_class, input_factory = modules[cfg_dict['type']]
        except KeyError:
            logger.warn('No module found for type %r. Skipping', cfg_dict['type'])
            continue
        module = module_class()
        infile = input_factory(pipeline)
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
                    msg = 'Additional Output:\n\n'
                    logger.info(msg + stdout)
    for job in pipeline:
        job.cleanup()
    logger.info('All done. Bye!')


class CenterPlotModule(CmdPipelineModule):

    def prepare(self, norm_table_file, outdir, prefix, cfg):
        cmd = [
            'stammp-makeCenterBothEnds',
            '%r' % norm_table_file,
            '%r' % outdir,
            '%r' % cfg['output_prefix'],
            '%r' % cfg['gff_file'],
            '-d %s' % cfg['downstream_bp'],
            '-u %s' % cfg['upstream_bp'],
            '--min %s' % cfg['min_trscr_size_bp'],
            '--max %s' % cfg['max_trscr_size_bp'],
            '-g %s' % cfg['gene_bp'],
            '--plotSmooth %s' % cfg['smoothing_window'],
            '--labelCenterA %r' % cfg['labelCenterA'],
            '--labelCenterB %r' % cfg['labelCenterB'],
            '--labelBody %r' % cfg['labelBody'],
            '--title %r' % '[%s] %s' % (cfg['output_prefix'], prefix),
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
        lo_negset_dir = os.path.join(outdir, 'lo_negset')
        width = 15
        self._data['files'].append(lo_negset_dir)
        negset_cmd = [
            'stammp-makeNegSets',
            '--number %s' % cfg['n_negative_seqs'],
            '%r' % cfg['negative_set_gff'],
            '%r' % cfg['genome_fasta'],
            '%r' % cfg['output_prefix'],
            '--width %s' % width,
            lo_negset_dir,
        ]
        self._cmds.append(negset_cmd)

        fmt_args = cfg['output_prefix'], cfg['n_negative_seqs'], width, cfg['kmer_k']
        negset_fname = 'rnd_sequences_%s_%s_w%s_%smer.table' % fmt_args
        negset_file = os.path.join(lo_negset_dir, negset_fname)
        cmd = [
            'stammp-makeKmerLogOdds',
            '%r' % norm_table_file,
            '%r' % outdir,
            '%r' % cfg['output_prefix'],
            '%r' % cfg['genome_fasta'],
            '%r' % negset_file,
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

        xx_negset_dir = os.path.join(outdir, 'xx_negset')
        self._data['files'].append(xx_negset_dir)
        negset_cmd = [
            'stammp-makeNegSets',
            '--number %s' % cfg['n_negative_seqs'],
            '--width %s ' % cfg['width'],
            '%r' % cfg['negative_set_gff'],
            '%r' % cfg['genome_fasta'],
            '%r' % cfg['output_prefix'],
            xx_negset_dir,
        ]
        self._cmds.append(negset_cmd)

        fmt_args = cfg['output_prefix'], cfg['n_negative_seqs'], cfg['width']
        negset_fname = 'rnd_sequences_%s_%s_w%s.fa' % fmt_args
        negset_file = os.path.join(xx_negset_dir, negset_fname)
        cmd = [
            'stammp-xxmotif',
            norm_table_file,
            '%r' % cfg['genome_fasta'],
            '%r' % outdir,
            '%r' % cfg['output_prefix'],
            '--plotPWM %s' % cfg['plot_top_n_pwm'],
            '--start %s' % cfg['first_index'],
            '--stop %s' % cfg['last_index'],
            '--width %s' % cfg['width'],
            '--key %s' % cfg['sort_key'],
            '--awidth %s' % cfg['gff_padding'],
            '--negSet %r' % negset_file,
        ]
        if cfg['gff_exclude_path']:
            cmd.append('--filterGFF %r' % cfg['gff_exclude_path'])
        if not cfg['remove_tmp_files']:
            cmd.append('--keep-tmp-files')
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


class GffFilterModule(CmdPipelineModule):

    def prepare(self, sites_file, outdir, prefix, cfg):
        file_name = os.path.basename(sites_file)
        bname, ext = os.path.splitext(file_name)
        out_bname = '%s_%s' % (bname, cfg['file_postfix'])
        out_file = os.path.join(outdir, out_bname + ext)
        self._data['output'] = out_file
        cmd = [
            'stammp-gffFilterSites',
            sites_file,
            out_file,
            cfg['filter_gff'],
            '--padding_bp %s' % cfg['padding_bp'],
        ]
        if len(cfg['features']) > 0:
            cmd.append('--filter_features')
            for feature in cfg['features']:
                cmd.append('%r' % feature)
        self._cmds.append(cmd)


class SSIndicatorModule(CmdPipelineModule):

    def prepare(self, sites_file, outdir, prefix, cfg):
        cmd = [
            'stammp-ss_indicator',
            '%r' % sites_file,
            '%r' % cfg['genome_fasta'],
            '%r' % cfg['sample_gff'],
            '%r' % outdir,
            '--start %s' % cfg['first_index'],
            '--stop %s' % cfg['last_index'],
            '--width %s' % cfg['width'],
            '--key %r' % cfg['sort_key'],
            '--prefix %r' % cfg['output_prefix'],
        ]
        if not cfg['remove_tmp_files']:
            cmd.append('--keep-tmp-files')
        self._cmds.append(cmd)


if __name__ == '__main__':
    main()
