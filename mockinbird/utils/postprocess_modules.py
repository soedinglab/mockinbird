import mockinbird.utils.pipeline as pl
import mockinbird.utils.config_validation as cv

import os
from functools import partial

sort_keys = ['occupancy', 'transitions', 'coverage', 'score']
sort_key_validator = partial(cv.in_set_validator, item_set=sort_keys)


def opt_file_validator(path, cfg_path):
    if path.strip() == '':
        return None
    return cv.rel_file_r_validator(path, cfg_path)


class CenterPlotModule(pl.CmdPipelineModule):

    def __init__(self, pipeline):
        relpath_conv = partial(cv.rel_file_r_validator, cfg_path=pipeline.cfg_path)
        cfg_fmt = [
            ('gff_file', cv.Annot(str, converter=relpath_conv)),
            ('output_prefix', cv.Annot(str)),
            ('downstream_bp', cv.Annot(int, default=1000, converter=cv.nonneg_integer)),
            ('upstream_bp', cv.Annot(int, default=1000, converter=cv.nonneg_integer)),
            ('gene_bp', cv.Annot(int, default=750, converter=cv.nonneg_integer)),
            ('min_trscr_size_bp', cv.Annot(int, default=1500, converter=cv.nonneg_integer)),
            ('max_trscr_size_bp', cv.Annot(int, default=100000, converter=cv.nonneg_integer)),
            ('smoothing_window', cv.Annot(int, default=20, converter=cv.nonneg_integer)),
            ('labelCenterA', cv.Annot(str)),
            ('labelCenterB', cv.Annot(str)),
            ('labelBody', cv.Annot(str)),
            ('remove_tmp_files', cv.Annot(bool, default=True)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']
        table_file = pipeline.get_curfile(fmt='table')

        cmd = [
            'mb-plot-metagene-nobs',
            '%r' % table_file,
            '%r' % output_dir,
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


class CenterPlotBSModule(pl.CmdPipelineModule):
    def __init__(self, pipeline):
        relpath_conv = partial(cv.rel_file_r_validator, cfg_path=pipeline.cfg_path)
        cfg_fmt = [
            ('gff_file', cv.Annot(str, converter=relpath_conv)),
            ('output_prefix', cv.Annot(str)),
            ('downstream_bp', cv.Annot(int, default=1000, converter=cv.nonneg_integer)),
            ('upstream_bp', cv.Annot(int, default=1000, converter=cv.nonneg_integer)),
            ('gene_bp', cv.Annot(int, default=750, converter=cv.nonneg_integer)),
            ('min_trscr_size_bp', cv.Annot(int, default=1500, converter=cv.nonneg_integer)),
            ('max_trscr_size_bp', cv.Annot(int, default=100000, converter=cv.nonneg_integer)),
            ('smoothing_window', cv.Annot(int, default=20, converter=cv.nonneg_integer)),
            ('labelCenterA', cv.Annot(str)),
            ('labelCenterB', cv.Annot(str)),
            ('labelBody', cv.Annot(str)),
            ('remove_tmp_files', cv.Annot(bool, default=True)),
            ('bootstrap_iter', cv.Annot(int, default=2500, converter=cv.nonneg_integer)),
            ('n_processes', cv.Annot(int, default=4, converter=cv.nonneg_integer)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']
        table_file = pipeline.get_curfile(fmt='table')

        cmd = [
            'mb-plot-metagene',
            '%r' % table_file,
            '%r' % output_dir,
            '%r' % cfg['output_prefix'],
            '%r' % cfg['gff_file'],
            '--downstream_bp %s' % cfg['downstream_bp'],
            '--upstream_bp %s' % cfg['upstream_bp'],
            '--min_ts_len %s' % cfg['min_trscr_size_bp'],
            '--max_ts_len %s' % cfg['max_trscr_size_bp'],
            '--gene_bp %s' % cfg['gene_bp'],
            '--smooth_window %s' % cfg['smoothing_window'],
            '--labelCenterA %r' % cfg['labelCenterA'],
            '--labelCenterB %r' % cfg['labelCenterB'],
            '--labelBody %r' % cfg['labelBody'],
            '--title %r' % '[%s] %s' % (cfg['output_prefix'], prefix),
            '--n_bs_iterations %s' % cfg['bootstrap_iter'],
            '--n_processes %s' % cfg['n_processes'],
        ]
        if cfg['remove_tmp_files']:
            cmd.append('--cleanup')
        self._cmds.append(cmd)


class KmerPerPositionModule(pl.CmdPipelineModule):
    def __init__(self, pipeline):
        optpath_conv = partial(opt_file_validator, cfg_path=pipeline.cfg_path)
        genome_conv = partial(cv.rel_genome_validator, cfg_path=pipeline.cfg_path)
        cfg_fmt = [
            ('genome_fasta', cv.Annot(str, converter=genome_conv)),
            ('output_prefix', cv.Annot(str)),
            ('kmer_k', cv.Annot(int, default=3, converter=cv.nonneg_integer)),
            ('first_index', cv.Annot(int, default=0, converter=cv.nonneg_integer)),
            ('last_index', cv.Annot(int, default=1500, converter=cv.nonneg_integer)),
            ('width', cv.Annot(int, default=50, converter=cv.nonneg_integer)),
            ('sort_key', cv.Annot(str, default='occupancy', converter=sort_key_validator)),
            ('gff_exclude_path', cv.Annot(str, default='', converter=optpath_conv,
                                          warn_if_missing=False)),
            ('gff_padding', cv.Annot(int, default=20, converter=cv.nonneg_integer)),
            ('remove_tmp_files', cv.Annot(bool, default=True)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']

        table_file = pipeline.get_curfile(fmt='table')
        cmd = [
            'mb-plot-kmer-enrichment',
            '%r' % table_file,
            '%r' % cfg['genome_fasta'],
            '%r' % output_dir,
            '%r' % cfg['output_prefix'],
            '--kmer %s' % cfg['kmer_k'],
            '--start %s' % cfg['first_index'],
            '--stop %s' % cfg['last_index'],
            '--width %s' % cfg['width'],
            '--key %s' % cfg['sort_key'],
            '--awidth %s' % cfg['gff_padding'],
        ]
        if cfg['gff_exclude_path']:
            cmd.append('--filterGFF %r' % cfg['gff_exclude_path'])

        if cfg['remove_tmp_files'] and not cfg['keep_all']:
            cmd.append('-r')
        self._cmds.append(cmd)


class KmerLogoddModule(pl.CmdPipelineModule):

    def __init__(self, pipeline):
        optpath_conv = partial(opt_file_validator, cfg_path=pipeline.cfg_path)
        relpath_conv = partial(cv.rel_file_r_validator, cfg_path=pipeline.cfg_path)
        genome_conv = partial(cv.rel_genome_validator, cfg_path=pipeline.cfg_path)
        cfg_fmt = [
            ('genome_fasta', cv.Annot(str, converter=genome_conv)),
            ('output_prefix', cv.Annot(str)),
            ('kmer_k', cv.Annot(int, default=3, converter=cv.nonneg_integer)),
            ('sort_key', cv.Annot(str, default='occ', converter=sort_key_validator)),
            ('gff_exclude_path', cv.Annot(str, default='', converter=optpath_conv,
                                          warn_if_missing=False)),
            ('use_quantiles', cv.Annot(bool, default=True)),
            ('negative_set_gff', cv.Annot(str, converter=relpath_conv)),
            ('n_negative_seqs', cv.Annot(int, default=20000, converter=cv.nonneg_integer)),
            ('remove_tmp_files', cv.Annot(bool, default=True)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']

        table_file = pipeline.get_curfile(fmt='table')
        lo_negset_dir = os.path.join(output_dir, 'lo_negset')
        width = 15
        self._tmp_files.append(lo_negset_dir)
        negset_cmd = [
            'mb-generate-negative-set',
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
            'mb-plot-kmer-logodds',
            '%r' % table_file,
            '%r' % output_dir,
            '%r' % cfg['output_prefix'],
            '%r' % cfg['genome_fasta'],
            '%r' % negset_file,
            '--kmer %s' % cfg['kmer_k'],
            '--key %s' % cfg['sort_key'],
        ]
        if cfg['gff_exclude_path']:
            cmd.append('--filterGFF %r' % cfg['gff_exclude_path'])
        if cfg['use_quantiles']:
            cmd.append('-q')
        if not cfg['remove_tmp_files'] or cfg['keep_all']:
            cmd.append('--keep-tmp-files')
        self._cmds.append(cmd)


class XXmotifModule(pl.CmdPipelineModule):
    def __init__(self, pipeline):
        optpath_conv = partial(opt_file_validator, cfg_path=pipeline.cfg_path)
        relpath_conv = partial(cv.rel_file_r_validator, cfg_path=pipeline.cfg_path)
        genome_conv = partial(cv.rel_genome_validator, cfg_path=pipeline.cfg_path)
        cfg_fmt = [
            ('genome_fasta', cv.Annot(str, converter=genome_conv)),
            ('output_prefix', cv.Annot(str)),
            ('negative_set_gff', cv.Annot(str, converter=relpath_conv)),
            ('n_negative_seqs', cv.Annot(int, default=20000, converter=cv.nonneg_integer)),
            ('plot_top_n_pwm', cv.Annot(int, default=3, converter=cv.nonneg_integer)),
            ('first_index', cv.Annot(int, default=0, converter=cv.nonneg_integer)),
            ('last_index', cv.Annot(int, default=1500, converter=cv.nonneg_integer)),
            ('width', cv.Annot(int, default=12, converter=cv.nonneg_integer)),
            ('sort_key', cv.Annot(str, default='occupancy', converter=sort_key_validator)),
            ('gff_exclude_path', cv.Annot(str, default='', converter=optpath_conv,
                                          warn_if_missing=False)),
            ('gff_padding', cv.Annot(int, default=20, converter=cv.nonneg_integer)),
            ('remove_tmp_files', cv.Annot(bool, default=True)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']
        table_file = pipeline.get_curfile(fmt='table')

        xx_negset_dir = os.path.join(output_dir, 'xx_negset')
        self._tmp_files.append(xx_negset_dir)
        negset_cmd = [
            'mb-generate-negative-set',
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
            'mb-xxmotif',
            table_file,
            '%r' % cfg['genome_fasta'],
            '%r' % output_dir,
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


class HeatmapPlotModule(pl.CmdPipelineModule):
    def __init__(self, pipeline):
        relpath_conv = partial(cv.rel_file_r_validator, cfg_path=pipeline.cfg_path)
        cfg_fmt = [
            ('gff_file', cv.Annot(str, converter=relpath_conv)),
            ('output_prefix', cv.Annot(str)),
            ('downstream_bp', cv.Annot(int, default=4000, converter=cv.nonneg_integer)),
            ('upstream_bp', cv.Annot(int, default=1000, converter=cv.nonneg_integer)),
            ('min_trscr_size_bp', cv.Annot(int, default=0, converter=cv.nonneg_integer)),
            ('max_trscr_size_bp', cv.Annot(int, default=5000, converter=cv.nonneg_integer)),
            ('xbins', cv.Annot(int, default=500, converter=cv.nonneg_integer)),
            ('ybins', cv.Annot(int, default=500, converter=cv.nonneg_integer)),
            ('x_pixels', cv.Annot(int, default=500, converter=cv.nonneg_integer)),
            ('y_pixels', cv.Annot(int, default=500, converter=cv.nonneg_integer)),
            ('remove_tmp_files', cv.Annot(bool, default=True)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']
        table_file = pipeline.get_curfile(fmt='table')

        cmd = [
            'mb-plot-heatmap',
            '%r' % table_file,
            '%r' % output_dir,
            '%r' % cfg['output_prefix'],
            '%r' % cfg['gff_file'],
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


class HeatmapSmallPlotModule(pl.CmdPipelineModule):
    def __init__(self, pipeline):
        relpath_conv = partial(cv.rel_file_r_validator, cfg_path=pipeline.cfg_path)
        cfg_fmt = [
            ('gff_file', cv.Annot(str, converter=relpath_conv)),
            ('output_prefix', cv.Annot(str)),
            ('downstream_bp', cv.Annot(int, default=500, converter=cv.nonneg_integer)),
            ('upstream_bp', cv.Annot(int, default=1000, converter=cv.nonneg_integer)),
            ('min_trscr_size_bp', cv.Annot(int, default=0, converter=cv.nonneg_integer)),
            ('max_trscr_size_bp', cv.Annot(int, default=5000, converter=cv.nonneg_integer)),
            ('xbins', cv.Annot(int, default=500, converter=cv.nonneg_integer)),
            ('ybins', cv.Annot(int, default=500, converter=cv.nonneg_integer)),
            ('x_pixels', cv.Annot(int, default=500, converter=cv.nonneg_integer)),
            ('y_pixels', cv.Annot(int, default=500, converter=cv.nonneg_integer)),
            ('remove_tmp_files', cv.Annot(bool, default=True)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']
        table_file = pipeline.get_curfile(fmt='table')

        cmd = [
            'mb-plot-heatmap-small',
            '%r' % table_file,
            '%r' % output_dir,
            '%r' % cfg['output_prefix'],
            '%r' % cfg['gff_file'],
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


class GffFilterModule(pl.CmdPipelineModule):
    def __init__(self, pipeline):
        relpath_conv = partial(cv.rel_file_r_validator, cfg_path=pipeline.cfg_path)
        cfg_fmt = [
            ('file_postfix', cv.Annot(str, default='fil')),
            ('padding_bp', cv.Annot(int, default=10, converter=cv.nonneg_integer)),
            ('features', cv.Annot(list, default=[])),
            ('filter_gff', cv.Annot(str, converter=relpath_conv)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']
        table_file = pipeline.get_curfile(fmt='table')

        file_name = os.path.basename(table_file)
        bname, ext = os.path.splitext(file_name)
        out_bname = '%s_%s' % (bname, cfg['file_postfix'])
        out_file = os.path.join(output_dir, out_bname + ext)

        cmd = [
            'mb-filter-sites',
            '%r' % table_file,
            '%r' % out_file,
            '%r' % cfg['filter_gff'],
            '--padding_bp %s' % cfg['padding_bp'],
        ]
        if len(cfg['features']) > 0:
            cmd.append('--filter_features')
            for feature in cfg['features']:
                cmd.append('%r' % feature)
        self._cmds.append(cmd)

        self._intermed_files.append(out_file)
        pipeline.upd_curfile(fmt='table', filepath=out_file)
