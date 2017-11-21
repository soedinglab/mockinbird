import os
from mockinbird.utils import pipeline as pl
from mockinbird.utils import config_validation as cv


class NaiveBSFinderModule(pl.CmdPipelineModule):

    def __init__(self, pipeline):
        cfg_fmt = [
            ('min_transitions', cv.Annot(cv.nonneg_integer, default=2)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        read_cfg = pipeline.get_config('reads')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']
        pileup_file = pipeline.get_curfile(fmt='mpileup')

        table_file = os.path.join(output_dir, prefix + '.pre_table')
        cmd = [
            'mb-naive-bsfinder',
            '-r %s' % read_cfg['reference_nucleotide'],
            '-m %s' % read_cfg['mutation_nucleotide'],
            pileup_file,
            table_file
        ]
        self._cmds.append(cmd)
        self._intermed_files.append(table_file)
        pipeline.upd_curfile(fmt='table', filepath=table_file)


class BSFinderModule(pl.CmdPipelineModule):

    def __init__(self, pipeline):
        cfg_fmt = [
            ('pval_threshold', cv.Annot(float, default=0.005)),
            ('min_cov', cv.Annot(int, default=2, converter=cv.nonneg_integer)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        read_cfg = pipeline.get_config('reads')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']
        pileup_file = pipeline.get_curfile(fmt='mpileup')

        table_file = os.path.join(output_dir, prefix + '.pre_table')
        cmd = [
            'mb-bsfinder',
            '-p %s' % cfg['pval_threshold'],
            '-c %s' % cfg['min_cov'],
            '-r %s' % read_cfg['reference_nucleotide'],
            '-m %s' % read_cfg['mutation_nucleotide'],
            pileup_file,
            table_file
        ]
        self._cmds.append(cmd)
        self._intermed_files.append(table_file)
        pipeline.upd_curfile(fmt='table', filepath=table_file)


class MockinbirdModule(pl.CmdPipelineModule):

    def __init__(self, pipeline):
        cfg_fmt = [
            ('plot_dir', cv.Annot(str, default='mockinbird_plots')),
            ('max_k_mock', cv.Annot(int, default=10, converter=cv.nonneg_integer)),
            ('min_k', cv.Annot(int, default=2, converter=cv.nonneg_integer)),
            ('min_post', cv.Annot(float, default=0.1)),
            ('extra_args', cv.Annot(list, default=[])),
            ('null_fraction', cv.Annot(float, default=1)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline

        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']

        prefix = general_cfg['prefix']
        trtable_file = pipeline.get_curfile(fmt='trtable')
        mock_model = pipeline.get_curfile(fmt='mock_model')
        post_file = os.path.join(output_dir, prefix + '.post')

        stat_file = pipeline.get_curfile(fmt='stat_file')
        plot_dir = os.path.join(output_dir, cfg['plot_dir'])

        cmd = [
            'mb-calculate-posterior',
            '--plot_dir %r' % plot_dir,
            '--max_k_mock %s' % cfg['max_k_mock'],
            '--bam_statistics_json %r' % stat_file,
            '--null_fraction %s' % cfg['null_fraction'],
            '%r' % trtable_file,
            '%r' % mock_model,
            '%r' % post_file,
        ]
        if cfg['extra_args']:
            cmd.extend(cfg['extra_args'])
        self._cmds.append(cmd)

        pipeline.upd_curfile(fmt='post', filepath=post_file)
        self._intermed_files.append(post_file)

        table_file = os.path.join(output_dir, prefix + '.table')
        cmd = [
            'mb-mockinbird2table',
            '%r' % post_file,
            '%r' % table_file,
            '--post_thresh %s' % cfg['min_post'],
            '--k_thresh %s' % cfg['min_k'],
        ]
        self._cmds.append(cmd)

        pipeline.upd_curfile(fmt='table', filepath=table_file)
        self._intermed_files.append(table_file)
