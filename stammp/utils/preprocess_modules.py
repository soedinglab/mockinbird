import os
import time
from functools import partial

from stammp.utils import pipeline as pl
from stammp.utils import config_validation as cv


class STARMapModule(pl.CmdPipelineModule):

    def __init__(self, pipeline):
        relgen_conv = partial(cv.rel_mapindex_validator, cfg_path=pipeline.cfg_path)
        cfg_fmt = [
            ('genome_index', cv.Annot(str, converter=relgen_conv)),
            ('n_mismatch', cv.Annot(int, default=1, converter=cv.nonneg_integer)),
            ('n_multimap', cv.Annot(int, default=1, converter=cv.nonneg_integer)),
            ('extra_flags', cv.Annot(list, default=[])),
            ('allow_soft_clipping', cv.Annot(cv.boolean, default=True)),
            ('outdir_name', cv.Annot(str, default='star_out')),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']
        star_dirname = cfg['outdir_name']
        fastq_file = pipeline.get_curfile(fmt='fastq')

        star_output_dir = os.path.join(output_dir, star_dirname)
        if not os.path.exists(star_output_dir):
            os.makedirs(star_output_dir)
        star_output_prefix = os.path.join(star_output_dir, prefix + '_')
        bam_file = star_output_prefix + 'Aligned.out.bam'
        unmapped = None
        self._tmp_files.append(bam_file)
        cmd = [
            'STAR',
            '--readFilesIn %s' % fastq_file,
            '--genomeDir %s' % cfg['genome_index'],
            '--outFilterMultimapNmax %s' % cfg['n_multimap'],
            '--outFilterMismatchNmax %s' % cfg['n_mismatch'],
            '--outSAMtype BAM Unsorted',
            '--outFileNamePrefix %s' % star_output_prefix,
            '--runThreadN %s' % general_cfg['n_threads'],
            # we do not want any insertions or deletions
            '--scoreDelOpen -10000',
            '--scoreInsOpen -10000',
            # with all the unmappable reads it's way too risky to call new junctions
            '--alignSJoverhangMin 10000',
            # it is not allowed to use mismatches to build splice sites
            '--alignSJstitchMismatchNmax 0 0 0 0',
            '--outSAMattributes NH HI AS NM MD',
            '--outReadsUnmapped Fastx',
        ]
        if not cfg['allow_soft_clipping']:
            cmd.append('--alignEndsType EndToEnd')
        for extra_flag in cfg['extra_flags']:
            cmd.append(extra_flag)
        self._cmds.append(cmd)

        sort_bam = star_output_prefix + 'sorted.bam'
        sort_cmd = [
            'samtools sort',
            bam_file,
            '-@ %s' % general_cfg['n_threads'],
            '| samtools view -h -b > %r' % sort_bam,
        ]
        self._cmds.append(sort_cmd)

        if pipeline.has_curfile(fmt='bam'):
            self._tmp_files.append(sort_bam)
            merged_bam = os.path.join(output_dir, prefix + '_merged%s.bam' % int(time.time()))
            merge_cmd = [
                'samtools',
                'merge',
                '%r' % merged_bam,
                '%r' % pipeline.get_curfile('bam'),
                '%r' % sort_bam,
            ]
            self._cmds.append(merge_cmd)
            out_bam = merged_bam
        else:
            out_bam = sort_bam

        self._intermed_files.append(out_bam)
        pipeline.upd_curfile(fmt='bam', filepath=out_bam)
        pipeline.upd_curfile(fmt='fastq', filepath=unmapped)

        index_cmd = [
            'samtools',
            'index',
            out_bam
        ]
        self._cmds.append(index_cmd)

        self._intermed_files.append(out_bam + '.bai')
        self._intermed_files.append(out_bam)
        pipeline.upd_curfile(fmt='bam', filepath=out_bam)


class BowtieMapModule(pl.CmdPipelineModule):

    def __init__(self, pipeline):
        relgen_conv = partial(cv.rel_mapindex_validator, cfg_path=pipeline.cfg_path)
        cfg_fmt = [
            ('genome_index', cv.Annot(str, converter=relgen_conv)),
            ('n_mismatch', cv.Annot(int, default=1, converter=cv.nonneg_integer)),
            ('n_mismatch', cv.Annot(int, default=1, converter=cv.nonneg_integer)),
            ('n_multimap', cv.Annot(int, default=1, converter=cv.nonneg_integer)),
            ('extra_flags', cv.Annot(str, default=[], converter=cv.comma_sep_args)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']
        fastq_file = pipeline.get_curfile(fmt='fastq')

        sam_file = os.path.join(output_dir, prefix + '.sam')
        unmapped = os.path.join(output_dir, prefix + '_unmapped.fq')

        cmd = [
            'bowtie',
            '-m %s' % cfg['n_multimap'],
            '-v %s' % cfg['n_mismatch'],
            '%r' % cfg['genome_index'],
            '%r' % fastq_file,
            '%r' % sam_file,
            '-p %s' % general_cfg['n_threads'],
            '-S',
            '--un %r' % unmapped,
            '2>&1',
        ]
        for extra_flag in cfg['extra_flags']:
            cmd.append(extra_flag)
        self._cmds.append(cmd)

        self._tmp_files.append(sam_file)
        bam_file = os.path.join(output_dir, prefix + '_raw_sorted.bam')
        self._tmp_files.append(bam_file)

        sort_cmd = [
            'samtools sort',
            sam_file,
            '-@ %s' % general_cfg['n_threads'],
            '| samtools view -h -b > %r' % bam_file,
        ]
        self._cmds.append(sort_cmd)

        calmd_file = os.path.join(output_dir, prefix + '.bam')
        calmd_cmd = [
            'samtools calmd',
            '%r' % bam_file,
            '%r' % general_cfg['genomefasta'],
            '-b',
            '> %r' % calmd_file,
        ]
        self._cmds.append(calmd_cmd)

        if pipeline.has_curfile(fmt='bam'):
            self._tmp_files.append(calmd_cmd)
            merged_bam = os.path.join(output_dir, prefix + '_merged%s.bam' % int(time.time()))
            merge_cmd = [
                'samtools',
                'merge',
                '%r' % merged_bam,
                '%r' % pipeline.get_curfile('bam'),
                '%r' % calmd_cmd,
            ]
            self._cmds.append(merge_cmd)
            out_bam = merged_bam
        else:
            out_bam = calmd_file

        self._intermed_files.append(out_bam)
        pipeline.upd_curfile(fmt='bam', filepath=out_bam)
        pipeline.upd_curfile(fmt='fastq', filepath=unmapped)

        index_cmd = [
            'samtools',
            'index',
            out_bam
        ]
        self._cmds.append(index_cmd)
        self._intermed_files.append(out_bam + '.bai')


class FastQCModule(pl.CmdPipelineModule):

    def __init__(self, pipeline):
        cfg_fmt = [
            ('kmer_length', cv.Annot(int, default=7, converter=cv.nonneg_integer)),
            ('extra_flags', cv.Annot(list, default=[])),
            ('outdir_name', cv.Annot(str, 'fastQC')),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']
        fastqc_cfg = cfg
        output_dir = os.path.join(output_dir, fastqc_cfg['outdir_name'])
        fastq_file = pipeline.get_curfile(fmt='fastq')

        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        adapter_file = os.path.join(output_dir, 'adapters.tmp')
        self._tmp_files.append(adapter_file)
        with open(adapter_file, 'w') as fc:
            print('adapter5prime\t %s' % general_cfg['adapter5prime'], file=fc)
            print('adapter3prime\t %s' % general_cfg['adapter3prime'], file=fc)

        cmd = [
            'fastqc',
            '-o %s' % output_dir,
            '-f fastq',
            '--threads %s' % general_cfg['n_threads'],
            '--kmers %s' % fastqc_cfg['kmer_length'],
            '--adapters %s' % adapter_file,
            '-d %s' % output_dir,  # temp directory
            fastq_file,
        ]
        for extra_flag in fastqc_cfg['extra_flags']:
            cmd.append(extra_flag)

        self._cmds.append(cmd)


class FastXStatisticsModule(pl.CmdPipelineModule):

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        read_cfg = pipeline.get_config('reads')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']
        fastq_file = pipeline.get_curfile(fmt='fastq')

        qual_stats_file = os.path.join(output_dir, 'fastxstats_raw.txt')
        self._tmp_files.append(qual_stats_file)

        cmds = []
        qual_stats_toks = [
            'fastx_quality_stats',
            '-i %s' % fastq_file,
            '-o %s' % qual_stats_file,
        ]
        if read_cfg['fx_Q33']:
            qual_stats_toks.append('-Q33')
        cmds.append(qual_stats_toks)

        qual_bp_toks = [
            'fastq_quality_boxplot_graph.sh',
            '-i %s' % qual_stats_file,
            '-o %s' % os.path.join(output_dir, prefix + '_raw_quality.pdf'),
            '-t %s' % prefix,  # title
        ]
        cmds.append(qual_bp_toks)

        nucdistr_toks = [
            'fastx_nucleotide_distribution_graph.sh',
            '-i %s' % qual_stats_file,
            '-o %s' % os.path.join(output_dir, prefix + '_raw_nuc.pdf'),
            '-t %s' % prefix,  # title
        ]
        cmds.append(nucdistr_toks)
        self._cmds.extend(cmds)


class DuplicateRemovalModule(pl.CmdPipelineModule):

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']
        fastq_file = pipeline.get_curfile(fmt='fastq')

        rmdup_file = os.path.join(output_dir, prefix + '_nodup.fastq')
        self._tmp_files.append(rmdup_file + '.hist')

        cmd = [
            'stammp-removePCRduplicates',
            fastq_file,
            rmdup_file,
        ]
        self._intermed_files.append(rmdup_file)
        pipeline.upd_curfile(fmt='fastq', filepath=rmdup_file)
        self._cmds.append(cmd)


class UmiToolsExtractModule(pl.CmdPipelineModule):

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        read_cfg = pipeline.get_config('reads')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']
        fastq_file = pipeline.get_curfile(fmt='fastq')
        log_file = os.path.join(output_dir, 'umi_extr.log')

        extr_file = os.path.join(output_dir, prefix + '_extracted.fastq')

        cmd = [
            'umi_tools',
            'extract',
            '-I %r' % fastq_file,
            '--bc-pattern=%s' % ('N' * read_cfg['bc_5prime']),
            '-L %r' % log_file,
            '-E %r' % log_file,
            '-S %r' % extr_file,
        ]
        self._intermed_files.append(extr_file)
        pipeline.upd_curfile(fmt='fastq', filepath=extr_file)
        self._cmds.append(cmd)


class UmiToolsDedupModule(pl.CmdPipelineModule):

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']

        bam_file = pipeline.get_curfile(fmt='bam')
        stat_file = os.path.join(output_dir, 'umi_dedup.stat')
        dedup_bam = os.path.join(output_dir, prefix + '_dedup.bam')

        cmd = [
            'umi_tools',
            'dedup',
            '-I %r' % bam_file,
            '--output-stats=%r' % stat_file,
            '-S %r' % dedup_bam,
        ]
        self._intermed_files.append(dedup_bam)
        pipeline.upd_curfile(fmt='bam', filepath=dedup_bam)
        self._cmds.append(cmd)


class SkewerAdapterClippingModule(pl.CmdPipelineModule):

    def __init__(self, pipeline):
        cfg_fmt = [
            ('extra_args', cv.Annot(list, default=[])),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        read_cfg = pipeline.get_config('reads')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']
        fastq_file = pipeline.get_curfile(fmt='fastq')

        adapter_clipped_file = os.path.join(output_dir, prefix + '_skewer.clipped')

        cmd = [
            'skewer',
            fastq_file,
            '-x %s' % general_cfg['adapter3prime'],
            '-m tail',
            '--min %s' % read_cfg['min_len'],
            '--quiet',
            '--stdout',
            '> %r' % adapter_clipped_file,
        ]

        if cfg['extra_args']:
            cmd.extend(cfg['extra_args'])
        self._cmds.append(cmd)

        self._intermed_files.append(adapter_clipped_file)
        pipeline.upd_curfile(fmt='fastq', filepath=adapter_clipped_file)


class ClippyAdapterClippingModule(pl.CmdPipelineModule):

    def __init__(self, pipeline):
        cfg_fmt = [
            ('clip_len', cv.Annot(int, default=10, converter=cv.nonneg_integer)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        read_cfg = pipeline.get_config('reads')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']
        fastq_file = pipeline.get_curfile(fmt='fastq')

        adapter_clipped_file = os.path.join(output_dir, prefix + '_adapter.clipped')

        cmd = [
            'stammp-adapter-clipper',
            fastq_file,
            adapter_clipped_file,
            general_cfg['adapter5prime'],
            general_cfg['adapter3prime'],
            '--clip_len %s' % cfg['clip_len'],
            '--min_len %s' % read_cfg['min_len'],
            '--nt_barcode_5prime %s' % read_cfg['bc_5prime'],
            '--nt_barcode_3prime %s' % read_cfg['bc_3prime'],
            '--plot_dir %s' % output_dir,
            '--verbose'
        ]
        self._cmds.append(cmd)

        self._intermed_files.append(adapter_clipped_file)
        pipeline.upd_curfile(fmt='fastq', filepath=adapter_clipped_file)


class FastxQualityTrimmingModule(pl.CmdPipelineModule):

    def __init__(self, pipeline):
        cfg_fmt = [
            ('quality_cutoff', cv.Annot(int, default=30, converter=cv.nonneg_integer)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        read_cfg = pipeline.get_config('reads')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']
        fastq_file = pipeline.get_curfile(fmt='fastq')

        qualtrim_file = os.path.join(output_dir, prefix + '_fstxqtrim.fastq')

        qt_cmd = [
            'fastq_quality_trimmer',
            '-i %s' % fastq_file,
            '-o %s' % qualtrim_file,
            '-t %s' % cfg['quality_cutoff'],
            '-l %s' % read_cfg['min_len'],
        ]
        if read_cfg['fx_Q33']:
            qt_cmd.append('-Q33')

        self._cmds.append(qt_cmd)

        self._intermed_files.append(qualtrim_file)
        pipeline.upd_curfile(fmt='fastq', filepath=qualtrim_file)


class SortIndexModule(pl.CmdPipelineModule):

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']
        bam_file = pipeline.get_curfile(fmt='bam')

        bam_sorted = os.path.join(output_dir, prefix + '_sim_sorted.bam')
        sort_cmd = [
            'samtools',
            'sort',
            bam_file,
            '> %s' % bam_sorted,
        ]
        self._cmds.append(sort_cmd)

        index_cmd = [
            'samtools',
            'index',
            bam_sorted,
        ]
        self._cmds.append(index_cmd)
        self._intermed_files.append(bam_sorted)
        self._intermed_files.append(bam_sorted + '.bai')
        pipeline.upd_curfile(fmt='bam', filepath=bam_sorted)


class PileupModule(pl.CmdPipelineModule):

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']
        bam_file = pipeline.get_curfile(fmt='bam')

        pileup_file = os.path.join(output_dir, prefix + '.mpileup')
        pileup_cmd = [
            'samtools',
            'mpileup',
            '-C 0',  # disable adjust mapping quality
            '-d 100000',  # max depth
            '-q 0',  # minimum alignment quality
            '-Q 0',  # minimum base quality
            '-f %s' % general_cfg['genomefasta'],
            bam_file,
            '> %s' % pileup_file,
        ]
        self._cmds.append(pileup_cmd)
        self._intermed_files.append(pileup_file)
        pipeline.upd_curfile(fmt='mpileup', filepath=pileup_file)


class BamPPModule(pl.CmdPipelineModule):

    def __init__(self, pipeline):
        cfg_fmt = [
            ('remove_n_edge_mut', cv.Annot(int, default=0, converter=cv.nonneg_integer)),
            ('max_mut_per_read', cv.Annot(int, default=1, converter=cv.nonneg_integer)),
            ('min_base_quality', cv.Annot(int, default=0, converter=cv.nonneg_integer)),
            ('min_avg_ali_quality', cv.Annot(int, default=20, converter=cv.nonneg_integer)),
            ('min_mismatch_quality', cv.Annot(int, default=20, converter=cv.nonneg_integer)),
            ('dump_raw_data', cv.Annot(cv.boolean, default=False)),
            ('outdir_name', cv.Annot(str, default='bam_analysis')),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        read_cfg = pipeline.get_config('reads')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']
        bam_file = pipeline.get_curfile(fmt='bam')

        out_bam_file = os.path.join(output_dir, prefix + '_pp.bam')
        pp_plot_dir = os.path.join(output_dir, cfg['outdir_name'])
        if not os.path.exists(pp_plot_dir):
            os.makedirs(pp_plot_dir)

        transition = read_cfg['reference_nucleotide'] + read_cfg['mutation_nucleotide']
        cmd = [
            'stammp-bam-postprocess',
            bam_file,
            out_bam_file,
            pp_plot_dir,
            '--min-length %s' % read_cfg['min_len'],
            '--mut_edge_bp %s' % cfg['remove_n_edge_mut'],
            '--max_transitions %s' % cfg['max_mut_per_read'],
            '--min_base_quality %s' % cfg['min_base_quality'],
            '--min_mismatch_quality %s' % cfg['min_mismatch_quality'],
            '--transition_of_interest %s' % transition,
        ]
        if cfg['dump_raw_data']:
            cmd.append('--dump_raw_data')

        self._cmds.append(cmd)
        self._intermed_files.append(out_bam_file)
        pipeline.upd_curfile(fmt='bam', filepath=out_bam_file)


class SoftclipAnalysisModule(pl.CmdPipelineModule):

    def __init__(self, pipeline):
        cfg_fmt = [
            ('outdir_name', cv.Annot(str, 'bam_analysis')),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']
        bam_file = pipeline.get_curfile(fmt='bam')

        pp_plot_dir = os.path.join(output_dir, 'bam_analysis')
        if not os.path.exists(pp_plot_dir):
            os.makedirs(pp_plot_dir)

        cmd = [
            'stammp-softclip-analyzer',
            bam_file,
            pp_plot_dir,
        ]
        self._cmds.append(cmd)


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
            'stammp-bsfinder',
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


class NormalizationModule(pl.CmdPipelineModule):
    def __init__(self, pipeline):
        cfg_fmt = [
            ('mut_snp_ratio', cv.Annot(float, default=0.75)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']
        table_file = pipeline.get_curfile(fmt='table')

        normed_table_file = os.path.join(output_dir, prefix + '.normed_table')
        cmd = [
            'stammp-normalize',
            table_file,
            normed_table_file,
            general_cfg['normalization_pileup'],
            '--mut_snp_ratio %s' % cfg['mut_snp_ratio'],
        ]
        self._cmds.append(cmd)

        self._intermed_files.append(normed_table_file)
        pipeline.upd_curfile(fmt='table', filepath=normed_table_file)


class MaxQuantileModule(pl.CmdPipelineModule):

    def __init__(self, pipeline):
        cfg_fmt = [
            ('max_quantile', cv.Annot(float, default=0.95)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']
        prefix = general_cfg['prefix']
        table_file = pipeline.get_curfile(fmt='table')

        maxq_file = os.path.join(output_dir, prefix + '.table')
        cmd = [
            'stammp-convert2quantile',
            table_file,
            maxq_file,
            '-q %s' % cfg['max_quantile'],
        ]
        self._cmds.append(cmd)
        self._intermed_files.append(maxq_file)
        pipeline.upd_curfile(fmt='table', filepath=maxq_file)
