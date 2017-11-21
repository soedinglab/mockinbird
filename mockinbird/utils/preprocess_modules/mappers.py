import os
import time
from functools import partial
from mockinbird.utils import pipeline as pl
from mockinbird.utils import config_validation as cv


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
        unmapped = star_output_prefix + 'Unmapped.out.mate1'
        self._intermed_files.append(unmapped)

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

        log_file = star_output_prefix + 'Log.final.out'
        print_cmd = ['cat', '%r' % log_file]
        self._cmds.append(print_cmd)

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
        self._intermed_files.append(unmapped)

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
