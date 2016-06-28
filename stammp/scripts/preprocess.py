"""
wrapper to convert raw fastq files from sequencing files to mpileup files. A
fastq-file is adapterclipped, qualityfiltered, mapped and converted.

Command line:
-------------
**Usage:** stammp-preprocess [-h] inputfile outputdir prefix configfile

**Positional arguments:**
  ==========   =================================================
  inputfile    Input fastq file
  outputdir    Output directory
  prefix       Set prefix for all outputfiles
  configfile   Configuration file containing additional options
  ==========   =================================================

**Optional arguments:**
  -h, --help  show this help message and exit

Module usage:
-------------
"""
import argparse
import os
import sys
import configparser
import glob
import logging
from functools import partial
from collections import OrderedDict

from stammp import LOG_DEFAULT_FORMAT, LOG_LEVEL_MAP
from stammp.utils import prepare_output_dir
from stammp.utils import config_validation as cv
from stammp.utils import pipeline as pl

logger = logging.getLogger()
cur_dir = os.path.dirname(os.path.realpath(__file__))
scriptPath = os.path.join(cur_dir, 'utils')


def prepare_dir_or_die(dir_path):
    try:
        prepare_output_dir(dir_path)
    except ValueError as e:
        logger.error('Error while creating output directory: %s', dir_path)
        logger.error(e)
        sys.exit(1)


def main(inputfile, outputdir, prefix, configfile):
    """
    Wrapper method to run 3rd-party programs necessary to pre-process and convert raw fastq files into pileup files for subsequent analysis steps. Detailed parameter settings for the used programs can be modified in the configuration file (:ref:`ref_pre-processing_config`).

    Args:
        inputfile (str): Input fastq file
        outputdir (str): Path to your output directory
        prefix (str): Prefix which is added to all output files e.g. a unique experiment name or protein name
        configfile (str): configuration file where additional options can be set

    Can be accessed directly::

        $ stammp-preprocess.py /path/to/input.fastq /path/output/ prefix /path/to/configfile


    """

    config = configparser.ConfigParser(
        inline_comment_prefixes=';',
        interpolation=configparser.ExtendedInterpolation()
    )
    try:
        config.read(configfile)
    except configparser.Error as e:
        logger.error(e)
        sys.exit(1)

    def mapindex_validator(genome_index):
        genome_index_glob = "%s*" % genome_index
        if len(glob.glob(genome_index_glob)) == 0:
            raise ValueError('genome index %r does not exist' % genome_index)
        return genome_index

    def genomefasta_validator(genome_fasta):
        if not os.path.isfile(genome_fasta):
            raise ValueError('genome fasta file %r does not exist' % genome_fasta)
        return genome_fasta

    mapper_validator = partial(cv.in_set_validator, item_set=set(['STAR']))

    def rel_file_r_validator(path):
        if not os.path.isabs(path):
            parent_path = os.path.dirname(configfile)
            path = os.path.join(os.path.abspath(parent_path), path)
        return cv.file_r_validator(path)

    cfg_format = {
        'general': OrderedDict([
            ('adapter5prime', cv.Annot(str, None, cv.dnastr_validator)),
            ('adapter3prime', cv.Annot(str, None, cv.dnastr_validator)),
            ('genomeindex', cv.Annot(str, None, mapindex_validator)),
            ('genomefasta', cv.Annot(str, None, genomefasta_validator)),
            ('normalization_pileup', cv.Annot(str, None, rel_file_r_validator)),
            ('rmTemp', cv.Annot(bool, True, cv.id_converter)),
            ('n_threads', cv.Annot(int, 2, cv.id_converter)),
        ]),
        'reads': OrderedDict([
            ('fx_Q33', cv.Annot(bool, True, cv.id_converter)),
            ('bc_5prime', cv.Annot(int, 0, cv.nonneg_integer)),
            ('bc_3prime', cv.Annot(int, 0, cv.nonneg_integer)),
            ('min_len', cv.Annot(int, 20, cv.nonneg_integer)),
            ('reference_nucleotide', cv.Annot(str, 'T', cv.dnanuc_validator)),
            ('mutation_nucleotide', cv.Annot(str, 'C', cv.dnanuc_validator)),
        ]),
        'pipeline': OrderedDict([
            ('remove_duplicates', cv.Annot(bool, True, cv.id_converter)),
            ('fastqc_statistics', cv.Annot(bool, True, cv.id_converter)),
            ('quality_trimming', cv.Annot(bool, False, cv.id_converter)),
            ('adapter_clipping', cv.Annot(bool, True, cv.id_converter)),
            ('quality_filtering', cv.Annot(bool, True, cv.id_converter)),
            ('mapping', cv.Annot(str, 'STAR', mapper_validator)),
        ]),
        'fastQC': OrderedDict([
            ('kmer_length', cv.Annot(int, 7, cv.nonneg_integer)),
            ('extra_flags', cv.Annot(str, [], cv.comma_sep_args)),
        ]),
        'clippyAdapterClipper': OrderedDict([
            ('clip_len', cv.Annot(int, 10, cv.nonneg_integer)),
        ]),
        'fastxQualityTrimmer': OrderedDict([
            ('quality_cutoff', cv.Annot(int, 30, cv.nonneg_integer)),
        ]),
        'lafugaQualityFilter': OrderedDict([
            ('quality_cutoff', cv.Annot(int, 30, cv.nonneg_integer)),
            ('chastity', cv.Annot(bool, False, cv.id_converter)),
            ('remove_n', cv.Annot(bool, True, cv.id_converter)),
        ]),
        'STAR': OrderedDict([
            ('n_mismatch', cv.Annot(int, 1, cv.nonneg_integer)),
            ('n_multimap', cv.Annot(int, 1, cv.nonneg_integer)),
            ('extra_flags', cv.Annot(str, [], cv.comma_sep_args)),
            ('allow_soft_clipping', cv.Annot(bool, True, cv.id_converter)),
        ]),
        'PostProcessing': OrderedDict([
            ('plot_transition_profiles', cv.Annot(bool, True, cv.id_converter)),
            ('remove_n_edge_mut', cv.Annot(int, 0, cv.nonneg_integer)),
            ('max_mut_per_read', cv.Annot(int, 1, cv.nonneg_integer)),
            ('min_base_quality', cv.Annot(int, 0, cv.nonneg_integer)),
            ('min_mismatch_quality', cv.Annot(int, 20, cv.nonneg_integer)),
            ('dump_raw_data', cv.Annot(bool, False, cv.id_converter)),
        ]),
        'bsfinder': OrderedDict([
            ('pval_threshold', cv.Annot(float, 0.005, cv.id_converter)),
            ('min_cov', cv.Annot(int, 2, cv.nonneg_integer)),
        ]),
        'normalizer': OrderedDict([
            ('mut_snp_ratio', cv.Annot(float, 0.75, cv.id_converter)),
        ]),
        'max_quantile': OrderedDict([
            ('max_quantile', cv.Annot(float, 0.95, cv.id_converter)),
        ]),
    }

    try:
        cfg_dict = cv.mand_config(config, cfg_format)
    except cv.ConfigError:
        sys.exit(1)

    pipeline = pl.Pipeline(inputfile)
    pipeline_cfg = cfg_dict['pipeline']

    # star pipeline
    if pipeline_cfg['mapping'] == 'STAR':

        # fastqc analysis
        if pipeline_cfg['fastqc_statistics']:
            fastqc_raw_dir = os.path.join(outputdir, 'fastQC_raw')
            prepare_dir_or_die(fastqc_raw_dir)
            fastqc_mod = FastQCModule()
            fastqc_mod.prepare(pipeline.cur_output, fastqc_raw_dir, prefix, cfg_dict)
            fastqc_mod.msg = 'started FastQC analysis of raw reads'
            pipeline.schedule(fastqc_mod)

        # duplicate removal
        if pipeline_cfg['remove_duplicates']:
            dupfil_mod = DuplicateRemovalModule()
            dupfil_mod.prepare(pipeline.cur_output, outputdir, prefix, cfg_dict)
            dupfil_mod.msg = 'started removing duplicated reads'
            pipeline.schedule(dupfil_mod)

        # adapter clipping
        if pipeline_cfg['adapter_clipping']:
            clip_mod = ClippyAdapterClippingModule()
            clip_mod.prepare(pipeline.cur_output, outputdir, prefix, cfg_dict)
            clip_mod.msg = 'started clipping adapter sequences'
            pipeline.schedule(clip_mod)

        # quality trimming
        if pipeline_cfg['quality_trimming']:
            qt_mod = FastxQualityTrimmingModule()
            qt_mod.prepare(pipeline.cur_output, outputdir, prefix, cfg_dict)
            qt_mod.msg = 'started quality trimming'
            pipeline.schedule(qt_mod)

        # quality filtering
        if pipeline_cfg['quality_filtering']:
            qf_mod = LafugaQualityFilterModule()
            qf_mod.prepare(pipeline.cur_output, outputdir, prefix, cfg_dict)
            qf_mod.msg = 'started quality filtering'
            pipeline.schedule(qf_mod)

        # fastqc analysis
        if pipeline_cfg['fastqc_statistics']:
            fastqc_fil_dir = os.path.join(outputdir, 'fastQC_filtered')
            prepare_dir_or_die(fastqc_fil_dir)
            fastqc_mod = FastQCModule()
            fastqc_mod.prepare(pipeline.cur_output, fastqc_fil_dir, prefix, cfg_dict)
            fastqc_mod.msg = 'started FastQC analysis of filtered reads'
            pipeline.schedule(fastqc_mod)

        # mapping with STAR
        star_mod = STARMapModule()
        star_mod.prepare(pipeline.cur_output, outputdir, prefix, cfg_dict)
        star_mod.msg = 'started mapping with STAR'
        pipeline.schedule(star_mod)

        clip_analysis_mod = SoftclipAnalysisModule()
        clip_analysis_mod.prepare(pipeline.cur_output, outputdir, prefix, cfg_dict)
        clip_analysis_mod.msg = 'extracting common soft-clipped sequences'
        pipeline.schedule(clip_analysis_mod)

        # bam postprocessing
        bampp_mod = BamPPModule()
        bampp_mod.prepare(pipeline.cur_output, outputdir, prefix, cfg_dict)
        bampp_mod.msg = 'started postprocessing of the bam file'
        pipeline.schedule(bampp_mod)

        # sort bam file
        sort_mod = SortIndexModule(remove_files=False)
        sort_mod.prepare(pipeline.cur_output, outputdir, prefix)
        sort_mod.msg = 'started sorting the bam file'
        pipeline.schedule(sort_mod)

        # generating the pileup
        pileup_mod = PileupModule()
        pileup_mod.prepare(pipeline.cur_output, outputdir, prefix, cfg_dict)
        pileup_mod.msg = 'started generating the pileup file'
        pipeline.schedule(pileup_mod)

        # finding potential binding sites
        finder = BSFinderModule()
        finder.prepare(pipeline.cur_output, outputdir, prefix, cfg_dict)
        finder.msg = 'scanning for potential binding sites'
        pipeline.schedule(finder)

        # calculate occupancies
        norm_mod = NormalizationModule(remove_files=False)
        norm_mod.prepare(pipeline.cur_output, outputdir, prefix, cfg_dict)
        norm_mod.msg = 'started normalizing with RNAseq data'
        pipeline.schedule(norm_mod)

        # set maximum to quantile
        quantile_mod = MaxQuantileModule(remove_files=False)
        quantile_mod.prepare(pipeline.cur_output, outputdir, prefix, cfg_dict)
        quantile_mod.msg = 'setting maximum quantile'
        pipeline.schedule(quantile_mod)

    for job in pipeline:
        if hasattr(job, 'msg'):
            logger.info(job.msg)
        res = job.execute()
        if res:
            for stdout, stderr in res:
                if stdout.strip() != '':
                    msg = 'Additional Output:\n\n'
                    logger.info(msg + stdout)

    if cfg_dict['general']['rmTemp']:
        logger.info('Started cleaning up temporary files')
        for job in pipeline:
            job.cleanup()

    logger.info('all done. See you soon!')


def run():
    parser = argparse.ArgumentParser(
        description=(
            'Wrapper to convert raw fastq files from sequencing files to mpileup files. '
            'A fastq-file is adapterclipped, qualityfiltered, mapped and converted.'
        ),
    )
    parser.add_argument('inputfile', help='Input fastq file')
    parser.add_argument('outputdir', help='Output directory')
    parser.add_argument('prefix', help='Set prefix for all outputfile')
    parser.add_argument('configfile', help='Config file containing arguments')
    parser.add_argument('--log_level', choices=LOG_LEVEL_MAP.keys(), default='info')
    args = parser.parse_args()
    prepare_dir_or_die(args.outputdir)

    # activate logging
    logging_file = os.path.join(args.outputdir, 'preprocess.log')
    logger = logging.getLogger()
    logger.setLevel(LOG_LEVEL_MAP[args.log_level])
    formatter = logging.Formatter(LOG_DEFAULT_FORMAT)

    console_handler = logging.StreamHandler(stream=sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    file_handler = logging.FileHandler(logging_file, mode='w')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    logger.info('started preprocessing via %r', ' '.join(sys.argv))
    if not os.path.isfile(args.inputfile):
        logger.error('Input fastq file: %r does not exist', args.inputfile)
        sys.exit(1)
    if not os.path.isfile(args.configfile):
        logger.error('Config file: %r does not exist', args.configfile)
        sys.exit(1)
    main(args.inputfile, args.outputdir, args.prefix, args.configfile)


class STARMapModule(pl.CmdPipelineModule):

    def prepare(self, fastq_file, output_dir, prefix, params):
        star_cfg = params['STAR']
        general_cfg = params['general']
        star_output_dir = os.path.join(output_dir, 'STAR_output')
        if not os.path.exists(star_output_dir):
            os.makedirs(star_output_dir)
        star_output_prefix = os.path.join(star_output_dir, prefix + '_')
        sorted_bam_file = star_output_prefix + 'Aligned.sortedByCoord.out.bam'
        self._data['files'].append(sorted_bam_file)
        cmd = [
            'STAR',
            '--readFilesIn %s' % fastq_file,
            '--genomeDir %s' % general_cfg['genomeindex'],
            '--outFilterMultimapNmax %s' % star_cfg['n_multimap'],
            '--outFilterMismatchNmax %s' % star_cfg['n_mismatch'],
            '--outSAMtype BAM SortedByCoordinate',
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
        ]
        if not star_cfg['allow_soft_clipping']:
            cmd.append('--alignEndsType EndToEnd')

        for extra_flag in star_cfg['extra_flags']:
            cmd.append(extra_flag)

        self._cmds.append(cmd)
        self._data['output'] = sorted_bam_file

        index_cmd = [
            'samtools',
            'index',
            sorted_bam_file
        ]
        self._cmds.append(index_cmd)


class FastQCModule(pl.CmdPipelineModule):

    def prepare(self, fastq_file, output_dir, prefix, params):
        general_cfg = params['general']
        fastqc_cfg = params['fastQC']

        adapter_file = os.path.join(output_dir, 'adapters.tmp')
        self._data['files'].append(adapter_file)
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

    def prepare(self, fastq_file, output_dir, prefix, params):
        read_cfg = params['reads']
        qual_stats_file = os.path.join(output_dir, 'fastxstats_raw.txt')
        self._data['files'].append(qual_stats_file)

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
            '-o %s' % os.path.join(output_dir, prefix + '_raw_quality.png'),
            '-t %s' % prefix,  # title
        ]
        cmds.append(qual_bp_toks)

        nucdistr_toks = [
            'fastx_nucleotide_distribution_graph.sh',
            '-i %s' % qual_stats_file,
            '-o %s' % os.path.join(output_dir, prefix + '_raw_nuc.png'),
            '-t %s' % prefix,  # title
        ]
        cmds.append(nucdistr_toks)
        self._cmds.extend(cmds)


class DuplicateRemovalModule(pl.CmdPipelineModule):

    def prepare(self, fastq_file, output_dir, prefix, params):
        rmdup_file = os.path.join(output_dir, prefix + '_nodup.fastq')
        self._data['files'].append(rmdup_file)
        self._data['files'].append(rmdup_file + '.hist')

        cmd = [
            'stammp-removePCRduplicates',
            fastq_file,
            rmdup_file,
        ]
        self._data['output'] = rmdup_file
        self._cmds.append(cmd)


class ClippyAdapterClippingModule(pl.CmdPipelineModule):

    def prepare(self, fastq_file, output_dir, prefix, params):
        general_cfg = params['general']
        read_cfg = params['reads']
        clippy_cfg = params['clippyAdapterClipper']

        adapter_clipped_file = os.path.join(output_dir, prefix + '_adapter.clipped')
        self._data['files'].append(adapter_clipped_file)

        cmd = [
            'stammp-adapter_clipper',
            fastq_file,
            adapter_clipped_file,
            general_cfg['adapter5prime'],
            general_cfg['adapter3prime'],
            '--clip_len %s' % clippy_cfg['clip_len'],
            '--min_len %s' % read_cfg['min_len'],
            '--nt_barcode_5prime %s' % read_cfg['bc_5prime'],
            '--nt_barcode_3prime %s' % read_cfg['bc_3prime'],
            '--verbose'
        ]
        self._cmds.append(cmd)
        self._data['output'] = adapter_clipped_file


class LafugaAdapterClippingModule(pl.CmdPipelineModule):

    def prepare(self, fastq_file, output_dir, prefix, params):
        general_cfg = params['general']
        read_cfg = params['reads']
        clipper_cfg = params['lafugaAdapterClipper']

        bc_5prime = read_cfg['bc_5prime']
        min_len = read_cfg['min_len']
        seed_len = clipper_cfg['seed']
        adapter5prime = general_cfg['adapter5prime']
        adapter3prime = general_cfg['adapter3prime']

        adapter_clipped_file = os.path.join(output_dir, prefix + '_adapter.clipped')
        self._data['files'] = adapter_clipped_file
        clipper_cmd = [
            'java -Xmx4g',
            '-jar %s' % os.path.join(scriptPath, 'AdaptorClipper.jar'),
            '-i %s' % fastq_file,
            '-o %s' % adapter_clipped_file,
            '-a %s,%s' % (adapter5prime, adapter3prime),
            '-m %s' % (min_len + bc_5prime),
            '-seed %s' % seed_len,
        ]
        self._cmds.append(clipper_cmd)
        self._data['output'] = adapter_clipped_file

        if bc_5prime > 0:
            adapter_clipped_bc_file = os.path.join(output_dir, prefix + '_adapter_bc.clipped')
            self._data['files'] = adapter_clipped_bc_file
            bc_trim_cmd = [
                'fastx_trimmer',
                '-i %s' % adapter_clipped_bc_file,
                '-o %s' % adapter_clipped_bc_file,
                '-f %s' % (bc_5prime + 1),  # first base to keep
            ]
            self._cmds.append(bc_trim_cmd)
            self._data['output'] = adapter_clipped_bc_file


class LafugaQualityTrimmingModule(pl.CmdPipelineModule):

    def prepare(self, fastq_file, output_dir, prefix, params):
        read_cfg = params['reads']
        qtrim_cfg = params['lafugaQualityTrimmer']

        qualtrim_file = os.path.join(output_dir, prefix + '_lafqtrim.fastq')
        self._data['files'].append(qualtrim_file)

        qualtrim_toks = [
            'java -Xmx4g',
            '-jar %s' % os.path.join(scriptPath, 'FastqQualityFilter.jar'),
            '-i %s' % fastq_file,
            '-o %s' % qualtrim_file,
            '-m %s' % read_cfg['min_len'],
            '-q %s' % qtrim_cfg['quality_cutoff'],
            '-3',
            '-5',
        ]
        self._cmds = [qualtrim_toks]
        self._data['output'] = qualtrim_file


class FastxQualityTrimmingModule(pl.CmdPipelineModule):

    def prepare(self, fastq_file, output_dir, prefix, params):
        qt_cfg = params['fastxQualityTrimmer']
        read_cfg = params['reads']
        qualtrim_file = os.path.join(output_dir, prefix + '_fstxqtrim.fastq')
        self._data['files'].append(qualtrim_file)

        qt_cmd = [
            'fastq_quality_trimmer',
            '-i %s' % fastq_file,
            '-o %s' % qualtrim_file,
            '-t %s' % qt_cfg['quality_cutoff'],
            '-l %s' % read_cfg['min_len'],
        ]
        if read_cfg['fx_Q33']:
            qt_cmd.append('-Q33')

        self._cmds.append(qt_cmd)
        self._data['output'] = qualtrim_file


class LafugaQualityFilterModule(pl.CmdPipelineModule):

    def prepare(self, fastq_file, output_dir, prefix, params):
        read_cfg = params['reads']
        qf_cfg = params['lafugaQualityFilter']
        qualfil_lafuga_file = os.path.join(output_dir, prefix + '_qfil_lafuga.fastq')
        self._data['files'].append(qualfil_lafuga_file)

        qualfil_toks = [
            'java -Xmx4g',
            '-jar %s' % os.path.join(scriptPath, 'FastqQualityFilter.jar'),
            '-i %s' % fastq_file,
            '-o %s' % qualfil_lafuga_file,
            '-m %s' % read_cfg['min_len'],
            '-q %s' % qf_cfg['quality_cutoff'],
        ]
        if qf_cfg['chastity']:
            qualfil_toks.append('-c Y')
        if qf_cfg['remove_n']:
            qualfil_toks.append('-n')

        self._cmds.append(qualfil_toks)
        self._data['output'] = qualfil_lafuga_file


class SortIndexModule(pl.CmdPipelineModule):

    def prepare(self, bam_file, output_dir, prefix):
        bam_sorted = os.path.join(output_dir, prefix + '_sorted.bam')
        self._data['files'].append(bam_sorted)
        self._data['files'].append(bam_sorted + '.bai')
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
        self._data['output'] = bam_sorted


class PileupModule(pl.CmdPipelineModule):

    def prepare(self, bam_file, output_dir, prefix, cfg_dict):
        general_cfg = cfg_dict['general']
        pileup_file = os.path.join(output_dir, prefix + '.mpileup')
        self._data['output'] = pileup_file
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


class BamPPModule(pl.CmdPipelineModule):

    def prepare(self, bam_file, output_dir, prefix, cfg_dict):
        read_cfg = cfg_dict['reads']
        pp_cfg = cfg_dict['PostProcessing']
        out_bam_file = os.path.join(output_dir, prefix + '_pp.bam')
        pp_plot_dir = os.path.join(output_dir, 'bam_analysis')
        if not os.path.exists(pp_plot_dir):
            os.makedirs(pp_plot_dir)
        self._data['files'].append(out_bam_file)

        transition = read_cfg['reference_nucleotide'] + read_cfg['mutation_nucleotide']
        cmd = [
            'stammp-bam_postprocess',
            bam_file,
            out_bam_file,
            pp_plot_dir,
            '--min-length %s' % read_cfg['min_len'],
            '--mut_edge_bp %s' % pp_cfg['remove_n_edge_mut'],
            '--max_transitions %s' % pp_cfg['max_mut_per_read'],
            '--min_base_quality %s' % pp_cfg['min_base_quality'],
            '--min_mismatch_quality %s' % pp_cfg['min_mismatch_quality'],
            '--transition_of_interest %s' % transition,
        ]
        if pp_cfg['dump_raw_data']:
            cmd.append('--dump_raw_data')

        self._cmds.append(cmd)
        self._data['output'] = out_bam_file


class SoftclipAnalysisModule(pl.CmdPipelineModule):

    def prepare(self, bam_file, output_dir, prefix, cfg_dict):
        pp_plot_dir = os.path.join(output_dir, 'bam_analysis')
        if not os.path.exists(pp_plot_dir):
            os.makedirs(pp_plot_dir)

        cmd = [
            'stammp-softclip_analyzer',
            bam_file,
            pp_plot_dir,
        ]
        self._cmds.append(cmd)


class BSFinderModule(pl.CmdPipelineModule):

    def prepare(self, pileup_file, outdir, prefix, cfg_dict):
        cfg = cfg_dict['bsfinder']
        read_cfg = cfg_dict['reads']
        table_file = os.path.join(outdir, prefix + '.table')
        self._data['output'] = table_file
        self._data['files'].append(table_file)
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


class NormalizationModule(pl.CmdPipelineModule):

    def prepare(self, table_file, outdir, prefix, cfg_dict):
        cfg = cfg_dict['normalizer']
        normed_table_file = os.path.join(outdir, prefix + '.normed_table')
        self._data['output'] = normed_table_file
        self._data['files'].append(normed_table_file)
        cmd = [
            'stammp-normalize',
            table_file,
            normed_table_file,
            cfg_dict['general']['normalization_pileup'],
            '--mut_snp_ratio %s' % cfg['mut_snp_ratio'],
        ]
        self._cmds.append(cmd)


class MaxQuantileModule(pl.CmdPipelineModule):

    def prepare(self, norm_table_file, outdir, prefix, cfg_dict):
        cfg = cfg_dict['max_quantile']
        maxq_file = os.path.join(outdir, prefix + '.qtable')
        self._data['output'] = maxq_file
        self._data['files'].append(maxq_file)
        cmd = [
            'stammp-convert2quantile',
            norm_table_file,
            maxq_file,
            '-q %s' % cfg['max_quantile'],
        ]
        self._cmds.append(cmd)


if __name__ == '__main__':
    run()
