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
from collections import namedtuple
import re
from operator import attrgetter
from functools import partial

from stammp import LOG_DEFAULT_FORMAT, LOG_LEVEL_MAP
from stammp.utils import prepare_output_dir, execute

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

    config = configparser.ConfigParser(inline_comment_prefixes=';')
    config.read(configfile)

    Annot = namedtuple('Annotation', ['type', 'default', 'converter'])

    def mapindex_validator(genome_index):
        genome_index_glob = "%s*" % genome_index
        if len(glob.glob(genome_index_glob)) == 0:
            raise ValueError('genome index %r does not exist' % genome_index)
        return genome_index

    def genomefasta_validator(genome_fasta):
        if not os.path.isfile(genome_fasta):
            raise ValueError('genome fasta file %r does not exist' % genome_fasta)
        return genome_fasta

    def nucstr_validator(nucleotide_string):
        nuc_pat = re.compile('[ACTGactg]+')
        if not nuc_pat.match(nucleotide_string):
            raise ValueError('adaptor sequence %r is invalid. Adaptor sequences may only '
                             'contain A, C, G or T nucleotides.' % nucleotide_string)
        return nucleotide_string

    def nonneg_integer(integer):
        if integer < 0:
            raise ValueError('Non-negative integer expected. Got %s.' % integer)
        return integer

    def in_set_validator(item, item_set):
        if item not in item_set:
            raise ValueError('%r is not in set %r' % (item, item_set))
        return item

    def is_subset_validator(item_str, item_set):
        items = []
        for item in item_str.split(','):
            in_set_validator(item)
            items.append(item)
        return items

    def comma_sep_args(item_str):
        if item_str.strip() == '':
            return []
        else:
            return item_str.split(',')

    trivial_converter = lambda x: x
    mapper_validator = partial(in_set_validator, item_set=set(['STAR']))

    config_format = {
        'general': {
            'adapter5prime': Annot(str, None, nucstr_validator),
            'adapter3prime': Annot(str, None, nucstr_validator),
            'genomeindex': Annot(str, None, mapindex_validator),
            'genomefasta': Annot(str, None, genomefasta_validator),
            'rmTemp': Annot(bool, True, trivial_converter),
            'n_threads': Annot(int, 2, trivial_converter),
        },
        'reads': {
            'fx_Q33': Annot(bool, True, trivial_converter),
            'bc_5prime': Annot(int, 0, nonneg_integer),
            'bc_3prime': Annot(int, 0, nonneg_integer),
            'min_len': Annot(int, 15, nonneg_integer),
        },
        'pipeline': {
            'remove_duplicates': Annot(bool, True, trivial_converter),
            'fastqc_statistics': Annot(bool, True, trivial_converter),
            'quality_trimming': Annot(bool, False, trivial_converter),
            'adapter_clipping': Annot(bool, True, trivial_converter),
            'quality_filtering': Annot(bool, True, trivial_converter),
            'mapping': Annot(str, 'STAR', mapper_validator),
        },
        'fastQC': {
            'kmer_length': Annot(int, 7, nonneg_integer),
            'extra_flags': Annot(str, [], comma_sep_args),
        },
        'clippyAdapterClipper': {
            'clip_len': Annot(int, 10, nonneg_integer),
        },
        'fastxQualityTrimmer': {
            'quality_cutoff': Annot(int, 30, nonneg_integer),
        },
        'lafugaQualityFilter': {
            'quality_cutoff': Annot(int, 30, nonneg_integer),
            'chastity': Annot(bool, False, trivial_converter),
            'remove_n': Annot(bool, True, trivial_converter),
        },
        'STAR': {
            'n_mismatch': Annot(int, 1, nonneg_integer),
            'n_multimap': Annot(int, 1, nonneg_integer),
            'extra_flags': Annot(str, [], comma_sep_args),
            'allow_soft_clipping': Annot(bool, True, trivial_converter),
        },
        'PostProcessing': {
            'plot_transition_profiles': Annot(bool, True, trivial_converter),
            'remove_n_edge_mut': Annot(int, 0, nonneg_integer),
            'max_mut_per_read': Annot(int, 1, nonneg_integer),
        }
    }

    type_getter = {
        int: attrgetter('getint'),
        str: attrgetter('get'),
        float: attrgetter('getfloat'),
        bool: attrgetter('getboolean'),
    }
    cfg_dict = {}
    for section, fields_dict in config_format.items():
        if section not in config:
            logger.warn('config section %r missing.', section)
            config.add_section(section)
        cfg_sec = config[section]
        sec_dict = {}
        for key, annot in fields_dict.items():
            if key not in cfg_sec:
                if annot.default is not None:
                    logger.warn('key %r missing in section %r. '
                                'Using default value %r.', key, section, annot.default)
                    sec_dict[key] = annot.default
                else:
                    logger.error('mandatory configuration key %r missing in section %r.',
                                 key, section)
                    sys.exit(1)
            else:
                try:
                    cfg_getter = type_getter[annot.type]
                    raw_data = cfg_getter(cfg_sec)(key)
                    sec_dict[key] = annot.converter(raw_data)
                    logger.debug('section %r: set %r to %r', section,
                                 key, raw_data)
                except ValueError as ex:
                    logger.error('invalid value for key %r in section %r', key, section)
                    logger.error(str(ex))
                    sys.exit(1)
        cfg_dict[section] = sec_dict

    pipeline = Pipeline()
    pipeline.schedule(DummyModule(inputfile))
    pipeline_cfg = cfg_dict['pipeline']

    # star pipeline
    if pipeline_cfg['mapping'] == 'STAR':

        # fastqc analysis
        if pipeline_cfg['fastqc_statistics']:
            fastqc_raw_dir = os.path.join(outputdir, 'fastQC_raw')
            prepare_dir_or_die(fastqc_raw_dir)
            fastqc_mod = FastQCModule()
            fastqc_mod.prepare(pipeline.cur_output(), fastqc_raw_dir, prefix, cfg_dict)
            fastqc_mod.msg = 'started FastQC analysis of raw reads'
            pipeline.schedule(fastqc_mod)

        # duplicate removal
        if pipeline_cfg['remove_duplicates']:
            dupfil_mod = DuplicateRemovalModule()
            dupfil_mod.prepare(pipeline.cur_output(), outputdir, prefix, cfg_dict)
            dupfil_mod.msg = 'started removing duplicated reads'
            pipeline.schedule(dupfil_mod)

        # adapter clipping
        if pipeline_cfg['adapter_clipping']:
            clip_mod = ClippyAdapterClippingModule()
            clip_mod.prepare(pipeline.cur_output(), outputdir, prefix, cfg_dict)
            clip_mod.msg = 'started clipping adapter sequences'
            pipeline.schedule(clip_mod)

        # quality trimming
        if pipeline_cfg['quality_trimming']:
            qt_mod = FastxQualityTrimmingModule()
            qt_mod.prepare(pipeline.cur_output(), outputdir, prefix, cfg_dict)
            qt_mod.msg = 'started quality trimming'
            pipeline.schedule(qt_mod)

        # quality filtering
        if pipeline_cfg['quality_filtering']:
            qf_mod = LafugaQualityFilterModule()
            qf_mod.prepare(pipeline.cur_output(), outputdir, prefix, cfg_dict)
            qf_mod.msg = 'started quality filtering'
            pipeline.schedule(qf_mod)

        # fastqc analysis
        if pipeline_cfg['fastqc_statistics']:
            fastqc_fil_dir = os.path.join(outputdir, 'fastQC_filtered')
            prepare_dir_or_die(fastqc_fil_dir)
            fastqc_mod = FastQCModule()
            fastqc_mod.prepare(pipeline.cur_output(), fastqc_fil_dir, prefix, cfg_dict)
            fastqc_mod.msg = 'started FastQC analysis of filtered reads'
            pipeline.schedule(fastqc_mod)

        # mapping with STAR
        star_mod = STARMapModule()
        star_mod.prepare(pipeline.cur_output(), outputdir, prefix, cfg_dict)
        star_mod.msg = 'started mapping with STAR'
        pipeline.schedule(star_mod)

        clip_analysis_mod = SoftclipAnalysisModule()
        clip_analysis_mod.prepare(pipeline.cur_output(), outputdir, prefix, cfg_dict)
        clip_analysis_mod.msg = 'extracting common soft-clipped sequences'
        pipeline.schedule(clip_analysis_mod)

        # bam postprocessing
        bampp_mod = BamPPModule()
        bampp_mod.prepare(pipeline.cur_output(), outputdir, prefix, cfg_dict)
        bampp_mod.msg = 'started postprocessing of the bam file'
        pipeline.schedule(bampp_mod)

        # sort bam file
        sort_mod = SortIndexModule(remove_files=False)
        sort_mod.prepare(pipeline.cur_output(), outputdir, prefix)
        sort_mod.msg = 'started sorting the bam file'
        pipeline.schedule(sort_mod)

        # generating the pileup
        pileup_mod = PileupModule()
        pileup_mod.prepare(pipeline.cur_output(), outputdir, prefix, cfg_dict)
        pileup_mod.msg = 'started generating the pileup file'
        pipeline.schedule(pileup_mod)

    for job in pipeline:
        if hasattr(job, 'msg'):
            logger.info(job.msg)
        res = job.execute()
        if res:
            for stdout, stderr in res:
                if stdout.strip() != '':
                    logger.info('\n\n' + stdout)

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


class PipelineModule:

    def __init__(self):
        self._data = {}

    def prepare(self):
        pass

    def execute(self):
        pass

    def cleanup(self):
        pass

    def get(self, key='output', default=None):
        return self._data.get(key, default)


class CmdPipelineModule(PipelineModule):

    def __init__(self, remove_files=True):
        super().__init__()
        self._remove_files = remove_files
        self._data['files'] = []
        self._cmds = []

    def execute(self):
        for cmd in self._cmds:
            yield execute(cmd)

    def cleanup(self):
        if self._remove_files:
            for rm_file in self._data['files']:
                execute('rm -rf %s' % rm_file, exit=False)


class DummyModule(PipelineModule):

    def __init__(self, infile):
        super().__init__()
        self._data['output'] = infile


class STARMapModule(CmdPipelineModule):

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


class FastQCModule(CmdPipelineModule):

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


class FastXStatisticsModule(CmdPipelineModule):

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


class DuplicateRemovalModule(CmdPipelineModule):

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


class ClippyAdapterClippingModule(CmdPipelineModule):

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


class LafugaAdapterClippingModule(CmdPipelineModule):

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


class LafugaQualityTrimmingModule(CmdPipelineModule):

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


class FastxQualityTrimmingModule(CmdPipelineModule):

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


class LafugaQualityFilterModule(CmdPipelineModule):

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


class SortIndexModule(CmdPipelineModule):

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


class PileupModule(CmdPipelineModule):

    def prepare(self, bam_file, output_dir, prefix, cfg_dict):
        general_cfg = cfg_dict['general']
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


class BamPPModule(CmdPipelineModule):

    def prepare(self, bam_file, output_dir, prefix, cfg_dict):
        read_cfg = cfg_dict['reads']
        pp_cfg = cfg_dict['PostProcessing']
        out_bam_file = os.path.join(output_dir, prefix + '_pp.bam')
        pp_plot_dir = os.path.join(output_dir, 'bam_analysis')
        if not os.path.exists(pp_plot_dir):
            os.makedirs(pp_plot_dir)
        self._data['files'].append(out_bam_file)

        cmd = [
            'stammp-bam_postprocess',
            bam_file,
            out_bam_file,
            pp_plot_dir,
            '--plot_transition_profiles',
            '--min-length %s' % read_cfg['min_len'],
            '--mut_edge_bp %s' % pp_cfg['remove_n_edge_mut'],
            '--max_transitions %s' % pp_cfg['max_mut_per_read'],
        ]
        self._cmds.append(cmd)
        self._data['output'] = out_bam_file


class SoftclipAnalysisModule(CmdPipelineModule):

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


class Pipeline:

    def __init__(self):
        self._jobs = []
        self._current = 0

    def schedule(self, module):
        self._jobs.append(module)

    def __iter__(self):
        self._current = 0
        return self

    def __next__(self):
        if self._current < len(self._jobs):
            self._current += 1
            return self._jobs[self._current - 1]
        else:
            raise StopIteration

    def cur_output(self):
        for job in self._jobs[::-1]:
            if job.get() is not None:
                return job.get()
        return None


if __name__ == '__main__':
    run()
