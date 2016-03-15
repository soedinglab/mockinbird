#! /usr/bin/python3
"""
Wrapper to convert raw fastq files from sequencing files to mpileup files. A
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
import time
import configparser
import subprocess
import glob

from stammp.utils import native_wordcount as wccount
from stammp.utils import prepare_output_dir


def _cmd(cmd, exit=True):
    try:
        if isinstance(cmd, list):
            cmd = ' '.join(cmd)

        proc = subprocess.Popen(args=cmd, shell=True, stderr=subprocess.PIPE,
                                stdout=subprocess.PIPE, universal_newlines=True)
        retcode = proc.wait()
        stdout, stderr = proc.communicate()
        if retcode != 0:
            print(stderr, file=sys.stderr)
            raise Exception
    except:
        print('Error at:\n %r \n' % cmd, file=sys.stderr)
        if exit:
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

    cur_dir = os.path.dirname(os.path.realpath(__file__))
    scriptPath = os.path.join(cur_dir, 'utils')

    prepare_dir_or_die(outputdir)

    if config['basic.options']['fx_Q33'] == 'Y':
        fx_Q33 = ' -Q33'
    else:
        fx_Q33 = ''

    bowtie_index = config['basic.options']['bowtieindex']
    bt_index_glob = "%s*" % bowtie_index
    if len(glob.glob(bt_index_glob)) == 0:
        print('bowtie index %r does not exist. Please check the configuration file' % bowtie_index,
              file=sys.stderr)
        sys.exit(1)

    genome_fasta_path = config['basic.options']['genomefasta']
    if not os.path.isfile(genome_fasta_path):
        print('genome fasta file %r does not exist. '
              'Please check the configuration file' % genome_fasta_path,
              file=sys.stderr)
        sys.exit(1)

    adapter5prime = config['basic.options']['adapter5prime']
    adapter3prime = config['basic.options']['adapter3prime']

    def print_tool_header(tool):
        time_format = time.strftime('[%Y-%m-%d %H:%M:%S]')
        print()
        print('%s #### %s ####' % (time_format, tool))
        print()

    # FastQC analysis of raw data
    if config['fastQC']['fqc_use'] == 'Y':
        print_tool_header('FastQC analysis of raw data')

        fastqc_raw_dir = os.path.join(outputdir, 'fastQC/raw')
        prepare_dir_or_die(fastqc_raw_dir)

        adapter_file = os.path.join(outputdir, 'fastQC/tmp_adapters.txt')
        with open(adapter_file, 'w') as fc:
            print('adapter5prime\t %s' % adapter5prime, file=fc)
            print('adapter3prime\t %s' % adapter3prime, file=fc)

        cmd_tokens = [
            'fastqc', '-o %s' % fastqc_raw_dir,
            '-f fastq',
            '--threads %s' % config['fastQC']['threads'],
            '--kmers %s' % config['fastQC']['kmers'],
            '--adapters %s' % adapter_file,
            '-d %s' % fastqc_raw_dir,  # temp directory
            inputfile,
        ]
        _cmd(cmd_tokens)

    # FastX-toolkit analysis of raw data
    if config['fastXstatistics']['use'] == 'Y':
        print_tool_header('FastX toolkit analysis of raw data')

        fastx_raw_dir = os.path.join(outputdir, 'fastXstats/raw')
        prepare_dir_or_die(fastx_raw_dir)

        qual_stats_file = os.path.join(fastx_raw_dir, 'fastxstats_raw.txt')
        qual_stats_toks = [
            'fastx_quality_stats',
            '-i %s' % inputfile,
            '-o %s' % qual_stats_file,
            fx_Q33,
        ]
        _cmd(qual_stats_toks)

        qual_bp_toks = [
            'fastq_quality_boxplot_graph.sh',
            '-i %s' % qual_stats_file,
            '-o %s' % os.path.join(fastx_raw_dir, prefix + '_raw_quality.png'),
            '-t %s' % prefix,  # title
        ]
        _cmd(qual_bp_toks)

        nucdistr_toks = [
            'fastx_nucleotide_distribution_graph.sh',
            '-i %s' % qual_stats_file,
            '-o %s' % os.path.join(fastx_raw_dir, prefix + '_raw_nuc.png'),
            '-t %s' % prefix,  # title
        ]
        _cmd(nucdistr_toks)

    # 5prime adapter removal
    print_tool_header('5prime adapter removal')
    line_count = wccount(inputfile)
    print('\tTotal raw reads: %s' % (line_count // 4))

    count_5prime_toks = ['grep', '-c', adapter5prime, inputfile]
    n_5prime_adapter = subprocess.getoutput(' '.join(count_5prime_toks))
    print('\tReads containing the given 5prime adapter %s: %s'
          % (adapter5prime, n_5prime_adapter))

    adapter_clipped_file = os.path.join(outputdir, prefix + '_5prime_adapter.clipped')
    adapter_clip_toks = [
        'stammp-remove5primeAdapter',
        inputfile,
        adapter_clipped_file,
        '--seed %s' % config['remove5primeAdapter']['rm5_seed'],
        '--adapter %s' % adapter5prime,
        '--barcode %s' % config['remove5primeAdapter']['rm5_barcode'],
    ]
    if config['remove5primeAdapter']['rm5_strict'].upper() == 'Y':
        adapter_clip_toks.append('--strict')
    if config['remove5primeAdapter']['rm5_clipanywhere'].upper() == 'Y':
        adapter_clip_toks.append('--clipanywhere')
    print('\t5prime adapter sequence is being clipped ...')
    _cmd(adapter_clip_toks)

    # Random barcode removal
    rndbc_clipped_file = os.path.join(outputdir, prefix + '_rndBC_clipped.fastq')
    fx_use = config['fastxTrimmer']['fx_use'].upper()
    rm_dup = config['removePCRduplicates']['rmPCR_use'].upper()
    if fx_use == 'Y' and rm_dup != 'Y':
        print_tool_header('Random barcode removal')
        trim_toks = [
            'fastx_trimmer',
            '-i %s' % adapter_clipped_file,
            '-o %s' % rndbc_clipped_file,
        ]

        if len(config['fastxTrimmer']['fx_f']) > 0:
            trim_toks.append('-f %s' % config['fastxTrimmer']['fx_f'])
        if len(config['fastxTrimmer']['fx_l']) > 0:
            trim_toks.append('-l %s' % config['fastxTrimmer']['fx_l'])
        _cmd(trim_toks)

    # Remove PCR duplicates by the random Barcode
    clipped_file = os.path.join(outputdir, prefix + '_5prime_adapter.clipped')
    rndbc_clipped_file = os.path.join(outputdir, prefix + '_rndBC_clipped.fastq')
    if fx_use != 'Y' and rm_dup == 'Y':
        print_tool_header('Remove PCR duplicates by the random Barcode')
        dup_rm_toks = [
            'stammp-removePCRduplicates',
            clipped_file,
            rndbc_clipped_file,
        ]
        _cmd(dup_rm_toks)
    if (fx_use == 'Y' and rm_dup == 'Y') or (fx_use != 'Y' and rm_dup != 'Y'):
        print_tool_header('No rnd barcode removal or PCR duplicate removal')
        os.rename(clipped_file, rndbc_clipped_file)

    # Quality Filtering [Graf]
    qfiltered_file = os.path.join(outputdir, prefix + '_qfiltered.fastq')
    if config['qualityFiltering']['qf_use'] == 'Y':
        print_tool_header('Quality filtering [Graf]')
        qual_fil_toks = [
            'java -Xmx4g',
            '-jar %s' % os.path.join(scriptPath, 'FastqQualityFilter.jar'),
            '-i %s' % rndbc_clipped_file,
            '-o %s' % qfiltered_file,
            '-m %s' % config['qualityFiltering']['qf_m'],
            '-q %s' % config['qualityFiltering']['qf_q'],
        ]
        if config['qualityFiltering']['qf_chastity'] == 'Y':
            qual_fil_toks.append('-c Y')
        if config['qualityFiltering']['qf_n'] == 'Y':
            qual_fil_toks.append('-n')
        _cmd(qual_fil_toks)

    # Quality Filtering [FastX]
    if config['fastxQualityFilter']['fxQ_use'] == 'Y':
        print_tool_header('Quality filtering [FastX]')
        qual_fil_toks = [
            'fastq_quality_filter',
            '-i %s' % rndbc_clipped_file,
            '-o %s' % qfiltered_file,
            '-q %s' % config['fastxQualityFilter']['fxQ_q'],
            '-p %s' % config['fastxQualityFilter']['fxQ_p'],
            fx_Q33,
        ]
        _cmd(qual_fil_toks)

    # FastQC analysis of filtered data
    if config['fastQC']['fqc_use'] == 'Y':
        print_tool_header('FastQC analysis of filtered data')

        fastqc_filtered_dir = os.path.join(outputdir, 'fastQC/filtered')
        prepare_dir_or_die(fastqc_filtered_dir)

        fastqc_fil_toks = [
            'fastqc',
            '-o %s' % fastqc_filtered_dir,
            '-f fastq',
            '--threads %s' % config['fastQC']['threads'],
            '--kmers %s' % config['fastQC']['kmers'],
            '--adapters %s' % adapter_file,
            '-d %s' % fastqc_filtered_dir,
            qfiltered_file,
        ]
        _cmd(fastqc_fil_toks)

    # FastX-toolkit analysis of filtered data
    if config['fastXstatistics']['use'] == 'Y':
        print_tool_header('FastX toolkit analysis of filtered data')

        fastx_fil_dir = os.path.join(outputdir, 'fastXstats/filtered')
        prepare_dir_or_die(fastx_fil_dir)

        qual_stats_file = os.path.join(fastx_fil_dir, 'fastxstats_filtered.txt')
        qual_stats_toks = [
            'fastx_quality_stats',
            '-i %s' % qfiltered_file,
            '-o %s' % qual_stats_file,
            fx_Q33,
        ]
        _cmd(qual_stats_toks)

        qual_bp_toks = [
            'fastq_quality_boxplot_graph.sh',
            '-i %s' % qual_stats_file,
            '-o %s' % os.path.join(fastx_fil_dir, prefix + '_filtered_quality.png'),
            '-t %s' % prefix,  # title
        ]
        _cmd(qual_bp_toks)

        nucdistr_toks = [
            'fastx_nucleotide_distribution_graph.sh',
            '-i %s' % qual_stats_file,
            '-o %s' % os.path.join(fastx_raw_dir, prefix + '_filtered_nuc.png'),
            '-t %s' % prefix,  # title
        ]
        _cmd(nucdistr_toks)

    # Bowtie mapping and SAM -> BAM -> mPileup conversion
    print_tool_header('Mapping filtered reads using Bowtie')
    sam_file = os.path.join(outputdir, prefix + '.sam')
    bowtie_toks = [
        'bowtie',
        '-t',
        '-q',
        '--threads %s' % config['bowtie']['bowtie_threads'],
        '-S',
        '-nohead',
        '-v %s' % config['bowtie']['bowtie_v'],
        '-y',
        '-m %s' % config['bowtie']['bowtie_m'],
        '--best',
        '--strata',
        bowtie_index,
        qfiltered_file,
        '> %s' % sam_file,
    ]
    _cmd(bowtie_toks)

    # SAM --> BAM conversion
    bam_file = os.path.join(outputdir, prefix + '.bam')
    print_tool_header('SAM --> BAM')
    sam2bam_toks = [
        'samtools',
        'view',
        '-@ 12',  # number of threads
        '-b',  # output bam
        '-T %s' % genome_fasta_path,
        '-o %s' % bam_file,
        sam_file,
    ]
    _cmd(sam2bam_toks)

    # sorting and indexing of bam file
    print_tool_header('Sorting and indexing of BAM file')
    sorted_bam_file = os.path.join(outputdir, prefix + '_sorted.bam')
    st_sort_toks = [
        'samtools',
        'sort',
        '-@ 12',  # number of threads
        bam_file,
        '-o %s' % sorted_bam_file,
    ]
    _cmd(st_sort_toks)

    st_index_toks = [
        'samtools',
        'index',
        sorted_bam_file,
    ]
    _cmd(st_index_toks)

    # generating pileup file
    print_tool_header('Generating mPileup')
    pileup_file = os.path.join(outputdir, prefix + '.mpileup')
    pileup_toks = [
        'samtools',
        'mpileup',
        '-C 0',  # disable adjust mapping quality
        '-d 100000',  # max depth
        '-q 0',  # minimum alignment quality
        '-Q 0',  # minimum base quality
        '-f %s' % genome_fasta_path,
        sorted_bam_file,
        '> %s' % pileup_file,
    ]
    _cmd(pileup_toks)

    print_tool_header('PreProcess completed')
    if config['basic.options']['rmTemp'] == 'Y':
        print('\tLet\'s remove temporary files!')
        _cmd('rm -f %s' % adapter_clipped_file, exit=False)
        _cmd('rm -f %s' % rndbc_clipped_file, exit=False)
        _cmd('rm -f %s' % sam_file, exit=False)
        _cmd('rm -f %s' % bam_file, exit=False)
        _cmd('rm -f %s' % qfiltered_file, exit=False)


def run():
    parser = argparse.ArgumentParser(
        description=(
            'Wrapper to convert raw fastq files from sequencing files to mpileup files. '
            'A fastq-file is adapterclipped, qualityfiltered, mapped and converted.'
        ),
        epilog='contact: torkler@genzentrum.lmu.de'
    )
    parser.add_argument('inputfile', help='Input fastq file')
    parser.add_argument('outputdir', help='Output directory')
    parser.add_argument('prefix', help='Set prefix for all outputfile')
    parser.add_argument('configfile', help='Config file containing arguments')
    args = parser.parse_args()

    if not os.path.isfile(args.inputfile):
        print('Input fastq file: %r does not exist' % args.inputfile)
        sys.exit(1)
    if not os.path.isfile(args.configfile):
        print('Config file: %r does not exist' % args.configfile)
        sys.exit(1)
    main(args.inputfile, args.outputdir, args.prefix, args.configfile)


def prepare_dir_or_die(dir_path):
    try:
        prepare_output_dir(dir_path)
    except ValueError as e:
        print('Error while preparing output directories: %s' % e, file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    run()
