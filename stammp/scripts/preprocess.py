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

from stammp.utils import prepare_output_dir


def main(inputfile, outputdir, prefix, configfile, verbose):
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

    # define helper functions

    def print_tool_header(tool):
        time_format = time.strftime('[%Y-%m-%d %H:%M:%S]')
        print()
        print('%s #### %s ####' % (time_format, tool))
        print()

    def _cmd(cmd, exit=True, verbose=verbose):
        try:
            if isinstance(cmd, list):
                cmd = ' '.join(cmd)

            if verbose:
                print()
                print('\t' + cmd)
                print()

            proc = subprocess.Popen(args=cmd, shell=True, stderr=subprocess.PIPE,
                                    stdout=subprocess.PIPE, universal_newlines=True)
            retcode = proc.wait()
            stdout, stderr = proc.communicate()
            if retcode != 0:
                print(stderr, file=sys.stderr)
                raise Exception
            return stdout, stderr
        except:
            print('Error at:\n %r \n' % cmd, file=sys.stderr)
            if exit:
                sys.exit(1)

    def prepare_dir_or_die(dir_path):
        try:
            prepare_output_dir(dir_path)
        except ValueError as e:
            print('Error while preparing output directories: %s' % e, file=sys.stderr)
            sys.exit(1)

    config = configparser.ConfigParser(inline_comment_prefixes=';')
    config.read(configfile)

    cur_dir = os.path.dirname(os.path.realpath(__file__))
    scriptPath = os.path.join(cur_dir, 'utils')

    general_cfg = config['general']
    pipeline_cfg = config['pipeline']
    read_cfg = config['reads']

    prepare_dir_or_die(outputdir)

    if config.getboolean('reads', 'fx_Q33'):
        fx_Q33 = '-Q33'
    else:
        fx_Q33 = ''

    bowtie_index = general_cfg['bowtieindex']
    bt_index_glob = "%s*" % bowtie_index
    if len(glob.glob(bt_index_glob)) == 0:
        print('bowtie index %r does not exist. Please check the configuration file' % bowtie_index,
              file=sys.stderr)
        sys.exit(1)

    genome_fasta_path = general_cfg['genomefasta']
    if not os.path.isfile(genome_fasta_path):
        print('genome fasta file %r does not exist. '
              'Please check the configuration file' % genome_fasta_path,
              file=sys.stderr)
        sys.exit(1)

    adapter5prime = general_cfg['adapter5prime']
    adapter3prime = general_cfg['adapter3prime']

    # Run pipeline
    fastq_file = inputfile

    # FastQC analysis of raw data
    if pipeline_cfg.getboolean('fastqc_statistics'):
        print_tool_header('FastQC analysis of raw data')

        fastqc_raw_dir = os.path.join(outputdir, 'fastQC/raw')
        prepare_dir_or_die(fastqc_raw_dir)

        adapter_file = os.path.join(outputdir, 'fastQC/tmp_adapters.txt')
        with open(adapter_file, 'w') as fc:
            print('adapter5prime\t %s' % adapter5prime, file=fc)
            print('adapter3prime\t %s' % adapter3prime, file=fc)

        cmd_tokens = [
            'fastqc',
            '-o %s' % fastqc_raw_dir,
            '-f fastq',
            '--threads %s' % general_cfg['n_threads'],
            '--kmers %s' % config['fastQC']['kmers'],
            '--adapters %s' % adapter_file,
            '-d %s' % fastqc_raw_dir,  # temp directory
            fastq_file,
        ]
        extra_flags = config['fastQC']['extra_flags'].split(',')
        cmd_tokens.extend(extra_flags)

        _cmd(cmd_tokens)

    # FastX-toolkit analysis of raw data
    if pipeline_cfg.getboolean('fastx_statistics'):
        print_tool_header('FastX toolkit analysis of raw data')

        fastx_raw_dir = os.path.join(outputdir, 'fastXstats/raw')
        prepare_dir_or_die(fastx_raw_dir)

        qual_stats_file = os.path.join(fastx_raw_dir, 'fastxstats_raw.txt')
        qual_stats_toks = [
            'fastx_quality_stats',
            '-i %s' % fastq_file,
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

    # Remove duplicates
    if pipeline_cfg.getboolean('remove_duplicates'):
        print_tool_header('Remove duplicated reads')
        rmdup_file = os.path.join(outputdir, prefix + '_nodup.fastq')
        dup_rm_toks = [
            'stammp-removePCRduplicates',
            fastq_file,
            rmdup_file,
        ]
        stdout, stderr = _cmd(dup_rm_toks)
        print(stdout)
    else:
        rmdup_file = fastq_file

    # Trim barcodes
    bc_5prime = read_cfg.getint('bc_5prime')
    if bc_5prime > 0:
        print_tool_header('Trim barcode')
        bc_trim_file = os.path.join(outputdir, prefix + '_trim.fastq')
        bc_trim_toks = [
            'fastx_trimmer',
            '-i %s' % rmdup_file,
            '-o %s' % bc_trim_file,
            '-f %s' % (bc_5prime + 1),  # first base to keep
            bc_trim_file,
        ]
        _cmd(bc_trim_toks)
    else:
        bc_trim_file = rmdup_file

    # quality trimming
    qualtrim_file = os.path.join(outputdir, prefix + '_qtrim.fastq')
    qual_trim_method = pipeline_cfg['quality_trimming']

    # trimming [lafuga]
    if qual_trim_method == 'lafuga':
        print_tool_header('Quality trimming [lafuga]')
        qualtrim_cfg = config['lafugaQualityTrimmer']

        qualtrim_toks = [
            'java -Xmx4g',
            '-jar %s' % os.path.join(scriptPath, 'FastqQualityFilter.jar'),
            '-i %s' % bc_trim_file,
            '-o %s' % qualtrim_file,
            '-m %s' % read_cfg['min_len'],
            '-q %s' % qualtrim_cfg['q'],
            '-3',
            '-5',
        ]
        _cmd(qualtrim_toks)

    # trimming fastx
    elif qual_trim_method == 'fastx':
        print_tool_header('Quality trimming [FastX]')
        qualfil_cfg = config['fastxQualityTrimmer']

        qualtrim_toks = [
            'fastq_quality_trimmer',
            '-i %s' % bc_trim_file,
            '-o %s' % qualtrim_file,
            '-t %s' % qualfil_cfg['q'],
            '-l %s' % read_cfg['min_len'],
            fx_Q33
        ]
        _cmd(qualtrim_toks)
    else:
        qualtrim_file = bc_trim_file

    # adapter removal
    clipping_method = pipeline_cfg['adapter_clipping']
    adapter_clipped_file = os.path.join(outputdir, prefix + '_adapter.clipped')
    if clipping_method == 'lafuga':
        print_tool_header('adapter clipping [lafuga]')
        clipper_cfg = config['lafugaAdapterClipper']
        adap_trim_bc_file = os.path.join(outputdir, prefix + '_adapter_bc.clipped')
        adapter_clip_toks = [
            'java -Xmx4g',
            '-jar %s' % os.path.join(scriptPath, 'AdaptorClipper.jar'),
            '-i %s' % qualtrim_file,
            '-o %s' % adap_trim_bc_file,
            '-a %s,%s' % (adapter5prime, adapter3prime),
            '-m %s' % (read_cfg.getint('min_len') + bc_5prime),
            '-seed %s' % clipper_cfg['seed'],
        ]
        _cmd(adapter_clip_toks)
        if bc_5prime > 0:
            bc_trim_toks = [
                'fastx_trimmer',
                '-i %s' % adap_trim_bc_file,
                '-o %s' % adapter_clipped_file,
                '-f %s' % (bc_5prime + 1),  # first base to keep
                bc_trim_file,
            ]
            _cmd(bc_trim_toks)
        else:
            adapter_clipped_file = adap_trim_bc_file

    elif clipping_method == 'clippy':
        print_tool_header('adapter clipping [clippy]')
        clipper_cfg = config['clippyAdapterClipper']
        adapter_clipped_file = os.path.join(outputdir, prefix + '_adapter.clipped')

        adapter_clip_toks = [
            'stammp-adapter_clipper',
            qualtrim_file,
            adapter_clipped_file,
            adapter5prime,
            adapter3prime,
            '--clip_len %s' % clipper_cfg['clip_len'],
            '--min_len %s' % read_cfg['min_len'],
            '--nt_barcode_5prime %s' % bc_5prime,
            '--verbose'
        ]
        if clipper_cfg.getboolean('aggressive'):
            adapter_clip_toks.append('--aggressive')
        stdout, stderr = _cmd(adapter_clip_toks)
        print(stdout)
    else:
        adapter_clipped_file = qualtrim_file

    # polyA removal
    if pipeline_cfg.getboolean('polyA_clipping'):
        print_tool_header('polyA clipping')

        polyA_clipped_file = os.path.join(outputdir, prefix + '_polyA.clipped')
        polyA_cfg = config['polyAClipper']
        polyA_toks = [
            'java',
            '-jar %s' % os.path.join(scriptPath, 'ClipPolyA.jar'),
            '-i %s' % adapter_clipped_file,
            '-o %s' % polyA_clipped_file,
            '-c %s' % '/dev/null',  # file of all clipped reads
            '-s %s' % polyA_cfg['min_len'],
        ]
        stdout, stderr = _cmd(polyA_toks)
        print(stdout)
    else:
        polyA_clipped_file = adapter_clipped_file

    # Quality filtering
    qualfil_methods = pipeline_cfg['quality_filtering'].split(',')
    qfiltered_file = os.path.join(outputdir, prefix + '_qfiltered.fastq')
    if 'lafuga' in qualfil_methods:
        print_tool_header('Quality filtering [lafuga]')
        qualfil_lafuga_file = os.path.join(outputdir, prefix + '_qfil_lafuga.fastq')
        qualfil_cfg = config['lafugaQualityFilter']

        qualfil_toks = [
            'java -Xmx4g',
            '-jar %s' % os.path.join(scriptPath, 'FastqQualityFilter.jar'),
            '-i %s' % polyA_clipped_file,
            '-o %s' % qualfil_lafuga_file,
            '-m %s' % read_cfg['min_len'],
            '-q %s' % qualfil_cfg['q'],
        ]
        if qualfil_cfg.getboolean('chastity'):
            qualfil_toks.append('-c Y')
        if qualfil_cfg.getboolean('n'):
            qualfil_toks.append('-n')
        _cmd(qualfil_toks)
    else:
        qualfil_lafuga_file = polyA_clipped_file

    if 'fastx' in qualfil_methods:

        print_tool_header('Quality filtering [fastx]')
        qualfil_cfg = config['fastxQualityFilter']
        qualfil_toks = [
            'fastq_quality_filter',
            '-i %s' % qualfil_lafuga_file,
            '-o %s' % qfiltered_file,
            '-q %s' % qualfil_cfg['q'],
            '-p %s' % qualfil_cfg['p'],
            fx_Q33,
        ]
        _cmd(qualfil_toks)
    else:
        qfiltered_file = qualfil_lafuga_file

    # FastQC analysis of filtered data
    if pipeline_cfg.getboolean('fastqc_statistics'):
        print_tool_header('FastQC analysis of filtered data')

        fastqc_filtered_dir = os.path.join(outputdir, 'fastQC/filtered')
        prepare_dir_or_die(fastqc_filtered_dir)

        fastqc_fil_toks = [
            'fastqc',
            '-o %s' % fastqc_filtered_dir,
            '-f fastq',
            '--threads %s' % general_cfg['n_threads'],
            '--kmers %s' % config['fastQC']['kmers'],
            '--adapters %s' % adapter_file,
            '-d %s' % fastqc_filtered_dir,
            qfiltered_file,
        ]
        extra_args = config['fastQC']['extra_flags'].split(',')
        fastqc_fil_toks.extend(extra_args)
        _cmd(fastqc_fil_toks)

    # FastX-toolkit analysis of filtered data
    if pipeline_cfg.getboolean('fastx_statistics'):
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
    bowtie_cfg = config['bowtie']
    sam_file = os.path.join(outputdir, prefix + '.sam')
    bowtie_toks = [
        'bowtie',
        '-t',
        '-q',
        '--threads %s' % general_cfg['n_threads'],
        '-S',
        '-nohead',
        '-v %s' % bowtie_cfg['v'],
        '-y',
        '-m %s' % bowtie_cfg['m'],
        '--best',
        '--strata',
        bowtie_index,
        qfiltered_file,
        '> %s' % sam_file,
    ]
    extra_flags = bowtie_cfg['extra_flags'].split(',')
    bowtie_toks.extend(extra_flags)
    stdout, stderr = _cmd(bowtie_toks)
    print(stderr)

    # SAM --> BAM conversion
    bam_file = os.path.join(outputdir, prefix + '.bam')
    print_tool_header('SAM --> BAM')
    sam2bam_toks = [
        'samtools',
        'view',
        '-@ %s' % general_cfg['n_threads'],
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
        '-@ %s' % general_cfg['n_threads'],
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
    if general_cfg.getboolean('rmTemp'):
        print('\tLet\'s remove temporary files!')
        if pipeline_cfg.getboolean('remove_duplicates'):
            _cmd('rm -f %s' % rmdup_file, exit=False)
        if bc_5prime:
            _cmd('rm -f %s' % bc_trim_file, exit=False)
        if pipeline_cfg['quality_trimming'] in ('lafuga', 'fastx'):
            _cmd('rm -f %s' % qualtrim_file, exit=False)
        if clipping_method == 'lafuga':
            if bc_5prime > 0:
                _cmd('rm -f %s' % adap_trim_bc_file, exit=False)
            _cmd('rm -f %s' % adapter_clipped_file, exit=False)
        elif clipping_method == 'clippy':
            _cmd('rm -f %s' % adapter_clipped_file, exit=False)
        if pipeline_cfg.getboolean('polyA_clipping'):
            _cmd('rm -f %s' % polyA_clipped_file, exit=False)
        if 'lafuga' in qualfil_methods:
            _cmd('rm -f %s' % qualfil_lafuga_file, exit=False)
        if 'fastx' in qualfil_methods:
            _cmd('rm -f %s' % qfiltered_file, exit=False)
        _cmd('rm -f %s' % sam_file, exit=False)
        _cmd('rm -f %s' % bam_file, exit=False)


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
    parser.add_argument('--verbose', '-v', help='more verbose output',
                        action='store_true')
    args = parser.parse_args()

    if not os.path.isfile(args.inputfile):
        print('Input fastq file: %r does not exist' % args.inputfile)
        sys.exit(1)
    if not os.path.isfile(args.configfile):
        print('Config file: %r does not exist' % args.configfile)
        sys.exit(1)
    main(args.inputfile, args.outputdir, args.prefix, args.configfile, args.verbose)


if __name__ == '__main__':
    run()
