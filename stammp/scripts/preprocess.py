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
import glob

from stammp.utils import prepare_output_dir, execute
from stammp.scripts.utils.plot_mutation_statistics import create_transition_plots
from stammp.scripts.filter_edge_mutations import filter_edge_mutations


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
        print('%s #### %s ####' % (time_format, tool), flush=True)
        print()

    def _cmd(cmd, exit=True, verbose=verbose):
        return execute(cmd, exit, verbose)

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

    genome_index = general_cfg['genomeindex']
    genome_index_glob = "%s*" % genome_index
    if len(glob.glob(genome_index_glob)) == 0:
        print('genome index %r does not exist. Please check the configuration file' % genome_index,
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

    # check if the mapping method is valid
    if pipeline_cfg['mapper'] not in ['bowtie', 'STAR']:
        print('unknown mapper %r. Please check the configuration file' % pipeline_cfg['mapper'],
              file=sys.stderr)
        sys.exit(1)

    # list of files that can be cleaned up at the end
    intermediate_files = []

    # Run pipeline
    fastq_file = inputfile

    # FastQC analysis of raw data
    if pipeline_cfg.getboolean('fastqc_statistics'):
        print_tool_header('FastQC analysis of raw data')

        fastqc_raw_dir = os.path.join(outputdir, 'fastQC/raw')
        prepare_dir_or_die(fastqc_raw_dir)

        adapter_file = os.path.join(outputdir, 'fastQC/tmp_adapters.txt')
        intermediate_files.append(adapter_file)

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
        intermediate_files.append(qual_stats_file)

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
        intermediate_files.append(rmdup_file)
        dup_rm_toks = [
            'stammp-removePCRduplicates',
            fastq_file,
            rmdup_file,
        ]
        stdout, stderr = _cmd(dup_rm_toks)
        print(stdout)
    else:
        rmdup_file = fastq_file

    # adapter removal
    bc_5prime = read_cfg.getint('bc_5prime')
    clipping_method = pipeline_cfg['adapter_clipping']
    adapter_clipped_file = os.path.join(outputdir, prefix + '_adapter.clipped')
    intermediate_files.append(adapter_clipped_file)

    if clipping_method == 'lafuga':
        print_tool_header('adapter clipping [lafuga]')
        clipper_cfg = config['lafugaAdapterClipper']
        adap_trim_bc_file = os.path.join(outputdir, prefix + '_adapter_bc.clipped')
        adapter_clip_toks = [
            'java -Xmx4g',
            '-jar %s' % os.path.join(scriptPath, 'AdaptorClipper.jar'),
            '-i %s' % rmdup_file,
            '-o %s' % adap_trim_bc_file,
            '-a %s,%s' % (adapter5prime, adapter3prime),
            '-m %s' % (read_cfg.getint('min_len') + bc_5prime),
            '-seed %s' % clipper_cfg['seed'],
        ]
        stdout, stderr = _cmd(adapter_clip_toks)
        print(stdout)

        if bc_5prime > 0:
            intermediate_files.append(adap_trim_bc_file)
            bc_trim_toks = [
                'fastx_trimmer',
                '-i %s' % adap_trim_bc_file,
                '-o %s' % adapter_clipped_file,
                '-f %s' % (bc_5prime + 1),  # first base to keep
            ]
            _cmd(bc_trim_toks)
        else:
            adapter_clipped_file = adap_trim_bc_file

    elif clipping_method == 'clippy':
        print_tool_header('adapter clipping [clippy]')
        clipper_cfg = config['clippyAdapterClipper']

        adapter_clip_toks = [
            'stammp-adapter_clipper',
            rmdup_file,
            adapter_clipped_file,
            adapter5prime,
            adapter3prime,
            '--clip_len %s' % clipper_cfg['clip_len'],
            '--min_len %s' % read_cfg['min_len'],
            '--nt_barcode_5prime %s' % bc_5prime,
            '--verbose'
        ]
        stdout, stderr = _cmd(adapter_clip_toks)
        print(stdout)
    else:
        adapter_clipped_file = rmdup_file

    # quality trimming
    qualtrim_file = os.path.join(outputdir, prefix + '_qtrim.fastq')
    qual_trim_method = pipeline_cfg['quality_trimming']

    # trimming [lafuga]
    if qual_trim_method == 'lafuga':
        intermediate_files.append(qualtrim_file)
        print_tool_header('Quality trimming [lafuga]')
        qualtrim_cfg = config['lafugaQualityTrimmer']

        qualtrim_toks = [
            'java -Xmx4g',
            '-jar %s' % os.path.join(scriptPath, 'FastqQualityFilter.jar'),
            '-i %s' % adapter_clipped_file,
            '-o %s' % qualtrim_file,
            '-m %s' % read_cfg['min_len'],
            '-q %s' % qualtrim_cfg['q'],
            '-3',
            '-5',
        ]
        stdout, stderr = _cmd(qualtrim_toks)
        print(stdout)

    # trimming fastx
    elif qual_trim_method == 'fastx':
        intermediate_files.append(qualtrim_file)
        print_tool_header('Quality trimming [FastX]')
        qualfil_cfg = config['fastxQualityTrimmer']

        qualtrim_toks = [
            'fastq_quality_trimmer',
            '-i %s' % adapter_clipped_file,
            '-o %s' % qualtrim_file,
            '-t %s' % qualfil_cfg['q'],
            '-l %s' % read_cfg['min_len'],
            fx_Q33
        ]
        _cmd(qualtrim_toks)
    else:
        qualtrim_file = adapter_clipped_file

    # polyA removal
    if pipeline_cfg.getboolean('polyA_clipping'):
        print_tool_header('polyA clipping')

        polyA_clipped_file = os.path.join(outputdir, prefix + '_polyA.clipped')
        intermediate_files.append(polyA_clipped_file)
        polyA_cfg = config['polyAClipper']
        polyA_toks = [
            'java',
            '-jar %s' % os.path.join(scriptPath, 'ClipPolyA.jar'),
            '-i %s' % qualtrim_file,
            '-o %s' % polyA_clipped_file,
            '-c %s' % '/dev/null',  # file of all clipped reads
            '-s %s' % polyA_cfg['min_len'],
        ]
        stdout, stderr = _cmd(polyA_toks)
        print(stdout)
    else:
        polyA_clipped_file = qualtrim_file

    # Quality filtering
    qualfil_methods = pipeline_cfg['quality_filtering'].split(',')
    if 'lafuga' in qualfil_methods:
        print_tool_header('Quality filtering [lafuga]')
        qualfil_lafuga_file = os.path.join(outputdir, prefix + '_qfil_lafuga.fastq')
        intermediate_files.append(qualfil_lafuga_file)
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
        qfiltered_file = os.path.join(outputdir, prefix + '_qfiltered.fastq')
        intermediate_files.append(qfiltered_file)
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
        intermediate_files.append(qual_stats_file)
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

    sorted_bam_file = os.path.join(outputdir, prefix + '_sorted.bam')
    if pipeline_cfg['mapper'] == 'bowtie':
        # Bowtie mapping and SAM -> BAM -> mPileup conversion
        print_tool_header('Mapping filtered reads [bowtie]')
        bowtie_cfg = config['bowtie']
        sam_file = os.path.join(outputdir, prefix + '.sam')
        intermediate_files.append(sam_file)
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
            genome_index,
            qfiltered_file,
            '> %s' % sam_file,
        ]
        extra_flags = bowtie_cfg['extra_flags'].split(',')
        bowtie_toks.extend(extra_flags)
        stdout, stderr = _cmd(bowtie_toks)
        print(stderr)

        # SAM --> BAM conversion
        bam_file = os.path.join(outputdir, prefix + '.bam')
        intermediate_files.append(bam_file)
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
        st_sort_toks = [
            'samtools',
            'sort',
            '-@ %s' % general_cfg['n_threads'],
            bam_file,
            '-o %s' % sorted_bam_file,
        ]
        _cmd(st_sort_toks)

    elif pipeline_cfg['mapper'] == 'STAR':
        print_tool_header('Mapping filtered reads [STAR]')
        star_output_dir = os.path.join(outputdir, 'STAR_output')
        if not os.path.exists(star_output_dir):
            os.makedirs(star_output_dir)
        star_output_prefix = os.path.join(star_output_dir, prefix + '_')
        sorted_bam_file = star_output_prefix + 'Aligned.sortedByCoord.out.bam'
        star_toks = [
            'STAR',
            '--readFilesIn %s' % qfiltered_file,
            '--genomeDir %s' % genome_index,
            '--outFilterMultimapNmax %s' % 1,
            '--outFilterMismatchNmax %s' % 1,
            '--outSAMtype BAM SortedByCoordinate',
            '--outFileNamePrefix %s' % star_output_prefix,
            '--runThreadN %s' % general_cfg['n_threads'],
            # I think we cannot allow local alignments because the reads are generally
            # rather short. Remaining adapter sequences probably frequently will be
            # clipped incompletely which leads to mismatched reads.
            '--alignEndsType EndToEnd',
            # we do not want to see insertions or deletions
            '--scoreDelOpen -10000',
            '--scoreInsOpen -10000',
            # remark: maybe we should disallow reads with mismatches at the end.
        ]
        stdout, stderr = _cmd(star_toks)
        print(stdout)

    intermediate_files.append(sorted_bam_file)

    calmd_bam_file = os.path.join(outputdir, prefix + '_calmd.bam')
    calmd_toks = [
        'samtools',
        'calmd',
        sorted_bam_file,
        genome_fasta_path,
        '-b',
        '> %s' % calmd_bam_file
    ]
    _cmd(calmd_toks)

    st_index_toks = [
        'samtools',
        'index',
        calmd_bam_file,
    ]
    _cmd(st_index_toks)

    # filter mutations at the 5' and 3' end of the read
    edge_bp = read_cfg.getint('ignore_edge_mut_bp')
    if edge_bp > 0:
        print_tool_header('Filtering outermost mutations')
        calmd_mutfil_file = os.path.join(outputdir, prefix + '_fil%s.bam' % edge_bp)
        filter_edge_mutations(calmd_bam_file, calmd_mutfil_file, edge_bp)
        intermediate_files.append(calmd_bam_file)
    else:
        calmd_mutfil_file = calmd_bam_file

    # plot transition profiles
    print_tool_header('Analyzing mutation profiles')
    tr_plot_dir = os.path.join(outputdir, 'mapping_plots')
    if not os.path.exists(tr_plot_dir):
        os.makedirs(tr_plot_dir)
    create_transition_plots(calmd_bam_file, tr_plot_dir)

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
        calmd_mutfil_file,
        '> %s' % pileup_file,
    ]
    _cmd(pileup_toks)

    print_tool_header('PreProcess completed')
    if general_cfg.getboolean('rmTemp'):
        print('\tRemoving temporary files...')
        for tmpfile in intermediate_files:
            assert tmpfile != fastq_file
            _cmd('rm -f %s' % tmpfile, exit=False)


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
