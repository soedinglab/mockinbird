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
#from stammp.obj.functions import *

def _cmd(cmd, exit=True):
    try:
        proc     = subprocess.Popen(args=cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE )
        retcode  = proc.wait()
        com      = proc.communicate()
        #sys.stdout.write(com[0])
        if retcode != 0:
            sys.stderr.write(str(com[1].decode(encoding="utf-8", errors="strict")))
            raise Exception
    except:
        sys.stderr.write('Error at:\n  '+cmd+'\n\n')
        if exit:
            sys.exit(-1)

def _printConfig(c):
    for section in c.sections():
        print(section)
        for key in c[section]:
            print('\t'+key+'\t:'+c[section][key])

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
    OUTDIR  = outputdir #config['basic.options']['outputdir']
    PREFIX  = prefix
    config = configparser.ConfigParser(inline_comment_prefixes=';')
    config.read(configfile)
    #_printConfig(config)
    scriptPath = os.path.dirname(os.path.realpath(__file__))+'/utils/'
    #scriptPath = '../utils/'
    if os.path.exists(OUTDIR) == False:
        os.makedirs(OUTDIR)
    fx_Q33 = ''
    if config['basic.options']['fx_Q33'] == 'Y':
       fx_Q33 = ' -Q33'
    ##################################################
    ################################################## FastQC analysis of raw data
    ##################################################
    if config['fastQC']['fqc_use'] == 'Y':
        print('')
        print(time.strftime("[%Y-%m-%d %H:%M:%S]"), '##### FastQC analysis of raw data ###########')
        print('')
        if os.path.exists(OUTDIR+'/fastQC/raw/') == False or os.path.exists(OUTDIR+'/fastQC/filtered/') == False:
            os.makedirs(OUTDIR+'/fastQC/raw/')
            os.makedirs(OUTDIR+'/fastQC/filtered/')
        fc = open(OUTDIR+'/fastQC/tmp_adapters.txt','w')
        fc.write('adapter5prime\t'+config['basic.options']['adapter5prime']+'\n')
        fc.write('adapter3prime\t'+config['basic.options']['adapter3prime'])
        fc.close()
        _cmd('fastqc -o '+OUTDIR+'fastQC/raw/ -f fastq --threads '+config['fastQC']['threads']+' --kmers '+config['fastQC']['kmers']+' --adapters '+OUTDIR+'fastQC/tmp_adapters.txt -d '+OUTDIR+'fastQC/raw/ '+inputfile)
    
    ##################################################
    ################################################## FastX-toolkit analysis of raw data
    ##################################################
    if config['fastXstatistics']['use'] == 'Y':
        print('')
        print(time.strftime("[%Y-%m-%d %H:%M:%S]"), '##### FastX toolkit analysis of raw data ###########')
        print('')
        if os.path.exists(OUTDIR+'/fastXstats/raw/') == False or os.path.exists(OUTDIR+'/fastXstats/filtered/') == False:
            os.makedirs(OUTDIR+'/fastXstats/raw/')
            os.makedirs(OUTDIR+'/fastXstats/filtered/')
        _cmd('fastx_quality_stats -i '+inputfile+' -o '+OUTDIR+'/fastXstats/raw/fastxstats_raw.txt'+fx_Q33)
        _cmd('fastq_quality_boxplot_graph.sh -i '+OUTDIR+'/fastXstats/raw/fastxstats_raw.txt -o '+OUTDIR+'/fastXstats/raw/'+PREFIX+'_raw_quality.png -t '+PREFIX)
        _cmd('fastx_nucleotide_distribution_graph.sh -i '+OUTDIR+'/fastXstats/raw/fastxstats_raw.txt -o '+OUTDIR+'/fastXstats/raw/'+PREFIX+'_raw_nuc.png -t '+PREFIX)
    
    ##################################################
    ################################################## 5prime adapter removal
    ##################################################
    print('')
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), '##### 5prime adapter removal ###########')
    print('')
    wc_output = subprocess.check_output(['wc', '-l', inputfile], universal_newlines=True)
    lineno_str, *_ = wc_output.split()
    print('\tTotal raw reads: %s' % (int(lineno_str) // 4))
    print('\tReads containing the given 5prime adapter ['+config['basic.options']['adapter5prime']+']: '+subprocess.check_output(['grep', '-c', config['basic.options']['adapter5prime'], inputfile], universal_newlines=True))
    cmd_string = 'stammp-remove5primeAdapter '+inputfile+' '+OUTDIR+PREFIX+'_5prime_adapter.clipped --seed '+config['remove5primeAdapter']['rm5_seed']+' --adapter '+config['basic.options']['adapter5prime']+' --barcode '+config['remove5primeAdapter']['rm5_barcode']
    if config['remove5primeAdapter']['rm5_strict'].upper() == 'Y':
        cmd_string = cmd_string+' --strict'
    if config['remove5primeAdapter']['rm5_clipanywhere'].upper() == 'Y':
        cmd_string = cmd_string+' --clipanywhere'
    print('\t5prime adapter sequence is being clipped ...')
    _cmd(cmd_string)
    
    ##################################################
    ################################################## Random barcode removal
    ##################################################
    if config['fastxTrimmer']['fx_use'].upper() == 'Y' and config['removePCRduplicates']['rmPCR_use'].upper() != 'Y':
        print('')
        print(time.strftime("[%Y-%m-%d %H:%M:%S]"), '##### Random barcode removal ###########')
        print('')
        cmd_string = 'fastx_trimmer -i '+OUTDIR+PREFIX+'_5prime_adapter.clipped -o '+OUTDIR+PREFIX+'_rndBC_clipped.fastq'
        if len(config['fastxTrimmer']['fx_f']) > 0:
            cmd_string = cmd_string+' -f '+config['fastxTrimmer']['fx_f']
        if len(config['fastxTrimmer']['fx_l']) > 0:
            cmd_string = cmd_string+' -l '+config['fastxTrimmer']['fx_l']
        #_cmd('fastx_trimmer -i '+OUTDIR+PREFIX+'.5prime_adapter.clipped -o '+OUTDIR+PREFIX+'_rndBC_clipped.fastq '+fx_f+' '+fx_l+' '+fx_Q33)
        #_cmd(['fastx_trimmer', '-i '+OUTDIR+PREFIX+'.5prime_adapter.clipped', '-o '+OUTDIR+PREFIX+'_rndBC_clipped.fastq', fx_f, fx_l, fx_Q33])
        _cmd(cmd_string+fx_Q33)
    
    ##################################################
    ################################################## Remove PCR duplicates by the random Barcode
    ##################################################
    if config['fastxTrimmer']['fx_use'].upper() != 'Y' and config['removePCRduplicates']['rmPCR_use'].upper() == 'Y':
        print('')
        print(time.strftime("[%Y-%m-%d %H:%M:%S]"), '##### Remove PCR duplicates by the random Barcode  ###########')
        print('')
        _cmd('stammp-removePCRduplicates '+OUTDIR+PREFIX+'_5prime_adapter.clipped '+OUTDIR+PREFIX+'_rndBC_clipped.fastq')
    if (config['fastxTrimmer']['fx_use'].upper() == 'Y' and config['removePCRduplicates']['rmPCR_use'].upper() == 'Y') or (config['fastxTrimmer']['fx_use'].upper() != 'Y' and config['removePCRduplicates']['rmPCR_use'].upper() != 'Y'):
        print('')
        print(time.strftime("[%Y-%m-%d %H:%M:%S]"), '##### No rnd barcode removal or PCR duplicate removal ###########')
        print('')
        os.rename(OUTDIR+PREFIX+'_5prime_adapter.clipped', OUTDIR+PREFIX+'_rndBC_clipped.fastq')
    
    ##################################################
    ################################################## Quality Filtering [Graf]
    ##################################################
    if config['qualityFiltering']['qf_use'] == 'Y':
        print('')
        print(time.strftime("[%Y-%m-%d %H:%M:%S]"), '##### Quality filtering [Graf] ###########')
        print('')
        cmd_string = 'java -Xmx4g -jar '+scriptPath+'FastqQualityFilter.jar -i '+OUTDIR+PREFIX+'_rndBC_clipped.fastq -o '+OUTDIR+PREFIX+'_qfiltered.fastq -m '+config['qualityFiltering']['qf_m']+' -q '+config['qualityFiltering']['qf_q']
        if config['qualityFiltering']['qf_chastity'] == 'Y':
            cmd_string = cmd_string+' -c Y'
        if config['qualityFiltering']['qf_n'] == 'Y':
            cmd_string = cmd_string+' -n'
        #_cmd('java -Xmx4g -jar '+scriptPath+'FastqQualityFilter.jar -i '+OUTDIR+PREFIX+'_rndBC_clipped.fastq -o '+OUTDIR+PREFIX+'_qfiltered.fastq -m '+config['qualityFiltering']['qf_m']+' -q '+config['qualityFiltering']['qf_q']+' '+qf_n)
        #_cmd(['java', '-Xmx4g', '-jar', scriptPath+'FastqQualityFilter.jar', '-i '+OUTDIR+PREFIX+'_rndBC_clipped.fastq', '-o '+OUTDIR+PREFIX+'_qfiltered.fastq', '-m '+config['qualityFiltering']['qf_m'], ' -q '+config['qualityFiltering']['qf_q'], qf_n])
        _cmd(cmd_strig)
    
    ##################################################
    ################################################## Quality Filtering [FastX]
    ##################################################
    if config['fastxQualityFilter']['fxQ_use'] == 'Y':
        print('')
        print(time.strftime("[%Y-%m-%d %H:%M:%S]"), '##### Quality filtering [FastX] ###########')
        print('')
        _cmd('fastq_quality_filter -i '+OUTDIR+PREFIX+'_rndBC_clipped.fastq -o '+OUTDIR+PREFIX+'_qfiltered.fastq -q '+config['fastxQualityFilter']['fxQ_q']+' -p '+config['fastxQualityFilter']['fxQ_p']+fx_Q33)
    
    ##################################################
    ################################################## FastQC analysis of filtered data
    ##################################################
    if config['fastQC']['fqc_use'] == 'Y':
        print('')
        print(time.strftime("[%Y-%m-%d %H:%M:%S]"), '##### FastQC analysis of filtered data ###########')
        print('')
        _cmd('fastqc -o '+OUTDIR+'fastQC/filtered/ -f fastq --threads '+config['fastQC']['threads']+' --kmers '+config['fastQC']['kmers']+' --adapters '+OUTDIR+'fastQC/tmp_adapters.txt -d '+OUTDIR+'fastQC/filtered/ '+OUTDIR+PREFIX+'_qfiltered.fastq')
    ##################################################
    ################################################## FastX-toolkit analysis of filtered data
    ##################################################
    if config['fastXstatistics']['use'] == 'Y':
        print('')
        print(time.strftime("[%Y-%m-%d %H:%M:%S]"), '##### FastX toolkit analysis of filtered data ###########')
        print('')
        _cmd('fastx_quality_stats -i '+OUTDIR+PREFIX+'_qfiltered.fastq -o '+OUTDIR+'/fastXstats/filtered/fastxstats_filtered.txt'+fx_Q33)
        _cmd('fastq_quality_boxplot_graph.sh -i '+OUTDIR+'/fastXstats/filtered/fastxstats_filtered.txt -o '+OUTDIR+'/fastXstats/filtered/'+PREFIX+'_filtered_quality.png -t '+PREFIX)
        _cmd('fastx_nucleotide_distribution_graph.sh -i '+OUTDIR+'/fastXstats/filtered/fastxstats_filtered.txt -o '+OUTDIR+'/fastXstats/filtered/'+PREFIX+'_filtered_nuc.png -t '+PREFIX)
    
    
    ##################################################
    ################################################## Bowtie mapping and SAM -> BAM -> mPileup conversion
    ##################################################
    print('')
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), '##### Mapping filtered reads using Bowtie 1.1.2 #################')
    print('')
    if config['bowtie']['bowtie_local'] == 'Y':
       _cmd(scriptPath+'bowtie-1.1.2/bowtie -t -q -p '+config['bowtie']['bowtie_threads']+' -S -nohead -v '+config['bowtie']['bowtie_v']+' -y -m '+config['bowtie']['bowtie_m']+' --best --strata '+config['basic.options']['bowtieindex']+' '+OUTDIR+PREFIX+'_qfiltered.fastq > '+OUTDIR+PREFIX+'.sam')
    else:
        _cmd('bowtie -t -q -p '+config['bowtie']['bowtie_threads']+' -S -nohead -v '+config['bowtie']['bowtie_v']+' -y -m '+config['bowtie']['bowtie_m']+' --best --strata '+config['basic.options']['bowtieindex']+' '+OUTDIR+PREFIX+'_qfiltered.fastq > '+OUTDIR+PREFIX+'.sam')
    #os.system(scriptPath+'bowtie-1.1.2/bowtie --version')
    #_cmd('bowtie --version')
    #sys.exit()
    print('')
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), '##### SAM --> BAM ###################################')
    _cmd('samtools view -@ 12 -b -T '+config['basic.options']['genomefasta']+' -o '+OUTDIR+PREFIX+'.bam '+OUTDIR+PREFIX+'.sam')

    print('')
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), '##### Sorting and indexing of BAM file ###########################')
    _cmd('samtools sort -@ 12 '+OUTDIR+PREFIX+'.bam -o '+OUTDIR+PREFIX+'_sorted.bam')
    _cmd('samtools index '+OUTDIR+PREFIX+'_sorted.bam')

    print('')
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), '##### Generating mPileup ############################')
    _cmd('samtools mpileup -C 0 -d 100000 -q 0 -Q 0 -f '+config['basic.options']['genomefasta']+' '+OUTDIR+PREFIX+'_sorted.bam > '+OUTDIR+PREFIX+'.mpileup')

    print('')
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), '##### PreProcess completed ############################')
    if config['basic.options']['rmTemp'] == 'Y':
        print('\tLet\'s remove temporary files!')
        _cmd('rm -f '+OUTDIR+PREFIX+'_5prime_adapter.clipped', exit=False)
        _cmd('rm -f '+OUTDIR+PREFIX+'_rndBC_clipped.fastq', exit=False)
        _cmd('rm -f '+OUTDIR+PREFIX+'.sam', exit=False)
        _cmd('rm -f '+OUTDIR+PREFIX+'.bam', exit=False)
        _cmd('rm -f '+OUTDIR+PREFIX+'_qfiltered.fastq', exit=False)

def run():
    parser = argparse.ArgumentParser(description='Wrapper to convert raw fastq files from sequencing files to mpileup files. A fastq-file is adapterclipped, qualityfiltered, mapped and converted.', epilog="contact: torkler@genzentrum.lmu.de")
    parser.add_argument('inputfile', help='Input fastq file')
    parser.add_argument('outputdir', help='Output directory')
    parser.add_argument('prefix', help='Set prefix for all outputfile')
    parser.add_argument('configfile', help='Config file containing arguments')
    #parser.add_argument('-v','--verbose', dest='verbose', action="store_true", default=False, help='verbose output')
    args = parser.parse_args()
    if os.path.isfile(args.inputfile) == False:
        print('Input fastq file: '+args.inputfile+' does not exist')
        sys.exit(-1)
    if os.path.isfile(args.configfile) == False:
        print('Config file: '+args.configfile+' does not exist')
        sys.exit(-1)
    main(args.inputfile, args.outputdir, args.prefix, args.configfile)

if __name__ == '__main__':
    run()
