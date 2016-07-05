.. _ref_pre-processing:

Pre-Processing
==============
Fastq files obtained from sequencing have to be modified prior to subsequent analysis. If you have more than one protein per sequencing lane you have to demultiplex your data first to get one fastq file per factor. 

The pre-processing script is a wrapper for 3rd party programs to 
    * obtain fastq quality statistics [fastqc and/or fastx]
    * remove adapter sequences
    * remove low quality reads [fastx]
    * map sequences to the genome [bowtie]
    * convert sam->bam->pileup [samtools]

Before running the pre-processing script you have to adjust the pre-process configuration file to your specific experiment (:ref:`ref_pre-processing_config`)

The resulting pileup file of the pre-process serves as input for the following :ref:`ref_binding-site-detection`.

.. _ref_pre-processing_names:

How to name your chromosomes, files and genome identifiers
----------------------------------------------------------
Data analysis of PAR-CLIP data depends on using the identical names for all connected analysis steps. If there are any differences in the chromosome names of the chosen bowtie-index and genome compared to the annotations for post-processing the post-process steps and the normalization won't be correct or won't work at all. The easiest way to get this without any errors is to take a look how names are chosen in the tutorial data. Make sure that you use the same names for the mutliple fasta for bowtie, the filenames of the directory which contains one fasta per chromosome and the annotation files you'll use.

If you map your data back to the genome, bowtie uses the chromosome names as used in the mutiple fasta file which was the input for the index build. In the examplary yeast data you will see the following::
    
    $ grep -P '>' genome.fa
    >chrI
    >chrII
    >chrIII
    >chrIV
    >chrV
    >chrVI
    >chrVII
    >chrVIII
    >chrIX
    >chrX
    >chrXI
    >chrXII
    >chrXIII
    >chrXIV
    >chrXV
    >chrXVI
    >chrMT

These identifiers will be used in the resulting pileup file. Due to the fact, that we will look for binding sites in these files these identifiers will remain. In order to normalize properly and to be able to use post-process scripts the pileups you will use for data normalization have have the same chromosome identifiers. The same is true for the annotations.

.. _ref_pre-processing_config:

:mod:`preprocess.config` --- advanced parameter settings
--------------------------------------------------------
The basic options of the configuration file have to be adjusted to your specific experiment and organism. An example configuration file is located in the stammp directory at::
    
    ... /stammp/scripts/preprocess.config
    
Make a copy of this file and edit the options to your needs.

Lines starting with **#** are considered as comments. Inline comments start with an **;** so everything after **;** is considered as a comment. All options are **key = value** pairs. At first, it is sufficient to take a look at the options below **[basic.options]**. 

.. warning:: Make sure, that you changed the examplary options below **[basic.options]** according to your experiment.

After setting up all options you can pre-process your data::
    
    $ stammp-preprocess /path/to/input.fastq /path/output/ prefix /path/to/configfile
    

If you like, you can also adjust detailed options of used software packages.

:mod:`~stammp.scripts.preprocess`
---------------------------------
.. automodule:: stammp.scripts.preprocess
   :members:

