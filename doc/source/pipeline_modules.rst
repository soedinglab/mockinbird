Pipeline modules
================


Preprocessing modules
---------------------

Diagnostics
^^^^^^^^^^^

These modules can be used to assess the quality of the PAR-CLIP library.

FastQCModule
""""""""""""
Uses `FastQC <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ to create a diagnostic
report of various quality aspects of the fastq reads.

**Pipeline input:** ``fastq``

**Pipeline output:**

================  =================  ============================================================
Parameter         Default value      Description
================  =================  ============================================================
kmer_length       7                  length of the kmers used for detecting overrepresented kmers
outdir_name       fastQC             name of the output directory that stores the report
extra_flags       []                 list of extra flags passed to the command line tool
================  =================  ============================================================

BamPPModule
"""""""""""
Postprocesses a bam alignment file and creates statistics of mutation rates and mismatch positions.
Can clip mutations from the end of the reads.

**Pipeline input:** ``bam``

**Pipeline output:** ``bam``

====================  =================  ==========================================================
Parameter             Default value      Description
====================  =================  ==========================================================
remove_n_edge_mut     0                  clip mutations that are at the n outermost bases in each
                                         alignment
max_mut_per_read      1                  discard reads that have more mutations than the given
                                         integer
min_base_quality      0                  discard reads that have at least one base with a quality
                                         lower than the given integer
min_avg_ali_quality   20                 discard reads that have an average quality lower than the
                                         given integer
min_mismatch_quality  20                 discard reads with mismatches of quality lower than the
                                         given integer
dump_raw_data         False              dump raw data for diagnostic purposes
outdir_name           bam_analysis       name of the output directory storing diagnostics and plots
====================  =================  ==========================================================

SoftclipAnalysisModule
""""""""""""""""""""""

If the mapper allows soft-clipping this module shows the most frequently clipped sequences.
The output can be used to understand common unclipped contamination. Creates two files
``3prime_clipped.txt`` and ``3prime_clipped.txt``.

Interpretation of the output:

::

        clipped bases: 2 | clipped sequences: 152

        TG                                    128
        CG                                      7
        CT                                      4
        CC                                      3
        TT                                      3
        TC                                      2
        AT                                      1
        CA                                      1
        GG                                      1
        GT                                      1

The header line contains the number of clipped bases and the number of affected alignments.
Below the most commonly clipped sequences are enumerated along with the number of affected
alignments in descending order.


**Pipeline input:** ``bam``

**Pipeline output:**

====================  =================  ==========================================================
Parameter             Default value      Description
====================  =================  ==========================================================
outdir_name           bam_analysis       name of the output directory storing diagnostics and plots
====================  =================  ==========================================================

Duplicate removal
^^^^^^^^^^^^^^^^^

In order to make meaningful biological conclusions, PCR duplicates should be removed after
sequencing.
We offer two modules here. A naive approach implemented in :ref:`dupremov_module` and two modules
using `umi_tools <https://github.com/CGATOxford/UMI-tools/tree/master/umi_tools>`_, applicable
when the library contains unique molecular identifiers (UMIs).

.. _dupremov_module:

DuplicateRemovalModule
""""""""""""""""""""""

A poor man's approach to duplicate removal: this module removes fastq reads with the same sequence.
If the PAR-CLIP library uses UMIs, it is recommended to use the more sophisticated modules based on
`umi_tools <https://github.com/CGATOxford/UMI-tools/tree/master/umi_tools>`_. Please refer to
:ref:`umi_tools_dedup_module` and  :ref:`umi_tools_extract_module`.

**Pipeline input:** ``fastq``

**Pipeline output:** ``fastq``

.. _umi_tools_extract_module:

UmiToolsExtractModule
"""""""""""""""""""""

Together with :ref:`umi_tools_dedup_module` this is the preferred way to deduplicate PAR-CLIP
libraries that use UMIs. This module extracts the UMI and appends it to the read name. As such it
is recommended to run this module as one of the first steps in a pipeline. After mapping the
:ref:`umi_tools_dedup_module` module can be used to remove duplicated reads.

For further information, please also check the documentation of `umi_tools <https://github.com/CGATOxford/UMI-tools/tree/master/umi_tools>`_.

**Pipeline input:** ``fastq``

**Pipeline output:** ``fastq``

.. _umi_tools_dedup_module:

UmiToolsDedupModule
"""""""""""""""""""

Together with :ref:`umi_tools_extract_module` this is the preferred way to deduplicate PAR-CLIP
libraries that use UMIs. This module deduplicates bam files based on extracted UMIs. This module
has to be run after :ref:`umi_tools_extract_module`.

For further information, please also check the documentation of `umi_tools <https://github.com/CGATOxford/UMI-tools/tree/master/umi_tools>`_.

**Pipeline input:** ``bam``

**Pipeline output:** ``bam``

Adapter clipping
^^^^^^^^^^^^^^^^

Adapter clipping is an important step in PAR-CLIP libraries: in case of very small inserts,
the beginning of the 3' adapter is present in the reads.
:ref:`skewer_module` is the module of choice for removing these adapters.

:ref:`clippy_module` clips adapters with less sensitity, but also detects partial 5' adapters. It
can be useful when dealing with libraries that were generated from very low amounts of RNA.

.. _skewer_module:

SkewerAdapterClippingModule
"""""""""""""""""""""""""""
This module uses `skewer <https://github.com/relipmoc/skewer>`_ to trim the 3' adapter from the
PAR-CLIP reads. Additional arguments can be directly passed to skewer's command line call.
Please consult skewer's documentation for a detailed description of all available options.

**Pipeline input:** ``fastq``

**Pipeline output:** ``fastq``

================  =================  ============================================================
Parameter         Default value      Description
================  =================  ============================================================
extra_args        []                 list of extra flags passed to the command line tool
================  =================  ============================================================

.. _clippy_module:

ClippyAdapterClippingModule
"""""""""""""""""""""""""""
This module removes longer traces of adapters by looking for perfect matches of the adapter ends.
It also removes random barcodes and adapters. If you are using :ref:`umi_tools_extract_module` to
extract the UMIs don't forget to set ``clipped_5prime_bc`` to ``True``.

**Pipeline input:** ``fastq``

**Pipeline output:** ``fastq``

=================  =================  ============================================================
Parameter          Default value      Description
=================  =================  ============================================================
clip_len           10                 minimum base pairs required to be detected as adapter sequence
clipped_5prime_bc  False              UMIs already removed from the 5' end
=================  =================  ============================================================

Mapping
^^^^^^^

With `Bowtie <http://bowtie-bio.sourceforge.net/tutorial.shtml>`_ and
`STAR <https://github.com/alexdobin/STAR>`_ we offer two very different alignment strategies.
`STAR` can map spliced reads and use softclipping to remove contaminants automatically. Soft
clipping however prevents from mapping transitions at either end of the alignment. Mappers feed
unmapped reads back in the pipeline and thus can be chained.

STARMapModule
"""""""""""""
`STAR <https://github.com/alexdobin/STAR>`_ is a general purpose RNA-seq data mapper. Unmapped
reads are returned to the pipeline in fastq format. For detailed configuration options, please also
refer to STAR's user manual.

**Pipeline input:** ``fastq``

**Pipeline output:** ``bam``, ``fastq``

===================  =================  ============================================================
Parameter            Default value      Description
===================  =================  ============================================================
genome_index                            path to the directory containing the STAR genome index
n_mismatch           1                  maximum number of allowed mismatches
n_multimap           1                  maximum number of mapping positions. Maps uniquely by
                                        default
allow_soft_clipping  True               enable softclipping
outdir_name          star_out           name of the output directory
extra_flags          []                 additional commandline options passed to STAR
===================  =================  ============================================================

BowtieMapModule
"""""""""""""""
`Bowtie <http://bowtie-bio.sourceforge.net/tutorial.shtml>`_ is a genomic aligner and as such
cannot map spliced reads. Unmapped reads are re-queued in the pipeline.

**Pipeline input:** ``fastq``

**Pipeline output:** ``bam``, ``fastq``

===================  =================  ============================================================
Parameter            Default value      Description
===================  =================  ============================================================
genome_index                            prefix of bowtie's genome index
n_mismatch           1                  maximum number of allowed mismatches
n_multimap           1                  maximum number of mapping positions. Maps uniquely by
                                        default
extra_flags          []                 additional commandline options passed to bowtie
===================  =================  ============================================================


Binding site prediction
^^^^^^^^^^^^^^^^^^^^^^^

We offer two different strategies for predicing binding sites. :ref:`bsfinder_module` does not
require mock information and therefore cannot distinguish background binding from factor specific
binding events.
:ref:`mockinbird_module` is the recommended binding site predictor. It requires parameters trained
on a mock experiment.

.. _bsfinder_module:

BSFinderModule
""""""""""""""

BSFinder calculates p-values by learning a statistical model on non-specific conversion events.
This model cannot distinguish background binding from factor specific binding events.

**Pipeline input:** ``mpileup``

**Pipeline output:** ``table``

===================  =================  ============================================================
Parameter            Default value      Description
===================  =================  ============================================================
pval_threshold       0.005              only sites with a p-value smaller than this are reported
min_cov              2                  minimum coverage of reported binding sites
===================  =================  ============================================================

.. _mockinbird_module:

MockinbirdModule
""""""""""""""""
Mockinbird is the core module that predicts binding sites by harnessing information from a mock
experiments. Its input files are generated by the modules in the :ref:`mockinbird_rel_modules`
section.

**Pipeline input:** ``trtable``, ``mock_model``

**Pipeline output:** ``table``

===================  =================  ============================================================
Parameter            Default value      Description
===================  =================  ============================================================
plot_dir             mockinbird_plots   directory for writing out diagnostic plots
max_k_mock           10                 sites with more specific conversions than ``max_k_mock`` are
                                        discarded
extra_args           []                 additional arguments directly passed to the called script
===================  =================  ============================================================

Miscellaneous
^^^^^^^^^^^^^

A collection of miscellaneous helper modules.

.. _sort_index_module:

SortIndexModule
"""""""""""""""
The SortIndexModule sorts and indexes a bam files using
`samtools <http://samtools.sourceforge.net/>`_. This is generally required before generating a
pileup file.

**Pipeline input:** ``bam``

**Pipeline output:** ``bam``

PileupModule
""""""""""""

Uses `samtools <http://samtools.sourceforge.net/>`_ to create a pileup file. Pileup report
coverage and transitions per genomic base and are the input of our predictors.
Please be aware that the bam file has to be sorted. If in doubt, queue after the
:ref:`sort_index_module`.

**Pipeline input:** ``bam``

**Pipeline output:** ``mpileup``

NormalizationModule
"""""""""""""""""""

The normalization module calculates an occupancy by dividing the number of observed transitions by
the coverage of a reference experiment. The appropriate reference experiment should reflect the
pool of RNA the factor *sees* when *choosing* where to bind. Depending on the binding properties
of the protein of interest, an RNA-seq experiment under PAR-CLIP conditions, PAR-CLIP of the
RNA polymerase or protocols to capture transient binding such as 4SU-seq may be appropriate.

Additionally, SNPs are removed by detecting elevated conversion rates in the normalization pileup
file.

**Pipeline input:** ``table``

**Pipeline output:** ``table``

===================  =================  ============================================================
Parameter            Default value      Description
===================  =================  ============================================================
mut_snp_ratio        0.75               ratio of conversations to coverage in the normalization
                                        pileup for a site being detected as SNP
===================  =================  ============================================================

QuantileCapModule
"""""""""""""""""
Caps the occupancy value at a given quantile. This module can help removing the influence of
outliers on downstream analyses, such as the gene plot.

**Pipeline input:** ``table``

**Pipeline output:** ``table``

===================  =================  ============================================================
Parameter            Default value      Description
===================  =================  ============================================================
max_quantile         0.95               all occupancy values are capped to the value of this
                                        quantile
===================  =================  ============================================================

Table2FastaModule
"""""""""""""""""
Converts a binding site table file to fasta by extracting the genomic sequence around the binding
site.

**Pipeline input:** ``table``

**Pipeline output:** ``fasta``

===================  =================  ============================================================
Parameter            Default value      Description
===================  =================  ============================================================
genome_fasta                            path to the genome fasta file
===================  =================  ============================================================

.. _mockinbird_rel_modules:

Mockinbird related modules
^^^^^^^^^^^^^^^^^^^^^^^^^^
:ref:`mockinbird_module` requires as input a joint dataset of factor of interest and the mock
measurement and parameters estimated on the mock. The following modules can be used to create the
required input data.

.. _prediction_sites_module:

PredictionSitesModule
"""""""""""""""""""""
This module creates a file that contains all sites that are considered in the prediction of binding
sites.
By default these are all genomic sites that have the transition nucleotide on either strand. This
can be restricted by giving gff files of genomic regions of interest.

**Pipeline input:**

**Pipeline output:** ``sites``

=====================  =================  ============================================================
Parameter              Default value      Description
=====================  =================  ============================================================
sites_file                                path to sites file. Will be created if does not exist yet.
                                          Will not be recreated if already existing.
fasta_file                                path to genomic fasta file
gff_file               ''                 gff file for restricting predictions to specific regions
transition_nucleotide  T                  nucleotide that converts in the PAR-CLIP experiment
=====================  =================  ============================================================

.. _mock_table_module:

MockTableModule
"""""""""""""""

Converts the pileup file from the mock experiment to a mock table. Required by the
:ref:`trtable_module`.

**Pipeline input:**

**Pipeline output:** ``mocktable``

=====================  =================  ============================================================
Parameter              Default value      Description
=====================  =================  ============================================================
mock_pileup                               path to the mock pileup file
mock_table                                path to mock table. Will be created if does not exist yet.
                                          Will not be recreated if already existing.
=====================  =================  ============================================================

.. _trtable_module:

TransitionTableModule
"""""""""""""""""""""

Combines mock table and factor pileup file to the so called transition table. Depends on the outputs
of :ref:`mock_table_module`, :ref:`prediction_sites_module`.

**Pipeline input:** ``sites``, ``mock_table``, ``mpileup``

**Pipeline output:** ``trtable``

LearnMockModule
"""""""""""""""

Learns the model parameters from the transition table. Requires inputs from :ref:`trtable_module`
and :ref:`bam_stat_module`.

**Pipeline input:** ``trtable``

**Pipeline output:** ``mock_model``

=====================  =================  ============================================================
Parameter              Default value      Description
=====================  =================  ============================================================
mock_model                                path to the mock model pickle file. Will be created if does
                                          not exist. Will not be recreated if already existing
mock_statistics                           path to the mock bam statistics
n_mixture_components   5                  number of mixture components for fitting the geometric
                                          mixture models
em_iterations          250                number of iterations of the EM algorithm fitting the
                                          geometric mixture model
=====================  =================  ============================================================

.. _bam_stat_module:

BamStatisticsModule
"""""""""""""""""""

This module stores additional information from a ``bam`` file in a json file.
If the predicted sites are restrained by a gff file in :ref:`prediction_sites_module`, the same
gff file should be used for generating the statistics.

**Pipeline input:** ``bam``

**Pipeline output:** ``stat_file``

=====================  =================  ============================================================
Parameter              Default value      Description
=====================  =================  ============================================================
gff_file               ''                 gff file for restricting the prediction sites
=====================  =================  ============================================================

Postprocessing modules
----------------------

Plots
^^^^^

CenterPlotBSModule
""""""""""""""""""

KmerPerPositionModule
"""""""""""""""""""""

TransitionFrequencyModule
"""""""""""""""""""""""""

HeatmapPlotModule
"""""""""""""""""

HeatmapSmallPlotModule
""""""""""""""""""""""


Motif detection
^^^^^^^^^^^^^^^

XXmotifModule
"""""""""""""


Miscellaneous
^^^^^^^^^^^^^

GffFilterModule
"""""""""""""""

Writing your own modules
------------------------

