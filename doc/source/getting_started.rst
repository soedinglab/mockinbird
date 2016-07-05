Getting started
===============

Overview
--------

Installation
------------

.. _ref_adapt_configfiles:

Generating and adapting the config files
----------------------------------------

Both the preprocessing and the postprocessing pipeline of STAMMP use a config file to set all parameters.
This is very convenient as config files make it very easy to execute a unified set of scripts that can be easily shared and reproduced on other machines.

For most users the first step is to use :ref:`ref_script_configgen` to generate default config files and adapt them to their needs::

The following command will install config files into `stammp_cfg` in the user's home directory::

        stammp-preprocess ~/stammp_cfg

Having created the default config files, open ``~/stammp_cfg`` in your favorite text editor and set ``adapter5prime`` and ``adapter3prime`` to the adapter sequences used in the experiment.

Next set ``genomeindex`` to the path of a STAR generated genome of the organism of interest and ``genomefasta`` to the fasta file containing the complete genome.

Set ``normalization_pileup`` to a pileup file (created by ``samtools mpileup``) of reads suitable for occupancy normaliation (commonly RNAseq under PAR-CLIP conditions).

Finally, if the adapters contain barcodes, set ``bc_5prime`` and ``bc_3prime`` accordingly.



