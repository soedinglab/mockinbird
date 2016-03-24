.. STAMMP documentation master file, created by
   sphinx-quickstart on Thu Jan 28 10:26:36 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PAR-CLIP data analysis with STAMMP
==================================

STAMMP - a **sta**\tistical **m**\ixture **m**\odel for **P**\AR-CLIP data - offers an analysis pipeline for PAR-CLIP data obtained by single-end sequencing. `Download <http://wwwuser.gwdg.de/~compbiol/parclip_stammp/data/stammp_02_26.tar.gz>`_ STAMMP and follow the descriptions at :ref:`ref_install`. We recommand to try out the :ref:`ref_tutorial` first to get immediately familiar with the usage of STAMMP. You will find more detailed descriptions about each point within the corresonding sections. As a reference you will find an overview about the available analysis plots and their usage in :ref:`ref_plotting`. If you like to extend STAMMP with your own analysis and like to build upon this project you'll find an overview about the datatypes in :ref:`ref_general`.  

In general, the analysis of PAR-CLIP data can be separated into 4 sections

1. :ref:`ref_pre-processing`
2. :ref:`ref_binding-site-detection`
3. :ref:`ref_normalization`
4. Post-Processing (:ref:`ref_plotting`)

which are all covered by various individual modules. However, you do not have to use the complete pipeline for your data. For example, if you have your own pre-processing procedure or if you are only interested in a subset of analysis feel free to skip parts and select only the modules you need. Please make sure, that you have installed all necessary programs and libraries on your machine (see :ref:`ref_install`).

If you use STAMMP for your analysis please cite:
    * `Keine Publikation <http://www.mpibpc.mpg.de/soeding>`_

STAMMP was successfully used in the following publications:
    * `Transcriptome Maps of mRNP Biogenesis Factors Define Pre-mRNA Recognition. Molecular Cell 55 (5), pp. 745-757 (2014) <http://www.cell.com/molecular-cell/abstract/S1097-2765(14)00638-8>`_
    * `Transcriptome Surveillance by Selective Termination of Noncoding RNA Synthesis. Cell 155 (5), pp. 1075-1087 2013 <http://www.cell.com/cell/abstract/S0092-8674(13)01300-7>`_

Contents
========

.. toctree::
   :maxdepth: 1
   :numbered:

   pre-processing.rst
   binding-site-detection.rst
   normalization.rst
   plotting.rst
   tutorial.rst
   general.rst


.. _ref_install:

Installation & Dependencies
===========================

`Download current STAMMP version <http://wwwuser.gwdg.de/~compbiol/parclip_stammp/data/stammp_02_26.tar.gz>`_

STAMMP is still under development, mainly written in Python 3 and R and was developed and tested in a Linux (Ubuntu 14.04) environment. In order to use STAMMP properly please make sure, that you have installed the listed programs as well as the additional libraries/modules:

* `Python <https://www.python.org/>`_  > 3.2
    * `SciPy and Numpy <http://www.scipy.org/>`_
* `R <https://cran.r-project.org/>`_ > 3.1.2
    * `LSD <https://cran.r-project.org/web/packages/LSD/index.html>`_
* `Java JRE <http://www.oracle.com/technetwork/java/javase/downloads/index.html>`_
* `FASTX-toolkit <http://hannonlab.cshl.edu/fastx_toolkit/download.html>`_
* `samtools <http://www.htslib.org/>`_
* `Bowtie <http://bowtie-bio.sourceforge.net/index.shtml>`_
* `XXmotif <http://xxmotif.genzentrum.lmu.de/index.php?id=download>`_
* `FastQC <http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc>`_

Python3, SciPy and R can be easily installed with your package manager (e.g. synaptic). 

Download STAMMP and unpack it into a directory of your choice. Next, open a terminal and go to the directory where you unpacked STAMMP::
    
    $ cd /path/to/stammp
    
You should find a name *setup.py* in that folder. Install STAMMP with::
    
    $ python3 setup.py install
    
or::
    
    $ sudo python3 setup.py install
    

You need root access to install STAMMP. Please ask your administrator if you do not have root access. If the installation was successful you should be able to run STAMMP's scripts. If you like you can check if you see the help page from one of the scripts by typing::
    
    $ stammp-bsfinder --help
    

If you see the help page of stammp-bsfinder you installed the python package correctly. Start with the :ref:`ref_tutorial` to see how your data can be analyzed.

The analysis of NGS data also requires additional data sets like genomic sequences, mapping indices, genomic annotations etc.. You'll find additional information about required data sets and their specific format type within the specific sections. Download our example data () and take the :ref:`ref_tutorial` to see how data sets are used during analysis.

Installing 3rd party programs
-----------------------------

Bowtie and XXmotif offer pre-compiled versions for downloading and FastQC is a java program. Download the pre-compiled versions and unpack them. These programs must be accessible via the commandline with the commands::
    
    $ bowtie
    $ fastqc
    $ XXmotif
    
To do so, set a symlinks to the starting scripts of these programs in the unpacked directories. Bowtie has a python-script named *bowtie* in the directory you just unpacked. You can set a symlink with something similar like this::
    
    $ sudo ln --force -s /path/to/bowtie-1.1.2/bowtie /usr/local/bin/bowtie
    
Afterwards, you need to change the permissions via::
    
    $ sudo chmod 755 /usr/local/bin/bowtie
    
Repeat this step for XXmotif and fastqc::
    
    $ sudo ln --force -s /path/to/XXmotif/XXmotif /usr/local/bin/XXmotif
    $ sudo ln --force -s /path/to/FastQC/fastqc /usr/local/bin/fastqc
    $ sudo chmod 755 /usr/local/bin/XXmotif
    $ sudo chmod 755 /usr/local/bin/fastqc
    

Open a new terminal and check if you can use the programs, e.g::
    
    $ fastqc --version
    FastQC v0.11.4
    $ bowtie --version
    bowtie version 1.1.2
    64-bit
    Built on localhost.localdomain
    Tue Jun 23 13:28:18 EDT 2015
    Compiler: gcc version 4.1.2 20080704 (Red Hat 4.1.2-54)
    Options: -O3 -m64  -Wl,--hash-style=both -DPOPCNT_CAPABILITY  
    Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}
    $ XXmotif -h | head
    
    ====================
    ==  run XXmotif   ==
    ====================
    
    =======================
    == XXmotif version 1.6
    =======================

To install FastX and samtools follow the installation guidelines for these programs.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

