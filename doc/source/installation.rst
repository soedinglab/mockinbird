Installation
============

Using docker
------------

The easiest way to obtain a working pipeline is by using docker: 
::

        docker pull soedinglab/mockinbird

A host folder with data can be mounted to the ``/data`` directory using the ``-v`` argument.

::

        docker run -v '/path/to/host/folder:/data' soedinglab/mockinbird preprocess \
                input.fastq
                output_folder
                prefix
                config.yaml

All data written to the ``/data`` directory in the docker container will be accessible from the host.

Via bioconda
------------

Users of the `anaconda <https://www.continuum.io/downloads>`_ python distribution can install
mockinbird and all dependencies using ``conda``:

::

        conda install -c bioconda -c conda-forge -c r mockinbird

Note that ``mockinbird`` requires python3 and thus has to be installed in a python3 environment


Manual installation
-------------------

In order to run mockinbird and all shipped modules, following tools and libraries have to be
installed. Please make sure that all command line tools are accessible through your ``PATH``
environment variable.

Python libraries
^^^^^^^^^^^^^^^^

- cython (tested: 0.25.2) :cite:`behnel2010cython`
- scipy (tested: 0.19.0)
- numpy (tested: 1.11.3) :cite:`walt2011numpy`
- pandas (tested: 0.19.2)
- matplotlib (tested: 2.0.0) :cite:`Hunter:2007`
- pysam (tested: 0.9.1.4)
- rpy2 (tested: 2.8.2)
- pyyaml (tested: 3.12)
- jinja2 (tested: 2.9.6)

R libraries
^^^^^^^^^^^

- VGAM (tested 1.0-2)
- fdrtool (tested 1.2.15) :cite:`strimmer2008fdrtool`
- LSD (tested 3.0)

Commandline tools
^^^^^^^^^^^^^^^^^

- `R <https://www.r-project.org/>`_ (tested 3.3.1)
- `samtools <http://samtools.sourceforge.net/>`_ (tested 1.4) :cite:`li2009sequence`
- `umi_tools <https://github.com/CGATOxford/UMI-tools>`_ (tested 0.4.3) :cite:`smith2017umi`
- `STAR <https://github.com/alexdobin/STAR>`_ (tested 2.5.3a) :cite:`dobin2013star`
- `bowtie <http://bowtie-bio.sourceforge.net/index.shtml>`_ (tested 1.2.0) :cite:`langmead2010aligning`
- `xxmotif <https://github.com/soedinglab/xxmotif>`_ (tested 1.6.0) :cite:`luehr2012xxmotif`
- `fastqc <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ (tested 0.11.5)
- `skewer <https://github.com/relipmoc/skewer>`_ (tested version 0.2.2) :cite:`jiang2014skewer`
