Postprocessing
==============

Getting started
---------------

The postprocessing script is a convenient way to automatically run a batch of downstream analyses.
Every downstream analysis consists of one or more modules which are defined in a postprocessing config file.

Each section in the config file defines one module::

        [Centerplot]                    <-- module name in square brackets
        skip = yes                      <-- skip module?
        type = center_plot              <-- module type
        ...

A module is defined by its name in squared brackets. The module name can be chosen arbitarily.
Optionally the ``skip`` attribute can be set in case the module should *not* be executed.
The ``type`` attribute defines the kind of module. All other arguments in the module body depend on the module type.

.. note::

        It is possible to run the same module type with different parameters in one downstream analysis. Modules define an ``output_prefix`` option that can be set to distinguish the created output files.

For a detailed description of the module types and their type keys see :ref:`ref_module_reference`.


Command line script
-------------------

.. argparse::
  :module: stammp.scripts.postprocess
  :func: create_parser
  :prog: stammp-postprocess


.. _ref_module_reference:

Module Reference
----------------

Centerplot
~~~~~~~~~~

.. image:: imgs/img_plotCenterBoth.png
  :align: center
  :height: 700px

Module type key: ``centerplot``

This module wraps the :ref:`ref_script_center_plot` command line script.

====================  =======  =============================================
Option                Default  Description
====================  =======  =============================================
gff_file              None     gff file with annotations
output_prefix         None     prefix of the output files
downstream_bp         1000     nt downstream of the end
upstream_bp           1000     nt upstream of the start
gene_bp               750      nt in the annotation body
min_trscr_size_bp     0        minimum size for transcripts
max_trscr_size_bp     5000     maximum size for transcripts
smoothing_window      20       size of windows for smoothing
labelCenterA          None     label of the start site
labelCenterB          None     label of the end site
labelBody             None     label of the annotation body
remove_tmp_files      yes      remove temporary files
====================  =======  =============================================

KmerCountPlot
~~~~~~~~~~~~~

.. image:: imgs/img_kmerPerPosition.png
  :align: center

Module type key: ``kmer_count_plot``

This module wraps the :ref:`ref_script_kmerperpos` command line script.

====================  =======  =============================================
Option                Default  Description
====================  =======  =============================================
genome_fasta          None     path to genome fasta file
output_prefix         None     prefix of the output files
kmer_k                3        kmer length
first_index           0        first index of PAR-CLIP sites
last_index            1500     last index of PAR-CLIP sites
width                 50       number of nt +/- of the crosslink site
sort_key              occ      key that is used for PAR-CLIP site ordering
gff_exclude_path               skip sites that overlap annotations in this  gff file
gff_padding           20       number of nt added to start/stop indices of the gff annotations
remove_tmp_files      yes      remove temporary files
====================  =======  =============================================


KmerLogOddPlot
~~~~~~~~~~~~~~

.. image:: imgs/img_plotKmerLogOdds.png
  :align: center
  :height: 700px

Module type key: ``kmer_logodd_plot``

This module wraps the :ref:`ref_script_kmerlogodds` command line script.

====================  =======  =============================================
Option                Default  Description
====================  =======  =============================================
genome_fasta          None     path to genome fasta file
output_prefix         None     prefix of the output files
kmer_k                3        kmer length
sort_key              occ      key that is used for PAR-CLIP site ordering
gff_exclude_path               skip sites that overlap annotations in this  gff file
use_quantiles         yes      use quantiles for binarization instead of fixed bin
negative_set_gff      None     gff file for negative sequence set sampling
n_negative_seqs       20000    number of negative sequences sampled
====================  =======  =============================================

XXmotif
~~~~~~~

Module type key: ``xxmotif``

This module wraps the :ref:`ref_script_xxmotif` command line script.

====================  =======  =============================================
Option                Default  Description
====================  =======  =============================================
genome_fasta          None     path to genome fasta file
output_prefix         None     prefix of the output files
n_negative_seqs       20000    number of negative sequences sampled
negative_set_gff      None     gff file for negative sequence set sampling
plot_top_n_pwm        3        plot PWMs of the top n motifs
first_index           0        first index of PAR-CLIP sites
last_index            1500     last index of PAR-CLIP sites
width                 12       number of nt +/- of the crosslink site
sort_key              occ      key that is used for PAR-CLIP site ordering
gff_exclude_path               skip sites that overlap annotations in this  gff file
gff_padding           20       number of nt added to start/stop indices of the gff annotations
remove_tmp_files      yes      remove temporary files
====================  =======  =============================================


Transition Frequency Plot
~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: imgs/img_nuc_probabilities.png
  :align: center

Module type key: ``tr_freq_plot``

This module wraps the :ref:`ref_script_nucprob` command line script.

====================  =======  =============================================
Option                Default  Description
====================  =======  =============================================
output_prefix         None     prefix of the output files
min_cov               5        minimum coverage
y_axis_limit          0        y-axis limit
remove_tmp_files      yes      remove temporary files
====================  =======  =============================================

Heatmap Plot
~~~~~~~~~~~~

.. image:: imgs/img_pub1_heatmap_sense.png
  :align: center


.. image:: imgs/img_pub1_heatmap_asense.png
  :align: center

Module type key: ``heatmap_plot``

This module wraps the :ref:`ref_script_heatmap` command line script.

====================  =======  =============================================
Option                Default  Description
====================  =======  =============================================
gff_file              None     gff file with annotations
output_prefix         None     prefix of the output files
downstream_bp         1000     nt downstream of the end
upstream_bp           1000     nt upstream of the start
min_trscr_size_bp     0        minimum size for transcripts
max_trscr_size_bp     5000     maximum size for transcripts
xbins                 500      number of bins in x direction
ybins                 500      number of bins in y direction
x_pixels              500      width of final plot in px
y_pixels              500      height of final plot in px
remove_tmp_files      yes      remove temporary files
====================  =======  =============================================

Heatmap Plot (small)
~~~~~~~~~~~~~~~~~~~~

.. image:: imgs/img_pub1_heatmap_sense_small.png
  :align: center


.. image:: imgs/img_pub1_heatmap_asense_small.png
  :align: center

Module type key: ``heatmap_small_plot``

This module wraps the :ref:`ref_script_heatmap_small` command line script.

====================  =======  =============================================
Option                Default  Description
====================  =======  =============================================
gff_file              None     gff file with annotations
output_prefix         None     prefix of the output files
downstream_bp         1000     nt downstream of the end
upstream_bp           1000     nt upstream of the start
min_trscr_size_bp     0        minimum size for transcripts
max_trscr_size_bp     5000     maximum size for transcripts
xbins                 500      number of bins in x direction
ybins                 500      number of bins in y direction
x_pixels              500      width of final plot in px
y_pixels              500      height of final plot in px
remove_tmp_files      yes      remove temporary files
====================  =======  =============================================

Filter sites
~~~~~~~~~~~~

Module type key: ``gff_filter``

This module wraps the :ref:`ref_script_filter_sites` command line script.

====================  =======  =============================================
Option                Default  Description
====================  =======  =============================================
filter_gff            None     gff file used for filtering
file_postfix          filt     suffix to be appended to the origial file name
padding_bp            10       bp to extend the annotation
features                       comma separated list of gff features
====================  =======  =============================================
