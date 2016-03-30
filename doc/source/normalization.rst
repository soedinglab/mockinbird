.. _ref_normalization:

Normalization
=============
PAR-CLIP data is dependent of transcript expression and the transcript lifecycle. For example. the higher the expression of a transcript the more signal is measured in CLIP and PAR-CLIP experiments. As a consequence, we suggest to normalize PAR-CLIP data to correct for transcript abundance in order to compare binding sites in a quantitative manner. 

To do so, RNA-seq data should be pre-processed with the same scripts that are used for PAR-CLIP data (e.g. the provided pre-processing script). Next, you can provide pileup files from RNA-seq experiments together with PAR-CLIP binding sites obtained by :mod:`~stammp.scripts.bsfinder` to :mod:`~stammp.scripts.normalize` to normalize your PAR-CLIP data. 

:mod:`~stammp.scripts.normalize` simply calculates the ratio of observed T-C mutations over the RNA-seq signal at PAR-CLIP binding site to correct for RNA expression levels. However, the RNA-seq data for PAR-CLIP data normalization should be generated under the same conditions to avoid expression biases introduced by the PAR-CLIP protocol.

If you do not have RNA-seq data for normalization you can replace the normal normalization script with :mod:`~stammp.scripts.normalizeFake`.

After normalization, we suggest to get rid of extreme outliers in your data with :mod:`~stammp.scripts.convert2quantile`.

:mod:`~stammp.scripts.normalize`
--------------------------------
.. automodule:: stammp.scripts.normalize

.. warning:: We recommend, that you use the same pre-processing for PAR-CLIP and RNA-seq data.
.. _ref_normalization_s:

What about **-s**?
^^^^^^^^^^^^^^^^^^
Historically we started with single-end sequencing for normalization which led to some problems regarding usability especially with paired end data. So this might be updated in the future. We highly recommend to check the orientation of your reads in a genome browser before you use your pileups for normalization. 

.. warning:: If you mess up this step your whole post-processing will be wrong, due to wrong calculated protein-RNA affinities. So please make sure that you checked your data!

For example in the image below the data of RNAseq_3 and RNAseq_4 is in correct orientation (compare reads and transcript annotation). On the other hand, RNAseq_1 and RNAseq_2 are wrong. A correct call of the normalization script regarding the *pileups* and *-s* settings would be::
    
    $ stammp-normalize parclip.table parclip_normalized.table RNAseq_1.pileup,RNAseq_2.pileup,RNAseq_3.pileup,RNAseq_4.pileup chr1,...,chrN -s 1,1,0,0
    

.. image:: img/img_rnaseq.png
    :align: center
    :height: 450px
    :alt: alternate text

:mod:`~stammp.scripts.normalizeFake`
------------------------------------
.. automodule:: stammp.scripts.normalizeFake

To get the most information out of your data, RNA-seq normalization is highly recommended. However, you can analyse your data without RNA-seq via :mod:`~stammp.scripts.normalizeFake`. The following plot shows the difference of RNA-seq normalization (left) vs. fake normalization (right).

.. image:: img/img_rnaseq_vs_fake.png
    :align: center
    :height: 450px
    :alt: alternate text

:mod:`~stammp.scripts.convert2quantile`
---------------------------------------
Normally, PAR-CLIP data suffers from having a couple of extreme outliers. Thus, we set the maximum value to a defined quantile (0.95 by default).

.. automodule:: stammp.scripts.convert2quantile

::
    
    $ stammp-convert2quantile /path/to/inpufile.table /path/to/outputfile.table -q 0.95
    

This is the final file you can use as input for various post-process analysis steps.



