.. _ref_normalization:

Normalization
=============
PAR-CLIP data depends on the expression and the transcript lifecycle of the underlying transcripts. The higher the expression of a transcript the more signal is measured in CLIP and PAR-CLIP experiments. Therefore, a normalization procedure to correct for transcript abundance is mandatory to compare binding sites in a quantitative manner. STAMMP simply calculates the ratio of observed T-C mutations over the RNA-seq signal at PAR-CLIP binding site to correct for mRNA expression levels. However, the RNA-seq data for PAR-CLIP data normalization should be generated under the same conditions to avoid expression biases introduced by the PAR-CLIP protocol.

.. warning:: We recommend, that you use the same pre-processing for PAR-CLIP and RNA-seq data.

.. automodule:: stammp.scripts.normalize

Example::
    
    $ stammp-normalize /path/parclipsites.table /path/parclipsites_normalized.table RNAseq1.pileup,...,RNAseqN.pileup chr1,...,chrN -s 0,0
    

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

Set max value
-------------
.. automodule:: stammp.scripts.convert2quantile

PAR-CLIP data often suffers from a few extremly bound outliers. Thus, we recommend to set the maximum occupancy value to the 0.95 quantile::
    
    $ stammp-convert2quantile /path/to/inpufile.table /path/to/outputfile.table -q 0.95
    

This is the final file you can use as input for various post-process analysis steps.



