.. _ref_binding-site-detection:

Binding Site Detection
======================
The identification of reliable binding sites is a crucial part of PAR-CLIP data analysis as it forms the basis for all subsequent analysis steps. STAMMP distinguishes between experimentally induced and non experimentally induced mutations. The underlying model calculates the probability that an observed number of mutations at a specific genomic position is caused by technical errors rather than from crosslink induced mutations. Therefore, a mixture model of a betabinomial and binomial distribution is fitted to non-T-C mutations in an PAR-CLIP experiment which leads to precise p-values for non-experimentally induced mutations. This negative model is used to distinguish true crosslink sites from false ones. 

:mod:`~stammp.scripts.bsfinder`
-------------------------------
.. automodule:: stammp.scripts.bsfinder
   :members:
