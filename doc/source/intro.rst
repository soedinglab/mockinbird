What is mockinbird?
###################

Introduction & Motivation
*************************

PAR-CLIP :cite:`hafner2010transcriptome` is a protocol derived from the general CLIP-seq technique and allows predicting protein-RNA binding sites at nucleotide resolution. The resolution is achieved by introducing a base-analogon that is misread by the reverse-transcriptase when crosslinked to a protein.

The complexity in the protocol requires several uncommon steps in the data preparation and quality control, such as analyses of the quality and distribution of the PAR-CLIP specific mutations.

The mockinbird PAR-CLIP pipeline offers modules to performs data processing starting from the raw fastq file down to biological binding analyses. While the predefined modules should make it possible to perform a full analysis, it can easily be extended by defining custom modules.

Last but not least, the pipeline offers a module for predicting PAR-CLIP binding sites by incorporating a measured mock PAR-CLIP experiment. To the best to our knowledge this is the first quantitative model that estimates the binding probability by integrating the coverage and transition properties of the mock experiment. 


How does it work?
*****************

The pipeline is divided into two phases. In each phase a predefined list of modules are executed sequentially. Each phase is controlled by a configuration file in `yaml <https://en.wikipedia.org/wiki/YAML>`_ format. The configuration files contain all information about the experiment design and the chain of modules with their additional arguments.

The input of the **preprocessing phase** is the raw reads in ``fastq`` format. mockinbird ships modules for common analyses and transformations such as quality control, adapter clipping, `UMI <https://en.wikipedia.org/wiki/Unique_molecular_identifier>`_-deduplication, read mapping, binding site prediction, and many more.

Modules of the **postprocessing phase** can be used to obtain analyses of the binding properties derived from the predicted binding sites. mockinbird ships default modules for meta-gene plots, binding heatmaps enriched sequence analysis.


Target users
************

The pipeline should be ideal for experimentalists when troubleshooting experimental PAR-CLIP conditions. It can help both bioinformaticians and experimentalists to check and generate first hypotheses about PAR-CLIP data in a reproducible way. The mockinbird pipeline certainly does not replace a careful bioinformatic analysis, but we hope it save time building up a pipeline from the scratch and help gaining deeper insights into the experimental quality and biological properties of PAR-CLIP'ed factors. 

Requirements
************

The only strict requirement of the mockinbird pipeline is the **raw fastq file** of the PAR-CLIP experiment.

Most downstream analyses in the post-processing phase rely on transforming transitions to factor occupancies. Occupancies can be estimated by sequencing **total RNA** in an RNA-seq experiment. To avoid biases due to stress response of the cells under PAR-CLIP conditions, it is recommended to sequence RNA in a PAR-CLIP setting, but before the immunoprecipitation step.

Optionally, if you want to try out our new mock-based prediction of PAR-CLIP binding sites, you need a **PAR-CLIP experiment under mock conditions**.

Why measuring a PAR-CLIP mock?
******************************

Previously published PAR-CLIP methods predict PAR-CLIP binding regions from read clustering and pinpoint binding sites by comparing the transition rates to a random model. Transition rates in the randommodel are estimated from the transions that are not specific to the PAR-CLIP experiment :cite:`corcoran2011paralyzer`, :cite:`comoglio2015sensitive`, :cite:`chen2014pipe`.

While these methods are well suited for distinguishing sites that arose from sequencing and mapping errors, they can not account for truely crosslinked RNA, which however is not bound to the PAR-CLIP factor of interest.

A recent study has shown that this so-called **background binding** is a major confounder in PAR-CLIP experimnents, despite of the stringent washing steps in the PAR-CLIP protocol that are supposed to select for interactions with the protein of interest :cite:`friedersdorf2014advancing`.

A mock experiment is a PAR-CLIP experiment with an antibody against a protein that is known for not being able to bind to RNA in vivo. As with regular PAR-CLIP experiments, the PAR-CLIP specific transition dominates all other transitions. The transitions observed in the mock experiment can be used to separate factor-specific binding sites from background binding sites.
