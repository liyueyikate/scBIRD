scBIRD
======

<!-- badges: start -->
<!-- badges: end -->
This package provides a graphical user interface to predict and
visualize chromatin accessibility using singe-cell RNA-seq. scBIRD is
also available online at

Installation
------------

scBIRD can be installed from github:

    devtools::install_github("liyueyikate/scBIRD")

Example
-------

To launch the graphical user interface, run the following command in R:

    library(scBIRD)
    scBIRDui()

scBIRD User Manual
------------------

scBIRD user manual is available on Github:

BIRD prediction models
----------------------

Four prebuilt models by Weiqiang Zhou for BIRD are available:

1.  RNA-seq model, current release (trained with 167 ENCODE samples):
    <a href="https://github.com/WeiqiangZhou/BIRD-model/releases/download/v1.3/human_hg19_model.bin.zip" class="uri">https://github.com/WeiqiangZhou/BIRD-model/releases/download/v1.3/human_hg19_model.bin.zip</a>

2.  RNA-seq model, previous release (trained with 70 Epigenome Roadmap
    samples):
    <a href="https://github.com/WeiqiangZhou/BIRD-model/releases/download/v1.2/RNAseq_model_file.bin.zip" class="uri">https://github.com/WeiqiangZhou/BIRD-model/releases/download/v1.2/RNAseq_model_file.bin.zip</a>

3.  RNA-seq model for 2 million loci, previous release (trained with 70
    Epigenome Roadmap samples):
    <a href="https://github.com/WeiqiangZhou/BIRD-model/releases/download/v1.0/RNAseq_model_file_2M.bin.zip" class="uri">https://github.com/WeiqiangZhou/BIRD-model/releases/download/v1.0/RNAseq_model_file_2M.bin.zip</a>

4.  Exon Array model:
    <a href="https://github.com/WeiqiangZhou/BIRD-model/releases/download/v1.1/Exonarray_model_file.bin.zip" class="uri">https://github.com/WeiqiangZhou/BIRD-model/releases/download/v1.1/Exonarray_model_file.bin.zip</a>

Contact
-------

Author: Kate(Yueyi) Li, Weiqiang Zhou, Hongkai Ji

Maintainer: Yueyi Li
(<a href="mailto:kateliyueyi2018@gmail.com" class="email">kateliyueyi2018@gmail.com</a>)
