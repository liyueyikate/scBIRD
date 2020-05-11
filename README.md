scBIRD
======

<!-- badges: start -->
<!-- badges: end -->
This package provides a graphical user interface to predict and
visualize chromatin accessibility using singe-cell RNA-seq data.

Installation
------------

scBIRD can be installed from github:

    devtools::install_github("liyueyikate/scBIRD")

To launch the graphical user interface, run the following command in R:

    library(scBIRD)
    scBIRDui()

scBIRD User Manual
------------------

scBIRD [user
manual](https://github.com/liyueyikate/scBIRD/blob/master/manual.pdf) is
available on Github. To get a quick overview of scBIRD, check our
[YouTube](https://youtu.be/wA4WAWnijIQ) video.

BIRD prediction models
----------------------

Four prebuilt models for BIRD are available:

-   [RNA-seq
    model](https://github.com/WeiqiangZhou/BIRD-model/releases/download/v1.3/human_hg19_model.bin.zip),
    current release (trained with 167 ENCODE samples).

-   [RNA-seq
    model](https://github.com/WeiqiangZhou/BIRD-model/releases/download/v1.2/RNAseq_model_file.bin.zip),
    previous release (trained with 70 Epigenome Roadmap samples).

-   [RNA-seq
    model](https://github.com/WeiqiangZhou/BIRD-model/releases/download/v1.0/RNAseq_model_file_2M.bin.zip)
    for 2 million loci, previous release (trained with 70 Epigenome
    Roadmap samples)

-   [Exon Array
    model](https://github.com/WeiqiangZhou/BIRD-model/releases/download/v1.1/Exonarray_model_file.bin.zip)

Contact
-------

Author: Kate(Yueyi) Li, Weiqiang Zhou, Runzhe Li, Hongkai Ji

Maintainer: Kate(Yueyi) Li
(<a href="mailto:kateliyueyi2018@gmail.com" class="email">kateliyueyi2018@gmail.com</a>)
