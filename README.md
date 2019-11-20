
# scBIRD

<!-- badges: start -->
<!-- badges: end -->

This package contains tools to predict and visualize chromatin accessibility using singe-cell RNA-seq.

## Installation

You can install the released version of scBIRD from [CRAN](https://CRAN.R-project.org) with:

``` r
devtools::install_github("liyueyikate/scBIRD")
```

## Example

Using BIRD to predict chromatin accessibility

``` r
library(scBIRD)
predict("/path/to/expression_matrix/","/path/to/output/file/","/path/to/model/","/path/to/expression_matrix_match/",
locus_model=0,double up_bound=14,match_mode=0, write_flag = 0)
```

