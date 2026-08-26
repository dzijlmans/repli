
<!-- README.md is generated from README.Rmd. Please edit that file -->

# repli

<!-- badges: start -->

<!-- badges: end -->

**Replicate-based imputation for mass spectrometry-based proteomics
data.**

`repli` implements replicate-based imputation for Missing At Radnom
(MAR) values in proteomics data. The method uses measurements available
in other biological or technical replicates of the same protein and
condition to estimate the mean for imputation, while using a global
standard deviation derived from the standard deviations of all
protein-condition pairs in the dataset.

The package supports both standard data frames and
`SummarizedExperiment` objects for use in workflows based on
[DEP](https://bioconductor.org/packages/DEP/).

## Installation

You can install the development version of `repli` from GitHub with:

``` r
# install.packages("pak")
pak::pak("dzijlmans/repli")
```

More information can be found in the vignette:

``` r
browseVignettes("repli")
#> No vignettes found by browseVignettes("repli")
```
