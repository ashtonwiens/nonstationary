
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nonstationary

<!-- badges: start -->

<!-- badges: end -->

Provides functions for specifying and fitting nonstationary covariance
models for Gaussian processes and Gaussian Markov random fields. See the
vignettes on climate ensemble emulation using a GMRF and a simple
example using a nonstationary generalization of the Matern.

## Installation

You can install the development version of this package from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ashtonwiens/nonstationary")
```

To build the vignettes (which takes about 11 minutes on my laptop)

``` r
  devtools::install_github("ashtonwiens/nonstationary", build_vignettes=TRUE)
```

Then

``` r
vignette('climate-ensemble-emulation', package='nonstationary')
vignette('nonstationary-matern-generalization', package='nonstationary')
```
