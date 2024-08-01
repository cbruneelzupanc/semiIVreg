
<!-- README.md is generated from README.Rmd. Please edit that file -->

# semiIVreg: R package for semi-IV regression

<!-- badges: start -->
<!-- badges: end -->

This package provides an implementation of semi-IV regressions as
described in [Bruneel-Zupanc (2024)](https://www.cbruneel.com/).

<div style="text-align: center;">

<img src="vignettes/images/causal_graph.png" alt="" style="width: 100%; max-width: 400px; height: auto;"/>

</div>

## Installation

You can find the development version of semiIVreg from
[GitHub](https://github.com/cbruneelzupanc/semiIVreg). You can download
it from there and then install it directly from the local source on your
computer:

``` r
# If the package is in a .tar.gz file
install.packages("/path/to/your/package.tar.gz", repos = NULL, type = "source")

# If the package is in a directory
install.packages("/path/to/your/package_directory", repos = NULL, type = "source")
```

Alternatively, you can directly download it from the GitHub repository:

``` r
# install.packages("devtools")
devtools::install_github("cbruneelzupanc/semiIVreg")
```

## Semi-IV Regression

This illustrates what the `semiivreg()`command reports for a semi-IV
regression. By default, it reports the common support plot of the
propensity score and the estimated marginal treatment effects (MTE).

``` r
library(semiIVreg)
## KernSmooth 2.23 loaded
## Copyright M. P. Wand 1997-2009
data(roydata) # load the data from a simulated Roy model

# semi-IV regression
semiiv = semiivreg(y~d|w0|w1, data=roydata, est_method="sieve") 
```

<img src="man/figures/README-mte-1.png" width="100%" />

One can also easily extract a plot for the marginal treatment responses
(MTR):

``` r
semiiv$plot$mtr
```

<img src="man/figures/README-mtr-1.png" width="480" style="display: block; margin: auto;" />

For more details, see the detailed documentation for each function and
the detailed vignettes.
