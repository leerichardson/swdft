<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
========

The swdft package implements the Sliding Window Discrete Fourier Transform (SWDFT) and provides some statistical and graphical tools for analyzing the output. The reference paper is available online at (<https://arxiv.org/abs/1807.07797>).

Installation
------------

``` r
# The swdft package is available on CRAN
install.packages("swdft")

# The development version is available on Github 
# install.packages("devtools")
devtools::install_github("leerichardson/swdft/r/swdft")
```

Usage
-----

The primary functions takes a 1D SWDFT of a vector:

``` r
library(swdft)

# The SWDFT of White-Noise 
x <- rnorm(n = 40)
a <- swdft(x, n = 2^4)
plot_swdft(a)
```

But the package is especially useful time-series w/ periodicities:

``` r
# The SWDFT of a Local Periodic Signal
x_periodic <- local_signal(N=40, A=1, Fr=8, phase=0, S=0, L=40)
a_periodic <- swdft(x, n=2^4)
plot_swdft(a_periodic)
```