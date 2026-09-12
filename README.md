# gdam -- Fast Robust Additive Models using Gamma Divergence

<!-- badges: start -->

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12194563.svg)](https://doi.org/10.5281/zenodo.12194563) -->

<!-- badges: end -->

`gdam` is an `R` package associated with the article "Fast robust additive models using gamma divergence" by [Hui](https://francishui.netlify.and Computing](https://doi.org/10.1007/s11222-026-10960-3). The package offers a computationally efficient method, at least at the time of development, for robust additive modeling based on the gamma divergence. The method assumes an identity link and normally distributed errors, and applies the gamma divergence to the resulting restricted maximum likelihood function so as to obtain a loss function that is less sensitive to outlying responses. Because `gdam` uses a model fitted via the [mgcv](10.32614/CRAN.package.mgcv) package as the starting point for optimization, and leverages many of its existing techniques for estimation and inference, users can take advantage of many of the smoothing options available in the `mgcv` package for constructing additive models.

# Installation

Currently, `gdam` is available and can be installed from Github with the help of `pak` package using:

```         
pak::pkg_install("fhui28/gdam")
```

A `CRAN` package is currently in the works, and (hopefully!) will be available in the near future.

<!-- Alternatively, or if the above does not work, you may download a (hopefully!) stable release of `COQUE` by choosing the latest release on the right hand side of this Github webpage. -->

# Getting started

Users are recommended to either:

1.  Check out the examples in the help files for the `gdam` functions in the package.

2.  Examine the `runsims.R` scripts inside the folders `manuscript/simulations/setting1` and `manuscript/simulations/setting2` and `manuscript/simulations/setting3`, which are template `R` scripts that can be adapted to run simulation settings I to III, respectively, in the associated manuscript. Please note that, given they are designed to reproduce simulation studies, then the `R` scripts also are not necessarily designed to be as user-friendly, and also implement a number of other robust and non-robust additive modeling approaches. Users are recommended to carefully read through the scripts, and the corresponding setting in the associated manuscript, before running them.

    - Note there are also `runsims.R` scripts inside the folders `manuscript/simulations/setting1a` and `manuscript/simulations/setting2a` and `manuscript/simulations/setting3a`, which are template `R` scripts that can be adapted to run supplementary simulation settings Ia to IIIa, respectively, in the associated manuscript, involving heavy tailed distributions.

# If you find any bugs and issues...

If you find something that looks like a bug/issue, please use Github issues and post it up there. As much as possible, please include in the issue:

1.  A description of the bug/issue;
2.  Paste-able code along with some comments that reproduces the problem e.g., using the [reprex](https://CRAN.R-project.org/package=reprex) package. If you also have an idea of how to fix the problem, then that is also much appreciated.
3.  Required data files etc...

Alternatively, please contact the corresponding author at [fhui28\@gmail.com](mailto:fhui28@gmail.com)
