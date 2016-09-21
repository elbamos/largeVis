largeVis
================

[![Travis-CI Build Status](https://travis-ci.org/elbamos/largeVis.svg?branch=master)](https://travis-ci.org/elbamos/largeVis) [![Coverage Status](https://img.shields.io/codecov/c/github/elbamos/largeVis/master.svg)](https://codecov.io/gh/elbamos/largeVis/branch/master) [![https://gitter.im/elbamos/largeVis](https://badges.gitter.im/elbamos/largeVis.svg)](https://gitter.im/elbamos/largeVis?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/elbamos/largeVis?branch=master&svg=true)](https://ci.appveyor.com/project/elbamos/largeVis?branch=master)

This is an implementation of the `largeVis` algorithm described in (<https://arxiv.org/abs/1602.00370>). It also incorporates:

-   A very fast algorithm for estimating k-nearest neighbors, implemented in C++ with `Rcpp` and `OpenMP`.
-   An efficient implementation of the `HDBSCAN` algorithm for unsupervised clustering, which leverages the nearest neighbor data genereated by `largeVis`.
-   Functions for visualizing manifolds like [this](http://cs.stanford.edu/people/karpathy/cnnembed/).

-   Version 0.1.9.1 has been accepted by CRAN. Much grattitude to Uwe Ligges and Kurt Hornik for their assistance, advice, and patience.
-   I am working on restoring the OPTICS and DBSCAN implementations that were removed from version 0.1.9 for CRAN submission. If you would like to experiment with these, check branch `features/opticsanddbscan`.
-   Anyone watching closely will have noticed that the git branches had become something of a mess. I have moved to git flow in an attempt to restore some sort of sanity. There may still be some inconsistencies in the git over the next few days. I am working on it.

#### Building Notes

-   The CRAN binaries are likely to have been compiled with 32-bit ARMA objects, and without OpenMP on OS X.
-   To get 64-bit ARMA objects, add to `~.R/Makevars`:

        R_XTRA_CXXFLAGS = -DARMA_64BIT_WORD

    and recompile from source. (`devtools::install_github("elbamos/largeVis")` will work).

-   Getting OpenMP support on OS X is a bit tricky. These instructions assume that your OS X installation is already setup to compile with xcode. Also, these directions are not precise, and you may need to fiddle.
    -   Use `homebrew` to install `llvm` version 3.8 or greater.
    -   `brew link --force llvm`. Note that this will make `clang` version 3.8 your default system compiler. This will probably improve your development experience overall, but might conceivably cause issues with xcode.
    -   Add the following to `~.R/Makevars`:

            SHLIB_OPENMP_CFLAGS = -fopenmp

    -   If you use Rstudio (which you should), and are running OS X Yosemite or later, add the following to `~/.Renviron`:

            PATH=/usr/local/bin:${PATH}

        This is necessary so that R will see the `llvm` compiler from within Rstudio.
    -   Recompile from source, as above.
