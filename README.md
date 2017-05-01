largeVis
================

[![Travis-CI Build Status](https://travis-ci.org/elbamos/largeVis.svg?branch=master)](https://travis-ci.org/elbamos/largeVis) [![Coverage Status](https://img.shields.io/codecov/c/github/elbamos/largeVis/master.svg)](https://codecov.io/gh/elbamos/largeVis/branch/master) [![https://gitter.im/elbamos/largeVis](https://badges.gitter.im/elbamos/largeVis.svg)](https://gitter.im/elbamos/largeVis?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/elbamos/largeVis?branch=master&svg=true)](https://ci.appveyor.com/project/elbamos/largeVis?branch=master)

This is an implementation of the `largeVis` algorithm described in (<https://arxiv.org/abs/1602.00370>). It also incorporates:

-   A very fast algorithm for estimating k-nearest neighbors, implemented in C++ with `Rcpp` and `OpenMP`. See the [Benchmarks](./benchmarks.md) file for performance details.
-   Efficient implementations of the clustering algorithms:
    -   `HDBSCAN`
    -   `OPTICS`
    -   `DBSCAN`
-   Functions for visualizing manifolds like [this](http://cs.stanford.edu/people/karpathy/cnnembed/).

### News Highlights

-   Version 0.1.10 re-adds clustering, and also adds momentum training to largeVis, as well as a host of other features and improvements.
-   Version 0.1.9.1 has been accepted by CRAN. Much grattitude to Uwe Ligges and Kurt Hornik for their assistance, advice, and patience.

### Some Examples

![MNIST](./README_files/figure-markdown_github/drawmnist-1.png)

![Wiki Words](./README_files/figure-markdown_github/drawwikiwords-1.png)

### Clustering With HDBSCAN

![](README_files/figure-markdown_github/clustering-1.png)

### Visualize Embeddings

![Visualize Embeddings](./README_files/figure-markdown_github/faceImages-1.png)

#### Building Notes

-   **Note on R 3.4**: Before R 3.4, the CRAN binaries were likely to have been compiled without OpenMP, and getting OpenMP to work on Mac OS X was somewhat tricky. This should all have changed (for the better) with R 3.4, which natively using `clang 4.0` by default. Since R 3.4 is new, I'm not able to provide advice, but am interested in hearing of any issues and any workarounds to issues that you may discover.
