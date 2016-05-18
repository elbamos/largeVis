largeVis
================

[![Travis-CI Build Status](https://travis-ci.org/elbamos/largeVis.svg?branch=master)](https://travis-ci.org/elbamos/largeVis) [![Coverage Status](https://img.shields.io/codecov/c/github/elbamos/largeVis/master.svg)](https://codecov.io/github/elbamos/largeVis?branch=master)[![https://gitter.im/elbamos/largeVis](https://badges.gitter.im/elbamos/largeVis.svg)](https://gitter.im/elbamos/largeVis?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

This is an implementation of the `largeVis` algorithm described in (<https://arxiv.org/abs/1602.00370>). It also incorporates code for a very fast algorithm for estimating k-nearest neighbors.

The inner loops for nearest-neighbor search and gradient descent are implemented in C++ using `Rcpp` and `RcppArmadillo`.

#### Project Status & Caveats

-   It works!
-   This project is under heavy development.
-   I am attempting to replicate the paper's results with larger and larger datasets. This takes time because my hardware is not as powerful as the authors'. If you have any to volunteer, please contact me!
-   The algorithm is memory intensive. Processing mnist, memory usage peaked at approximately 8GB. I would appreciate any reports using it with larger datasets.
-   Note that your installation of R must be configured to work with OpenMP. I have had a report that on Federa 22, even small datasets could not be processed because of exceeding the C stack space. If you experience any compilation issues or similar crashes, please create an issue.

#### Examples:

##### MNIST

``` r
load("./mnist.Rda")
dat <- mnist$images
dim(dat) <- c(42000, 28 * 28)
dat <- (dat / 255) - 0.5
coords <- vis(dat, check=FALSE,
                   n_tree = 50, tree_th = 200, 
                   K = 50, alpha = 2, max.iter = 4)
```

![](README_files/figure-markdown_github/drawmnist-1.png)
