largeVis
================

This is an implementation of the `largeVis` algorithm described in (<https://arxiv.org/abs/1602.00370>). It also incorporates code for a very fast algorithm for estimating k-nearest neighbors.

The inner loops for nearest-neighbor search and gradient descent are implemented in C++ using `Rcpp` and `RcppArmadillo`. (If you get an error that a `NULL value passed as symbol address` this relates to `Rcpp` and please open an issue here.)

This has been tested and confirmed to work in many circumstances. More extensive documentation and examples are being prepared.

Please note that this package is under development (the paper is only two weeks old) so it is likely that implementation bugs will be found and changes made to the api.

Some notes:

-   There may be a bug in one of the gradients.
-   This implementation uses C++ implementations of the neighbor-exploration and sgd phrases. While the implementations **should** be using OpenMP, I have not been able to determine that they do.
-   The random partition trees and sigma-estimation phases are implemented with `mclapply` from the `parallel` package. The number of cores that will be used may be set with `options(mc.cores = n)`

Examples:
---------

``` r
library(largeVis)
library(ggplot2)
data(iris)
dat <- as.matrix(iris[,1:4])
coords <- largeVis(dat, pca.first = F, 
                   max.iter = 5, sgd.batches = 2000000, 
                   gamma = 7, K = 40, M = 5, rho = 1,min.rho = 0, verbose = FALSE)
```

    ## Estimating sigmas

``` r
coords <- data.frame(coords$coords)
colnames(coords) <- c("X", "Y")
coords$Species <- iris$Species
ggplot(coords, aes(x = X, y = Y, color = Species)) + geom_point(size = 0.5)
```

![](README_files/figure-markdown_github/iris-1.png)<!-- -->
