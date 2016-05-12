largeVis
================

This is an implementation of the `largeVis` algorithm described in (<https://arxiv.org/abs/1602.00370>). It also incorporates code for a very fast algorithm for estimating k-nearest neighbors.

The inner loops for nearest-neighbor search and gradient descent are implemented in C++ using `Rcpp` and `RcppArmadillo`. (If you get an error that a `NULL value passed as symbol address` this relates to `Rcpp` and please open an issue here.)

This has been tested and confirmed to work in many circumstances. More extensive documentation and examples are being prepared.

Please note that this package is under development (the paper is only two weeks old) so it is likely that implementation bugs will be found and changes made to the api.

Some notes:

-   I see very occasional seg faults running in Rstudio (never in console) that appear to be related to OpenMP but I have not yet been able to track down.
-   Several phases are implemented with `parallel::mclapply`. The number of parallel threads can be adjusted with `options(mc.cores = n)`.

### Examples:

``` r
library(largeVis)
library(ggplot2)
data(iris)
dat <- as.matrix(iris[,1:4])
dat <- scale(dat)
dupes = which(duplicated(dat))
dat <- dat[-dupes,]
coords <- largeVis(dat, pca.first = F, 
                   max.iter = 5, sgd.batches = 2000000, 
                   gamma = 7, K = 40, M = 5, rho = 1,min.rho = 0, verbose = FALSE)
coords <- data.frame(coords$coords)
colnames(coords) <- c("X", "Y")
coords$algo <- "largeVis"

library(Rtsne)
more <- Rtsne(dat, initial_dims = 4, pca = F)$Y
more <- data.frame(X = more[,1], Y = more[,2], algo = "B-H t-SNE")
coords <- rbind(coords, more)
coords$Species <- rep(iris$Species[-dupes])

ggplot(coords, aes(x = X, y = Y, color = Species)) + 
  geom_point(size = 0.5) + 
  facet_grid(~ algo) + 
  ggtitle("The Iris Dataset largeVis vs. t-SNE")
```

![](README_files/figure-markdown_github/iris-1.png)<!-- -->
