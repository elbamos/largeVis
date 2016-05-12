largeVis
================

This is an implementation of the `largeVis` algorithm described in (<https://arxiv.org/abs/1602.00370>). It also incorporates code for a very fast algorithm for estimating k-nearest neighbors.

The inner loops for nearest-neighbor search and gradient descent are implemented in C++ using `Rcpp` and `RcppArmadillo`. (If you get an error that a `NULL value passed as symbol address` this relates to `Rcpp` and please open an issue here.)

This has been tested and confirmed to work in many circumstances. More extensive documentation and examples are being prepared.

Please note that this package is under development (the paper is only two weeks old) so it is likely that implementation bugs will be found and changes made to the api.

#### Project Status

Most of the computationally intensive code has been re-written in C++, and I am experimenting with larger datasets and the datasets from the original paper for testing and to attempt to reproduce the results.

#### Some Notes:

-   I see very occasional seg faults running in Rstudio (never in console) that appear to be related to OpenMP but I have not yet been able to track down.
-   Several phases are implemented with `parallel::mclapply`. The number of parallel threads can be adjusted with `options(mc.cores = n)`.

#### Examples:

``` r
library(largeVis)
data(iris)
dat <- as.matrix(iris[,1:4])
dat <- scale(dat)
dupes = which(duplicated(dat))
dat <- dat[-dupes,]
coords <- largeVis(dat, pca.first = F, 
                   max.iter = 10, sgd.batches = 2000000, 
                   K = 20, verbose = FALSE, rho = 1, gamma = 7, M = 10)
coords <- data.frame(coords$coords)
colnames(coords) <- c("X", "Y")
coords$Species <- iris$Species[-dupes]
coords$algo <- "largeVis"
```

``` r
library(Rtsne)
more <- Rtsne(dat, initial_dims = 4, pca = F)$Y
more <- data.frame(X = more[,1], Y = more[,2], algo = "B-H t-SNE")
more$Species <- iris$Species[-dupes]
coords <- rbind(coords, more)
coords$algo <- factor(coords$algo)
save(coords, file = "iriscoords.Rda")
```

``` r
ggplot(coords, aes(x = X, y = Y, color = Species)) + 
  geom_point(size = 0.5) + 
  facet_grid(~ algo) + 
  ggtitle("The Iris Dataset largeVis vs. t-SNE")
```

![](README_files/figure-markdown_github/showiris-1.png)<!-- -->
