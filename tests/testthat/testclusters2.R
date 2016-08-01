context("cluster")
library(largeVis)
set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
neighbors <- randomProjectionTreeSearch(dat, K = 20, verbose = FALSE)
edges <- buildEdgeMatrix(data = dat,
                          neighbors = neighbors,
                          verbose = FALSE)
test_that("optics doesn't crash on iris with neighbors and data", {
  expect_silent(optics(neighbors = neighbors, data = dat, eps = 10, minPts = 10, verbose = FALSE))
})

test_that("optics doesn't crash on iris with edges", {
  expect_silent(optics(edges = edges, eps = 10, minPts = 10, verbose = FALSE))
})

test_that("optics doesn't crash on iris with edges and data", {
  expect_silent(optics(edges = edges, data = dat, eps = 10, minPts = 10, verbose = FALSE))
})

test_that("dbscan doesn't crash on iris with edges", {
  expect_silent(dbscan(edges = edges, eps = 10, minPts = 10, verbose = FALSE, partition = FALSE))
})

test_that("dbscan doesn't crash on iris with partitions", {
  expect_silent(clusters <- dbscan(edges = edges, eps = 10, minPts = 10,
                                   verbose = FALSE, partition = TRUE))
})
