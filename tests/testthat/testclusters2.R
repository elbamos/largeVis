context("cluster")
library(dbscan, quietly = TRUE)
set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
K <- 20
neighbors <- randomProjectionTreeSearch(dat, K = K, verbose = FALSE)
edges <- buildEdgeMatrix(data = dat,
                          neighbors = neighbors,
                          verbose = FALSE)
test_that("optics doesn't crash on iris with neighbors and data", {
  expect_silent(largeVis:::optics(neighbors = neighbors, data = dat, eps = 10, minPts = 10, verbose = FALSE))
})

test_that("optics doesn't crash on iris with edges", {
  expect_silent(largeVis:::optics(edges = edges, eps = 10, minPts = 10, verbose = FALSE))
})

test_that("optics doesn't crash on iris with edges and data", {
  expect_silent(largeVis:::optics(edges = edges, data = dat, eps = 10, minPts = 10, verbose = FALSE))
})

test_that("dbscan doesn't crash on iris with edges", {
  expect_silent(largeVis:::dbscan(edges = edges, eps = 10, minPts = 10, verbose = FALSE, partition = FALSE))
})

test_that("dbscan doesn't crash on iris with partitions", {
  expect_silent(clusters <- largeVis:::dbscan(edges = edges, eps = 10, minPts = 10,
                                   verbose = FALSE, partition = TRUE))
})

test_that(paste("LOF is consistent", K), {
	truelof <- dbscan::lof(t(dat), k = K)
	ourlof <- largeVis:::lof(edges)
	expect_lt(sum(truelof - ourlof)^2 / ncol(dat), 0.4)
})

test_that("LOF is consistent 10", {
	edges <- buildEdgeMatrix(data = dat,
													 neighbors = neighbors[1:10,],
													 verbose = FALSE)
	truelof <- dbscan::lof(t(dat), k = 10)
	ourlof <- largeVis:::lof(edges)
	expect_lt(sum(truelof - ourlof)^2 / ncol(dat), 0.4)
})