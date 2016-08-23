context("cluster")

set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
K <- 20
neighbors <- randomProjectionTreeSearch(dat, K = K,  threads = 2, verbose = FALSE)
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

test_that(paste("LOF is consistent", 20), {
	load(system.file("extdata/truelof20.Rda", package = "largeVis"))
	ourlof <- largeVis:::lof(edges)
	expect_lt(sum(truelof20 - ourlof)^2 / ncol(dat), 0.4)
})

test_that("LOF is consistent 10", {
	edges <- buildEdgeMatrix(data = dat,
													 neighbors = neighbors[1:10,],
													 verbose = FALSE)
	load(system.file("extdata/truelof10.Rda", package = "largeVis"))
	ourlof <- largeVis:::lof(edges)
	expect_lt(sum(truelof10 - ourlof)^2 / ncol(dat), 0.4)
})

context("hdbscan")

test_that("hdbscan doesn't crash", {
  expect_silent(hdbscan(edges, minPts = 20, K = 3, threads = 2, verbose = FALSE))
})

test_that("hdbscan doesn't crash with neighbors", {
  expect_silent(hdbscan(edges, minPts = 20, neighbors = neighbors, K = 3, threads = 2,  FALSE))
})

test_that("hdbscan is correct", {
  clustering <- hdbscan(edges, minPts = 10, K = 3,  threads = 2, verbose = FALSE)
  expect_equal(length(unique(clustering$clusters)), 2)
})

test_that("hdbscan is correct with neighbors", {
  clustering <- hdbscan(edges, neighbors = neighbors, minPts = 10, K = 3,  threads = 2, FALSE)
  expect_equal(length(unique(clustering$clusters)), 2)
})

test_that("hdbscan doesn't crash on glass edges", {
	load(system.file("extdata/glassEdges.Rda", package = "largeVis"))
	expect_silent(clustering <- hdbscan(edges))
})