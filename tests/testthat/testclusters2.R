context("LOF")

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

test_that("hdbscan finds 3 clusters and outliers in spiral", {
	load(system.file("extdata/spiral.Rda", package = "largeVis"))
	clustering <- hdbscan(spiral, K = 3, minPts = 20, threads = 1)
	expect_equal(length(unique(clustering$clusters, 0)), 3)
})

set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
K <- 20
neighbors <- randomProjectionTreeSearch(dat, K = K,  threads = 2, verbose = FALSE)
edges <- buildEdgeMatrix(data = dat, neighbors = neighbors, verbose = FALSE)


test_that("hdbscan doesn't crash without neighbors", {
	expect_silent(hdbscan(edges, minPts = 10, K = 4, threads = 2, verbose = FALSE))
})

test_that("hdbscan doesn't crash with neighbors", {
	expect_silent(hdbscan(edges, minPts = 10, neighbors = neighbors, K = 4, threads = 2, verbose = FALSE))
})

test_that("hdbscan doesn't crash without 3 neighbors", {
	expect_silent(hdbscan(edges, minPts = 10, K = 3, threads = 2, verbose = FALSE))
})

clustering <- hdbscan(edges, minPts = 10, K = 3,  threads = 2, verbose = FALSE)

test_that("hdbscan is correct", {
	expect_equal(length(unique(clustering$clusters, 0)), 3)
})

test_that("gplot isn't broken", {
	expect_silent(plt <- gplot(clustering, matrix(rnorm(ncol(neighbors) * 2), ncol = 2)))
})

test_that("hdbscan doesn't crash on glass edges", {
	load(system.file("extdata/glassEdges.Rda", package = "largeVis"))
	clustering <- hdbscan(edges, threads = 2, verbose = FALSE)
	expect_equal(length(unique(clustering$clusters)), 3)
})

test_that("hdbscan doesn't crash on big bad edges", {
	skip("skipping big bad edges test because the data is too big for cran")
	load(system.file("extdata/kddneighbors.Rda", package = "largeVis"))
	load(system.file("extdata/kddedges.Rda", package = "largeVis"))
	expect_silent(clusters <- hdbscan(edges, neighbors = neighbors, threads = 2, verbose = FALSE))
})

test_that("hdbscan doesn't crash on big bad edges", {
	skip_on_cran()
	skip("skipping long test")
	load(system.file("extdata/badedges.Rda", package = "largeVis"))
	expect_silent(clustering <- hdbscan(badedges, threads = 2, verbose = FALSE))
})
