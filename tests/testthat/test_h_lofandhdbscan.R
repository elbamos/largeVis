context("LOF")

set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
K <- 80
neighbors <- randomProjectionTreeSearch(dat, K = K,  threads = 2, verbose = FALSE)
edges <- buildEdgeMatrix(data = dat,
												 neighbors = neighbors,
												 verbose = FALSE)

test_that(paste("LOF is consistent", 20), {
	load(system.file("testdata/truelof20.Rda", package = "largeVis"))
	ourlof <- lof(edges)
	expect_lt(sum(truelof20 - ourlof)^2 / ncol(dat), 0.4)
})

test_that("LOF is consistent 10", {
	edges <- buildEdgeMatrix(data = dat,
													 neighbors = neighbors[1:10,],
													 verbose = FALSE)
	load(system.file("testdata/truelof10.Rda", package = "largeVis"))
	ourlof <- lof(edges)
	expect_lt(sum(truelof10 - ourlof)^2 / ncol(dat), 0.4)
})

context("hdbscan")

test_that("hdbscan finds 3 clusters and outliers in spiral", {
	load(system.file("testdata/spiral.Rda", package = "largeVis"))
	clustering <- hdbscan(spiral, K = 3, minPts = 20, threads = 1)
	expect_equal(length(unique(clustering$clusters)), 3)
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
	expect_silent(hdbscan(edges, neighbors = neighbors, minPts = 10, K = 3, threads = 2, verbose = FALSE))
})

test_that("hdbscan accepts a largevis object", {
	vis <- largeVis(dat, sgd_batches = 1, threads = 2)
	expect_silent(clusters <- hdbscan(vis, minPts = 10, K = 3, threads = 2, verbose = FALSE))
	expect_equal(length(unique(clusters$clusters, 0)), 3)
})

clustering <- hdbscan(edges, minPts = 10, K = 3,  threads = 2, verbose = FALSE)

test_that("hdbscan is correct", {
	expect_equal(length(unique(clustering$clusters, 0)), 3)
})

test_that("hdbscan doesn't crash on glass edges", {
	load(system.file("testdata/glassEdges.Rda", package = "largeVis"))
	clustering <- hdbscan(edges, threads = 2, verbose = FALSE)
	expect_equal(length(unique(clustering$clusters)), 3)
})

test_that("hdbscan doesn't crash on big bad edges", {
	skip("skipping big bad edges test because the data is too big for cran")
	load(system.file("testdata/kddneighbors.Rda", package = "largeVis"))
	load(system.file("testdata/kddedges.Rda", package = "largeVis"))
	expect_silent(clusters <- hdbscan(edges, neighbors = neighbors, threads = 2, verbose = FALSE))
})

test_that("hdbscan doesn't crash on big bad edges", {
	skip_on_cran()
	skip("skipping long test")
	load(system.file("testdata/badedges.Rda", package = "largeVis"))
	expect_silent(clustering <- hdbscan(badedges, threads = 2, verbose = FALSE))
})

context("as.dendrogram")

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

test_that("as.dendrogram succeeds on iris4", {
	hdobj <- hdbscan(edges, neighbors = neighbors, minPts = 10, K = 4, threads = 2, verbose = FALSE)
	dend <- as_dendrogram_hdbscan(hdobj)
	expect_equal(length(dend[[1]]), sum(hdobj$hierarchy$nodemembership == 1) + sum(hdobj$hierarchy$parent == 1) - 1)
	expect_equal(sum(is.null(dend)), 0)
	expect_equal(class(dend), "dendrogram")
	expect_equal(nobs(dend), ncol(dat))
}	)

test_that("as.dendrogram succeeds on iris3", {
	hdobj <- hdbscan(edges, minPts = 10, K = 3, threads = 2, verbose = FALSE)
	dend <- as_dendrogram_hdbscan(hdobj)
	expect_equal(length(dend), sum(hdobj$hierarchy$nodemembership == 1) + sum(hdobj$hierarchy$parent == 1) - 1)
	expect_equal(sum(is.null(dend)), 0)
	expect_equal(class(dend), "dendrogram")
	expect_equal(nobs(dend), ncol(dat))
}	)

context("gplot")

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
clustering <- hdbscan(edges, neighbors, minPts = 10, K = 4,  threads = 2, verbose = FALSE)

test_that("gplot isn't broken", {
	expect_silent(plt <- gplot(clustering, t(dat)))
})

test_that("gplot isn't broken with text", {
	expect_silent(plt <- gplot(clustering, t(dat), text = TRUE))
})