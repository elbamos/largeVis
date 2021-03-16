cutoptics <- function(x) {
	clust <- 1
	eps <- x$eps
	minPts <- x$minPts
	clusters <- rep(0, length(x$order))
	for (i in 1:length(x$order)) {
		thisNode <- x$order[i]
		if (is.na(x$predecessor[thisNode]) | x$reachdist[thisNode] > eps) clust <- clust + 1
		clusters[thisNode] <- clust
	}
	badClusters <- which(tabulate(clusters) < minPts)
	clusters[clusters %in% badClusters] <- NA
	remainders <- unique(clusters)
	match(clusters, remainders)
}

context("dbscan-iris")

set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
K <- 147
neighbors <- randomProjectionTreeSearch(dat, K = K,
																				n_trees = 10,  max_iter = 4,
																				verbose = FALSE)
edges <- buildEdgeMatrix(data = dat,
												 neighbors = neighbors,
												 verbose = FALSE)

test_that("dbscan doesn't crash on iris", {
	expect_silent(lv_dbscan(edges = edges, neighbors = neighbors, eps = 1, minPts = 10, verbose = FALSE))
})

load(system.file(package = "largeVis", "testdata/irisdbscan.Rda"))

test_that("dbscan matches iris", {
	dbclusters <- lv_dbscan(edges = edges, neighbors = neighbors, eps = 1, minPts = 10, verbose = FALSE)
	expect_lte(sum(dbclusters$cluster != irisclustering$cluster), 1)
})

test_that("dbscan works with largeVis objects", {
	vis <- largeVis(dat, sgd_batches = 1)
	expect_silent(cl <- lv_dbscan(vis, eps = 1, minPts = 10))
	expect_lte(sum(cl$cluster != irisclustering$cluster), 1)
})

context("optics-iris")

set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
K <- 147
neighbors <- randomProjectionTreeSearch(dat, K = K, tree_threshold = 80, n_trees = 10,  max_iter = 4, verbose = FALSE)
edges <- buildEdgeMatrix(data = dat,
												 neighbors = neighbors,
												 verbose = FALSE)

test_that("optics doesn't crash on iris", {
	expect_silent(lv_optics(edges = edges, neighbors = neighbors, eps = 10, minPts = 10, useQueue = FALSE, verbose = FALSE))
	expect_silent(lv_optics(edges = edges, neighbors = neighbors, eps = 10, minPts = 10, useQueue = TRUE, verbose = FALSE))
})

load(system.file(package = "largeVis", "testdata/irisoptics.Rda"))
opclusters <- lv_optics(edges = edges, neighbors = neighbors, eps = 1, minPts = 10,  useQueue = FALSE, verbose = FALSE)

test_that("optics matches optics core infinities", {
	expect_equal(which(is.infinite(opclusters$coredist)), which(is.infinite(irisoptics$coredist)))
})

test_that("optics matches optics core dist not infinities", {
	expect_equal(opclusters$coredist[!is.infinite(opclusters$coredist)], irisoptics$coredist[!is.infinite(irisoptics$coredist)])
})

test_that("opticis iris cut to dbscan matches dbscan", {
	cl <- cutoptics(opclusters)
	dbclusters <- lv_dbscan(edges = edges, neighbors = neighbors, eps = 1, minPts = 10, verbose = FALSE)
	expect_equal(cl, dbclusters$cluster)
})

test_that("optics works with largeVis objects", {
	skip_on_travis()

	vis <- largeVis(dat, sgd_batches = 1)
	expect_silent(cl <- lv_optics(vis, eps = 1, minPts = 10))
	expect_equal(cl$coredist[!is.infinite(cl$coredist)], irisoptics$coredist[!is.infinite(irisoptics$coredist)])
})

test_that("optics works with dbscan", {
	skip_on_travis()

	vis <- largeVis(dat, sgd_batches = 1)
	expect_silent(cl <- lv_optics(vis, eps = 1, minPts = 10, eps_cl = .4, xi = .05))
	expect_equal(cl$coredist[!is.infinite(cl$coredist)], irisoptics$coredist[!is.infinite(irisoptics$coredist)])
})


context("optics-elki")

test_that("optics output format is correct", {
	skip_on_cran()
	skip_on_travis()
	load(system.file("testdata/opttest.Rda", package = "largeVis"))

	x <- opttest$test_data
	neighbors <- randomProjectionTreeSearch(t(opttest$test_data), K = 399, tree_threshold = 100, max_iter = 10, seed = 1974)
	edges <- buildEdgeMatrix(t(opttest$test_data), neighbors = neighbors)

	eps <- .1
	eps_cl <- .1
	minPts <- 10
	res <- lv_optics(edges, neighbors, eps = eps, useQueue = FALSE,  minPts = minPts)
	expect_identical(length(res$order), nrow(x))
	expect_identical(length(res$reachdist), nrow(x))
	expect_identical(length(res$coredist), nrow(x))
	expect_identical(res$eps, eps)
	expect_identical(res$minPts, minPts)
	expect_equal(res$coredist, opttest$elkiopt$coredist)
	optcut <- cutoptics(res)
	refcut <- cutoptics(opttest$elkiopt)
	expect_equal(optcut, refcut)
})
