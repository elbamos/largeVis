context("Edge Matrix")

test_that("Edge Matrix doesn't crash", {
	M <- 10
	data(quakes)
dat <- t(scale(as.matrix(quakes)))
	neighbors <- randomProjectionTreeSearch(dat,
																					K = 20,
																					tree_threshold = 30,
																					max_iter = 1,
																					n_trees = 10,
																					threads = 2,
																					verbose = FALSE)
	expect_silent(edges <- buildEdgeMatrix(dat, neighbors, verbose = FALSE))
	expect_equal(attr(edges, "method"), "euclidean")
})

context("wij")

test_that("wij doesn't crash", {
	data(iris)
	set.seed(1974)
	dat <- as.matrix(iris[, 1:4])
	dat <- scale(dat)
	dupes <- which(duplicated(dat))
	dat <- dat[-dupes, ]
	dat <- t(dat)
	neighbors <- randomProjectionTreeSearch(dat, K = 20, threads = 2)
	edges <- buildEdgeMatrix(dat, neighbors)
	expect_silent(wij <- buildWijMatrix(edges, threads = 2))
})

context("project knns")
data(iris)
set.seed(1974)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
neighbors <- randomProjectionTreeSearch(dat, K = 20, threads = 2)
edges <- buildEdgeMatrix(dat, neighbors)
wij <- buildWijMatrix(edges, threads = 2)

test_that("project knns doesn't crash", {
	coords <- projectKNNs(wij, sgd_batches = 100, verbose = FALSE, threads = 2)
})

test_that("project knns doesn't crash with momentum", {
	expect_silent(coords <- projectKNNs(wij, sgd_batches = 100, verbose = FALSE, momentum = 0.5, threads = 2))
})

test_that("project knns doesn't crash with useDegree", {
	expect_silent(coords <- projectKNNs(wij, sgd_batches = 100, verbose = FALSE, useDegree = TRUE, threads = 2))
})

test_that("project knns doesn't crash with useDegree and momentum", {
	expect_silent(coords <- projectKNNs(wij, sgd_batches = 100, verbose = FALSE, useDegree = TRUE, momentum = 0.5, threads = 2))
})

context("sgd batches")

test_that("sgd batches is linear with E", {
	expect_equal(sgdBatches(1e6, 10), sgdBatches(1e6, 1e4))
})

test_that("sgd batches has a small bump", {
	expect_lt(sgdBatches(5000), sgdBatches(10001))
	expect_lt(sgdBatches(9999) / 1.5, sgdBatches(10001))
})

context("dist")

data(iris)
set.seed(1974)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
do <- dist(t(dat))

test_that("build edge matrix as distance matches dist", {
	neighbors <- randomProjectionTreeSearch(dat, K = ncol(dat) - 1, max_iter = 10, threads = 2)
	edges <- buildEdgeMatrix(dat, neighbors)
	d2 <- as_dist_edgematrix(edges)
	expect_equal(as.matrix(do), as.matrix(d2))
	expect_equal(attr(d2, "method"), "euclidean")
	expect_equal(sum(is.na(as.matrix(d2))), 0)
})

test_that("build edge matrix as distance matches dist with nas", {
	neighbors <- randomProjectionTreeSearch(dat, K = 20, threads = 2)
	edges <- buildEdgeMatrix(dat, neighbors)
	d3 <- as_dist_edgematrix(edges)
	todelete <- !is.na(as.matrix(d3))
	expect_equal(as.matrix(do)[todelete], as.matrix(d3)[todelete])
	expect_equal(attr(d3, "method"), "euclidean")
})
