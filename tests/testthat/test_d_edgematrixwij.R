context("Edge Matrix")

test_that("Edge Matrix doesn't crash", {
	data(quakes)
	dat <- t(scale(as.matrix(quakes)))
	neighbors <- randomProjectionTreeSearch(dat,
																					K = 5,
																					max_iter = 4,
																					n_trees = 20,
																					verbose = FALSE)
	expect_silent(edges <- buildEdgeMatrix(dat, neighbors, verbose = FALSE))
	expect_equal(attr(edges, "Metric"), "euclidean")
	expect_equal(class(edges), "edgematrix")
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
	neighbors <- randomProjectionTreeSearch(dat, K = 20)
	edges <- buildEdgeMatrix(dat, neighbors)
	expect_silent(wij <- buildWijMatrix(edges))
})

context("project knns")
data(iris)
set.seed(1974)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
neighbors <- randomProjectionTreeSearch(dat, K = 20)
edges <- buildEdgeMatrix(dat, neighbors)
wij <- buildWijMatrix(edges)

test_that("project knns doesn't crash", {
	coords <- projectKNNs(wij, sgd_batches = 100, verbose = FALSE)
})

test_that("project knns doesn't crash with momentum", {
	expect_silent(coords <- projectKNNs(wij, sgd_batches = 100, verbose = FALSE, momentum = 0.5))
})

test_that("project knns doesn't crash with useDegree", {
	expect_silent(coords <- projectKNNs(wij, sgd_batches = 100, verbose = FALSE, useDegree = TRUE))
})

test_that("project knns doesn't crash with useDegree and momentum", {
	expect_silent(coords <- projectKNNs(wij, sgd_batches = 100, verbose = FALSE, useDegree = TRUE, momentum = 0.5))
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
	neighbors <- randomProjectionTreeSearch(dat, K = ncol(dat) - 1, max_iter = 10)
	edges <- buildEdgeMatrix(dat, neighbors)
	d2 <- as.dist(edges)
	expect_true(inherits(d2, "dist"))
	expect_true(inherits(d2, "dissimilarity"))
	expect_equal(as.matrix(do), as.matrix(d2))
	expect_equal(attr(d2, "Metric"), "euclidean")
	expect_equal(sum(is.na(as.matrix(d2))), 0)
})

test_that("build edge matrix as distance matches dist with nas", {
	neighbors <- randomProjectionTreeSearch(dat, K = 20)
	edges <- buildEdgeMatrix(dat, neighbors)
	d3 <- as.dist(edges)
	todelete <- !is.na(as.matrix(d3))
	expect_equal(as.matrix(do)[todelete], as.matrix(d3)[todelete])
	expect_equal(attr(d3, "Metric"), "euclidean")
})
