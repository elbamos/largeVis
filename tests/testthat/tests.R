context("Edge Matrix")

test_that("Edge Matrix doesn't crash", {
	M <- 10
	data (quakes)
	dat <- t(scale(as.matrix(quakes)))
	neighbors <- randomProjectionTreeSearch(dat,
																						K = 20,  threads = 2, verbose = FALSE)
	expect_silent(edges <- buildEdgeMatrix(dat, neighbors, verbose = FALSE))
})

test_that("build edge matrix works as expected", {
	data (iris)
	set.seed(1974)
	dat <- as.matrix(iris[, 1:4])
	dat <- scale(dat)
	dupes <- which(duplicated(dat))
	dat <- dat[-dupes, ]
	dat <- t(dat)
	neighbors <- randomProjectionTreeSearch(dat, K = 20, threads = 2)
	expect_silent(edges <- buildEdgeMatrix(dat, neighbors))
})

context("wij")

test_that("wij doesn't crash", {
	data (iris)
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

test_that("project knns doesn't crash with momentum", {
	data (iris)
	set.seed(1974)
	dat <- as.matrix(iris[, 1:4])
	dat <- scale(dat)
	dupes <- which(duplicated(dat))
	dat <- dat[-dupes, ]
	dat <- t(dat)
	neighbors <- randomProjectionTreeSearch(dat, K = 20, threads = 2)
	edges <- buildEdgeMatrix(dat, neighbors)
	wij <- buildWijMatrix(edges, threads = 2)
	expect_silent(coords <- projectKNNs(wij, sgd_batches = 100, verbose = FALSE, momentum = 0.5, threads = 2))
})

test_that("project knns doesn't crash", {
	data (iris)
	set.seed(1974)
	dat <- as.matrix(iris[, 1:4])
	dat <- scale(dat)
	dupes <- which(duplicated(dat))
	dat <- dat[-dupes, ]
	dat <- t(dat)
	neighbors <- randomProjectionTreeSearch(dat, K = 20, threads = 2)
	edges <- buildEdgeMatrix(dat, neighbors)
	wij <- buildWijMatrix(edges, threads = 2)
	coords <- projectKNNs(wij, sgd_batches = 100, verbose = FALSE, threads = 2)
})

context("sgd batches")

test_that("sgd batches is linear with E", {
	expect_equal(sgdBatches(1e6, 10), sgdBatches(1e6, 1e4))
})
test_that("sgd batches has a small bump", {
	expect_lt(sgdBatches(5000), sgdBatches(10001))
	expect_lt(sgdBatches(9999) / 1.5, sgdBatches(10001))
})

