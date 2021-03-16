context("vis")
set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)

RcppParallel::setThreadOptions(numThreads=2)

test_that("largeVis simple linux failure is fixed", {
	d2 <- t(as.matrix(iris[, 1:4]))
	expect_silent(vis <- largeVis(d2))
})

test_that("largeVis works", {
	visObject <- largeVis(dat, max_iter = 20, n_trees = 100,
											  sgd_batches = 1000,
												K = 20,  verbose = FALSE)
	expect_false(any(is.na(visObject$coords)))
	expect_false(any(is.nan(visObject$coords)))
	expect_false(any(is.infinite(visObject$coords)))
})

test_that("largeVis does not NaN on iris", {
	visObject <- largeVis(dat, max_iter = 20,
												coords = matrix(rnorm(ncol(dat) * 2), nrow = 2),
												K = 20,  verbose = FALSE,
												sgd_batches = 20000 * 150)
	expect_false(any(is.na(visObject$coords)))
	expect_false(any(is.nan(visObject$coords)))
	expect_false(any(is.infinite(visObject$coords)))
})

test_that("largeVis works when alpha == 0", {

	visObject <- largeVis(dat,
												max_iter = 20,
												sgd_batches = 10000,
												K = 10,
												alpha = 0,
												verbose = FALSE)
	expect_false(any(is.na(visObject$coords)))
	expect_false(any(is.nan(visObject$coords)))
	expect_false(any(is.infinite(visObject$coords)))
})

test_that("largeVis works with cosine", {

	visObject <- largeVis(dat, max_iter = 20,
												sgd_batches = 1000,
												K = 10, verbose = FALSE,
												distance_method = "Cosine")
	expect_false(any(is.na(visObject$coords)))
	expect_false(any(is.nan(visObject$coords)))
	expect_false(any(is.infinite(visObject$coords)))
})

test_that("largeVis graidents aren't off", {
	skip_on_cran()
	visObject <- largeVis(dat, K = 30, max_iter = 20, verbose = FALSE)
	expect_false(any(is.na(visObject$coords)))
	expect_false(any(is.nan(visObject$coords)))
	expect_false(any(is.infinite(visObject$coords)))
	expect_equal(sum(visObject$coords > 50), 0)
	expect_equal(sum(visObject$coords < -50), 0)
})

test_that("largeVis can eat a data.frame", {
	expect_silent(visObj <- largeVis(iris, K = 20, max_iter = 10, sgd_batches = 1, verbose = FALSE))
})

test_that("largeVis rejects wrong data types", {
	expect_error(visObj <- largeVis(letters, K = 2, max_iter = 10, sgd_batches = 1, verbose = FALSE))
	expect_error(visObj <- largeVis(matrix(sample(letters, 100, replace = T), nrow = 10, ncol = 10),
																	K = 2, max_iter = 10, sgd_batches = 1, verbose = FALSE))
})
