context("vis")
set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)

test_that("largeVis works", {
	visObject <- largeVis(dat, max_iter = 20, n_trees = 100,
												tree_threshold = 50, sgd_batches = 1000,  threads = 2,
												K = 20,  verbose = FALSE)
	expect_false(any(is.na(visObject$coords)))
	expect_false(any(is.nan(visObject$coords)))
	expect_false(any(is.infinite(visObject$coords)))
})

test_that("largeVis does not NaN on iris", {
	visObject <- largeVis(dat, max_iter = 20,
												coords = matrix(rnorm(ncol(dat) * 2), nrow = 2),  threads = 2,
												K = 20,  verbose = FALSE,
												sgd_batches = 20000 * 150)
	expect_false(any(is.na(visObject$coords)))
	expect_false(any(is.nan(visObject$coords)))
	expect_false(any(is.infinite(visObject$coords)))
})

test_that("largeVis works when alpha == 0", {
	visObject <- largeVis(dat,
												max_iter = 20,
												sgd_batches = 10000,  threads = 2,
												K = 10,
												alpha = 0,
												verbose = FALSE)
	expect_false(any(is.na(visObject$coords)))
	expect_false(any(is.nan(visObject$coords)))
	expect_false(any(is.infinite(visObject$coords)))
})

test_that("largeVis works with cosine", {
	visObject <- largeVis(dat, max_iter = 20,
												sgd_batches = 1000,  threads = 2,
												K = 10, verbose = FALSE,
												distance_method = "Cosine")
	expect_false(any(is.na(visObject$coords)))
	expect_false(any(is.nan(visObject$coords)))
	expect_false(any(is.infinite(visObject$coords)))
})

test_that("largeVis continues to work as it scales up", {
	visObject <- largeVis(dat, max_iter = 20, sgd_batches = 1000,  threads = 2,
												K = 10,  gamma = 0.5, verbose = FALSE)
	expect_false(any(is.na(visObject$coords)))
	expect_false(any(is.nan(visObject$coords)))
	expect_false(any(is.infinite(visObject$coords)))
	for (i in c(10000, 100000, 1000000, 20000 * length(visObject$wij@x))) {
		coords <- projectKNNs(visObject$wij, sgd_batches = i,  threads = 2,
													verbose = FALSE)
		expect_false(any(is.na(coords)))
		expect_false(any(is.nan(coords)))
		expect_false(any(is.infinite(coords)))
	}
})

context("specific issue tests")

test_that("dim064 issue is resolved", {
	load(system.file("extdata/badmat.Rda", package = "largeVis"))
	expect_silent(neighbors <- randomProjectionTreeSearch(badmat, K = 50, threads = 2))
	expect_silent(edges <- buildEdgeMatrix(badmat, neighbors = neighbors, threads = 2))
	expect_silent(wij <- buildWijMatrix(edges, threads = 2))
	expect_silent(vis <- largeVis(badmat, K = 50, sgd_batches = 1000))
})