context("specific issue tests")

test_that("dim064 issue is resolved (large numbers of duplicate points)", {
	load(system.file("testdata/badmat.Rda", package = "largeVis"))
	badmat <- badmat
	expect_silent(neighbors <- randomProjectionTreeSearch(x = badmat, K = 50))
	expect_silent(edges <- buildEdgeMatrix(data = badmat, neighbors = neighbors, threads = 2))
	expect_silent(vis <- largeVis(badmat, K = 50, sgd_batches = 1000))
})

test_that("dim512 issue (excessive distances in the edge matrix)", {
	load(system.file("testdata/dim512.Rda", package = "largeVis"))
	dim512 <- dim512
	dat <- t(scale(as.matrix(dim512)))
	expect_warning(vis <- largeVis(dat, K = 100, sgd_batches = 1000))
})

test_that("neighbors does not fail with 0 max iters when the neighborhood is complete", {
	set.seed(1974)
	data(iris)
	dat <- as.matrix(iris[, 1:4])
	dupes <- which(duplicated(dat))
	dat <- dat[-dupes, ]
	dat <- t(dat)
	K <- 148
	distances <- as.matrix(dist(t(dat)))
	expect_silent(neighbors <- randomProjectionTreeSearch(dat, K = K, tree_threshold = 40,
																												n_trees = 10, max_iter = 0, verbose = FALSE))
})

test_that("sparse division by zero is resolved", {
	skip("big test")
	load(system.file("testdata/zerotest.rda", package = "largeVis"))
	dat <- Matrix::Matrix(zerotest, sparse = TRUE)
	expect_silent(ted <- randomProjectionTreeSearch(dat, K = 50, max_iter = 1, n_trees = 10))
})