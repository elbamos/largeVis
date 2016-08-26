context("specific issue tests")

test_that("dim064 issue is resolved", {
	load(system.file("extdata/badmat.Rda", package = "largeVis"))
	expect_silent(neighbors <- randomProjectionTreeSearch(badmat, K = 50, threads = 2))
	expect_silent(edges <- buildEdgeMatrix(badmat, neighbors = neighbors, threads = 2))
	expect_silent(wij <- buildWijMatrix(edges, threads = 2))
	expect_silent(vis <- largeVis(badmat, K = 50, threads = 2, sgd_batches = 1000))
})

test_that("dim512 issue (excessive distances in the edge matrix)", {
	skip_on_cran()
	load(system.file("extdata/dim512.Rda", package = "largeVis"))
	dat <- t(scale(as.matrix(dim512)))
	expect_warning(vis <- largeVis(dat, K = 100, threads = 2, sgd_batches = 1000))
})