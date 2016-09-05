skip_old_windows <- function() {
#	testthat::skip_if_not(R.Version()$arch != "i386", "largeVis does not run on 32-bit Windows.")
}

context("Edge Matrix")

test_that("Edge Matrix doesn't crash", {
	skip_old_windows()
	M <- 10
	data (quakes)
	dat <- t(scale(as.matrix(quakes)))
	neighbors <- randomProjectionTreeSearch(dat,
																						K = 20,  threads = 2, verbose = FALSE)
	expect_silent(edges <- buildEdgeMatrix(dat, neighbors, verbose = FALSE))
})

test_that("build edge matrix works as expected", {
	skip_old_windows()
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
	skip_old_windows()
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

test_that("project knns doesn't crash", {
	skip_old_windows()
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
	expect_silent(coords <- projectKNNs(wij, sgd_batches = 100, verbose = FALSE, threads = 2))
})
