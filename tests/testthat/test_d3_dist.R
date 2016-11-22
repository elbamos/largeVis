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