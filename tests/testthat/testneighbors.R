context("neighbor sorting")
set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
K <- 20
distances <- as.matrix(dist(t(dat)))

check <- function(adjacency, distance) {
	good <- 0
	bad <- 0
	for (i in 1:ncol(adjacency)) {
		for (j in 2:nrow(adjacency)) {
			if (adjacency[j, i] >= 0) {
				if (distance[i, adjacency[j, i] + 1] < distance[i, adjacency[j - 1, i] + 1]) bad <- bad + 1
				else good <- good + 1
			}
		}
	}
	list(good = good, bad = bad)
}

test_that("20 neighbors are sorted when max_iter = 1", {
	neighbors <- randomProjectionTreeSearch(dat, K = K, tree_threshold = 20, n_trees = 10,  threads = 2, max_iter = 1, verbose = FALSE)
	res <- check(neighbors, distances)
	expect_equal(res$bad, 0)
})

K <- 40

test_that("40 neighbors are sorted when max_iter = 1", {
	neighbors <- randomProjectionTreeSearch(dat, K = K, tree_threshold = 40, n_trees = 10,  threads = 2, max_iter = 1, verbose = FALSE)
	res <- check(neighbors, distances)
	expect_equal(res$bad, 0)
})

K <- 120

test_that("120 neighbors are sorted when max_iter = 1", {
	neighbors <- randomProjectionTreeSearch(dat, K = K, tree_threshold = 60, n_trees = 10,  threads = 2, max_iter = 1, verbose = FALSE)
	res <- check(neighbors, distances)
	expect_equal(res$bad, 0)
})

K <- 147

test_that("147 neighbors are sorted when max_iter = 1", {
	neighbors <- randomProjectionTreeSearch(dat, K = K, tree_threshold = 70, n_trees = 10,  threads = 2, max_iter = 1, verbose = FALSE)
	res <- check(neighbors, distances)
	expect_equal(res$bad, 0)
})

test_that("147 neighbors are sorted when max_iter = 1 and threshodl = 140", {
	neighbors <- randomProjectionTreeSearch(dat, K = K, tree_threshold = 140, n_trees = 10,  threads = 2, max_iter = 1, verbose = FALSE)
	res <- check(neighbors, distances)
	expect_equal(res$bad, 0)
})

test_that("neighbors are sorted when max_iter = 0", {
	neighbors <- randomProjectionTreeSearch(dat, K = K, tree_threshold = 40, n_trees = 10,  threads = 2, max_iter = 0, verbose = FALSE)
	res <- check(neighbors, distances)
	expect_equal(res$bad, 0)
})

test_that("neighbors are sorted when max_iter = 2", {
	neighbors <- randomProjectionTreeSearch(dat, K = K, tree_threshold = 40, n_trees = 10,  threads = 2, max_iter = 1, verbose = FALSE)
	res <- check(neighbors, distances)
	expect_equal(res$bad, 0)
})

