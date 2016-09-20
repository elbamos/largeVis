context("neighbor sorting")
set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
K <- 20
distances <- as.matrix(dist(t(dat)))

convert <- function(x, edges) {
	apply(1:(nrow(x) - 1), MARGIN = 2, FUN = function(y) {
		edges[x[y, ], ] - edges[x[y + 1, ], ]
	})
}

test_that("20 neighbors are sorted when max_iter = 1", {
	neighbors <- randomProjectionTreeSearch(dat, K = K, tree_threshold = 20, n_trees = 10,  threads = 2, max_iter = 1, verbose = FALSE)
	good <- (neighbors != -1)[- nrow(neighbors), ]
	results <- convert(neighbors, distances)
	expect_equal(sum(results[good] < 0), 0)
})

K <- 40

test_that("40 neighbors are sorted when max_iter = 1", {
	neighbors <- randomProjectionTreeSearch(dat, K = K, tree_threshold = 40, n_trees = 10,  threads = 2, max_iter = 1, verbose = FALSE)
	good <- (neighbors != -1)[- nrow(neighbors), ]
	biggers <- neighbors[-1, ]
	smallers <- neighbors[- nrow(neighbors), ]
	results <- biggers - smallers
	expect_equal(sum(results[good] < 0), 0)
})

K <- 120

test_that("120 neighbors are sorted when max_iter = 1", {
	neighbors <- randomProjectionTreeSearch(dat, K = K, tree_threshold = 60, n_trees = 10,  threads = 2, max_iter = 1, verbose = FALSE)
	good <- (neighbors != -1)[- nrow(neighbors), ]
	biggers <- neighbors[-1, ]
	smallers <- neighbors[- nrow(neighbors), ]
	results <- biggers - smallers
	expect_equal(sum(results[good] < 0), 0)
})

K <- 147

test_that("147 neighbors are sorted when max_iter = 1", {
	neighbors <- randomProjectionTreeSearch(dat, K = K, tree_threshold = 70, n_trees = 10,  threads = 2, max_iter = 1, verbose = FALSE)
	good <- (neighbors != -1)[- nrow(neighbors), ]
	biggers <- neighbors[-1, ]
	smallers <- neighbors[- nrow(neighbors), ]
	results <- biggers - smallers
	expect_equal(sum(results[good] < 0), 0)
})

test_that("neighbors are sorted when max_iter = 0", {
	neighbors <- randomProjectionTreeSearch(dat, K = K, tree_threshold = 40, n_trees = 10,  threads = 2, max_iter = 0, verbose = FALSE)
	good <- (neighbors != -1)[- nrow(neighbors), ]
	biggers <- neighbors[-1, ]
	smallers <- neighbors[- nrow(neighbors), ]
	results <- biggers - smallers
	expect_equal(sum(results[good] < 0), 0)
})

test_that("neighbors are sorted when max_iter = 2", {
	neighbors <- randomProjectionTreeSearch(dat, K = K, tree_threshold = 40, n_trees = 10,  threads = 2, max_iter = 1, verbose = FALSE)
	good <- (neighbors != -1)[- nrow(neighbors), ]
	biggers <- neighbors[-1, ]
	smallers <- neighbors[- nrow(neighbors), ]
	results <- biggers - smallers
	expect_equal(sum(results[good] < 0), 0)
})

