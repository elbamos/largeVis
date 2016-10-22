context("neighbors")

data (iris)
set.seed(1974)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)

test_that("Trees does not error", {
	expect_silent(neighbors <- randomProjectionTreeSearch(dat,
																												K = 5,
																												n_trees = 10,
																												tree_threshold = 30,
																												max_iter = 0, threads = 1,
																												verbose = FALSE))

})

M <- 5
d_matrix <- as.matrix(dist(t(dat), method = "euclidean"))
bests <- apply(d_matrix, MARGIN = 1, FUN = function(x) order(x)[1:(M + 1)])
bests <- bests[-1,] - 1

test_that("max threshold is sufficient to find all neighbors", {
	neighbors <- randomProjectionTreeSearch(dat,
																					K = M,
																					n_trees = 1,
																					tree_threshold = ncol(dat),
																					max_iter = 0, threads = 2,
																					verbose = FALSE)
	scores <- lapply(1:ncol(dat), FUN = function(x) sum(neighbors[,x] %in% bests[,x]))
	score <- sum(as.numeric(scores))
	expect_gte(score, M * ncol(dat) - 1) # Two neighbors are equidistanct
})

test_that("Trees does not error if neighbors are explored once", {
	neighbors <- randomProjectionTreeSearch(dat,
																												K = 5,
																												n_trees = 50,
																												tree_threshold = 20,
																												max_iter = 1, threads = 2,
																												verbose = FALSE)
})

test_that("Trees does not error if neighbors are explored more than once", {
	expect_silent(neighbors <- randomProjectionTreeSearch(dat,
																												K = 5,
																												n_trees = 50,
																												tree_threshold = 20,
																												max_iter = 2, threads = 2,
																												verbose = FALSE))
})

test_that("exploring once after max threshold does not reduce accuracy", {
	neighbors <- randomProjectionTreeSearch(dat,
																					K = M,
																					n_trees = 1,
																					tree_threshold = ncol(dat),
																					max_iter = 0, threads = 2,
																					verbose = FALSE)
	scores <- lapply(1:ncol(dat), FUN = function(x) sum(neighbors[, x] %in% bests[, x]))
	score <- sum(as.numeric(scores))
	expect_gte(score, (M * ncol(dat)) - 1)
	oldscore <- score
	neighbors <- randomProjectionTreeSearch(dat,
																					K = M,
																					n_trees = 1,
																					tree_threshold = ncol(dat),
																					max_iter = 1, threads = 2,
																					verbose = FALSE)
	scores <- lapply(1:ncol(dat), FUN = function(x) sum(neighbors[, x] %in% bests[, x]))
	score <- sum(as.numeric(scores))
	expect_gte(score, oldscore)
})


test_that("exploring more than once after max threshold does not reduce accuracy", {
	neighbors <- randomProjectionTreeSearch(dat,
																					K = M,
																					n_trees = 1,
																					tree_threshold = ncol(dat),
																					max_iter = 0, threads = 2,
																					verbose = FALSE)
	scores <- lapply(1:ncol(dat), FUN = function(x) sum(neighbors[, x] %in% bests[, x]))
	score <- sum(as.numeric(scores))
	expect_gte(score, (M * ncol(dat)) - 1)
	oldscore <- score
	neighbors <- randomProjectionTreeSearch(dat,
																					K = M,
																					n_trees = 1,
																					tree_threshold = ncol(dat),
																					max_iter = 2, threads = 2,
																					verbose = FALSE)
	scores <- lapply(1:ncol(dat), FUN = function(x) sum(neighbors[, x] %in% bests[, x]))
	score <- sum(as.numeric(scores))
	expect_gte(score, oldscore)
})


test_that("Can determine iris neighbors", {
	neighbors <- randomProjectionTreeSearch(dat,
																					K = 5,
																					n_trees = 20,
																					tree_threshold = 30,
																					max_iter = 10,
																					threads = 2,
																					verbose = FALSE)
	expect_equal(nrow(neighbors), 5)
	expect_equal(ncol(neighbors), ncol(dat))
	expect_lt(sum(neighbors == -1), 20)
	expect_equal(sum(neighbors[, 1:40] > 50), 0)
	scores <- lapply(1:ncol(dat), FUN = function(x) sum(neighbors[,x] %in% bests[,x]))
	score <- sum(as.numeric(scores))
	expect_gte(score, M * ncol(dat) - 1) # Two neighbors are equidistanct
})

test_that("Can determine iris neighbors accurately, Euclidean", {
	neighbors <- randomProjectionTreeSearch(dat,
																					K = M,
																					n_trees = 20,
																					tree_threshold = 30,
																					max_iter = 12,
																					verbose = FALSE,  threads = 2,
																					seed = 1974)
	expect_lte(sum(neighbors != bests, na.rm = TRUE), 5)
})

M <- 10
data (quakes)
dat <- as.matrix(quakes)
quakes <- scale(dat)
d_matrix = as.matrix(dist(quakes, method = "euclidean"))

test_that("With a bigger dataset, increasing threshold improves result", {
	bests <- apply(d_matrix, MARGIN = 1, FUN = function(x) order(x)[1:(M + 1)])
	bests <- bests[-1, ] - 1

	oldscore <- nrow(quakes) * M

	for (t in c(10, 40, 80, 160)) {
		set.seed(1974)
		neighbors <- randomProjectionTreeSearch(t(quakes),
																						K = M,
																						n_trees = 20,
																						tree_threshold = t,
																						max_iter = 0,
																						verbose = FALSE, threads = 2,
																						seed = 1974)
		score <- sum(neighbors != bests, na.rm = TRUE)
		expect_lte(score, oldscore, label = t)
		oldscore <- score
	}
})

test_that("With a bigger dataset, increasing n_trees improves result", {
	bests <- apply(d_matrix, MARGIN=1, FUN = function(x) order(x)[1:(M + 1)])
	bests <- bests[-1,] - 1

	oldscore <- nrow(quakes) * M

	for (t in c(5, 20, 40, 90)) {
		neighbors <- randomProjectionTreeSearch(t(quakes),
																						K = M,
																						n_trees = t,
																						tree_threshold = 10,
																						max_iter = 0,
																						verbose = FALSE, threads = 2,
																						seed = 1974)
		score <- sum(neighbors != bests, na.rm = TRUE)
		expect_lte(score, oldscore, label = t)
		oldscore <- score
	}
})

test_that("With a bigger dataset, increasing iters improves result", {
	M <- 10
	bests <- apply(d_matrix, MARGIN = 1, FUN = function(x) order(x)[1:(M + 1)])
	bests <- bests[ - 1,] - 1

	oldscore <- nrow(quakes) * M

	for (t in c(0, 5, 20, 40)) {
		neighbors <- randomProjectionTreeSearch(t(quakes),
																						K = M,
																						n_trees = 5,
																						tree_threshold = 10,
																						max_iter = t,
																						verbose = FALSE, threads = 2,
																						seed = 1974)
		score <- max(0, sum(neighbors != bests, na.rm = TRUE))
		expect_lte(score, oldscore, label = t)
		oldscore <- score
	}
})

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

