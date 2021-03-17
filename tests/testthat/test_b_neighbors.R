context("neighbors")

data(iris)
set.seed(1974)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)

test_that("Trees does not error", {
	RcppParallel::setThreadOptions(numThreads = 1)
	expect_silent(randomProjectionTreeSearch(dat,
																												K = 5,
																												n_trees = 10,
																												max_iter = 0,
																												verbose = FALSE))

	RcppParallel::setThreadOptions(numThreads = 2)
	expect_silent(randomProjectionTreeSearch(dat,
																												K = 5,
																												n_trees = 10,
																												max_iter = 0,
																												verbose = FALSE))
	RcppParallel::setThreadOptions(numThreads = 1)
	expect_silent(randomProjectionTreeSearch(dat,
																												K = 5,
																												n_trees = 50,
																												max_iter = 1,
																												verbose = FALSE))
	expect_silent(randomProjectionTreeSearch(dat,
																												K = 5,
																												n_trees = 50,
																												max_iter = 1,
																												verbose = FALSE))
	RcppParallel::setThreadOptions(numThreads = 2)
	expect_silent(randomProjectionTreeSearch(dat,
																												K = 5,
																												n_trees = 50,
																												max_iter = 2,
																												verbose = FALSE))
})

RcppParallel::setThreadOptions(numThreads = 2)
test_that("Can use an on-disk index", {
	filename <- tempfile(pattern = "largevistest")
	expect_silent(randomProjectionTreeSearch(dat,
																					K = 5,
																					n_trees = 10,
																					max_iter = 1,
																					save_file = filename,
																					verbose = FALSE))
	expect_true(file.exists(filename))
})


M <- 5
d_matrix <- as.matrix(dist(t(dat), method = "euclidean"))
bests <- apply(d_matrix, MARGIN = 1, FUN = function(x) order(x)[1:(M + 1)])
bests <- bests[-1,] - 1


RcppParallel::setThreadOptions(numThreads = 2)
test_that("exploration is not negative", {
	neighbors <- randomProjectionTreeSearch(dat,
																					K = M,
																					n_trees = 1,
																					max_iter = 0,
																					verbose = FALSE)
	scores <- lapply(1:ncol(dat), FUN = function(x) sum(neighbors[, x] %in% bests[, x]))
	score <- sum(as.numeric(scores))
	oldscore <- score
	neighbors <- randomProjectionTreeSearch(dat,
																					K = M,
																					n_trees = 1,
																					max_iter = 1,
																					verbose = FALSE)
	scores <- lapply(1:ncol(dat), FUN = function(x) sum(neighbors[, x] %in% bests[, x]))
	score <- sum(as.numeric(scores))
	expect_gte(score, oldscore, label = "1 iteration")
	oldscore <- score
	neighbors <- randomProjectionTreeSearch(dat,
																					K = M,
																					n_trees = 1,
																					max_iter = 2,
																					verbose = FALSE)
	scores <- lapply(1:ncol(dat), FUN = function(x) sum(neighbors[, x] %in% bests[, x]))
	score <- sum(as.numeric(scores))
	expect_gte(score, oldscore, label = "2 iterations")
})

RcppParallel::setThreadOptions(numThreads = 1)
test_that("Can determine iris neighbors with iterations 1 thread", {
	neighbors <- randomProjectionTreeSearch(dat,
																					K = 5,
																					n_trees = 20,
																					max_iter = 10,
																					verbose = FALSE)
	expect_equal(sum(is.na(neighbors)), 0)
	expect_equal(nrow(neighbors), 5)
	expect_equal(ncol(neighbors), ncol(dat))
	expect_lt(sum(neighbors == -1), 20)
	expect_equal(sum(neighbors[, 1:40] > 50), 0)
	scores <- lapply(1:ncol(dat), FUN = function(x) sum(neighbors[,x] %in% bests[,x]))
	score <- sum(as.numeric(scores))
	expect_gte(score, M * ncol(dat) - 2) # Two neighbors are equidistanct
})

RcppParallel::setThreadOptions(numThreads = 2)
test_that("Can determine iris neighbors with iterations 2 threads", {
	neighbors <- randomProjectionTreeSearch(dat,
																					K = 5,
																					n_trees = 20,
																					max_iter = 10,
																					verbose = FALSE)

	expect_equal(sum(is.na(neighbors)), 0)
	expect_equal(nrow(neighbors), 5)
	expect_equal(ncol(neighbors), ncol(dat))
	expect_lt(sum(neighbors == -1), 20)
	expect_equal(sum(neighbors[, 1:40] > 50), 0)
	scores <- lapply(1:ncol(dat), FUN = function(x) sum(neighbors[,x] %in% bests[,x]))
	score <- sum(as.numeric(scores))
	expect_gte(score, M * ncol(dat) - 2) # Two neighbors are equidistanct
})

RcppParallel::setThreadOptions(numThreads = 2)
test_that("Can determine iris neighbors accurately, Euclidean", {
	neighbors <- randomProjectionTreeSearch(dat,
																					K = M,
																					n_trees = 20,
																					max_iter = 12,
																					verbose = FALSE,
																					seed = 1974)
	expect_lte(sum(neighbors != bests, na.rm = TRUE), 5)
})

RcppParallel::setThreadOptions(numThreads = 1)
test_that("With a bigger dataset, performance is as expected", {
	M <- 10
	data(quakes)
	dat <- as.matrix(quakes)
	quakes <- scale(dat)
	d_matrix = as.matrix(dist(quakes, method = "euclidean"))


	bests <- apply(d_matrix, MARGIN = 1, FUN = function(x) order(x)[1:(M + 1)])
	bests <- bests[-1, ] - 1

	oldscore <- nrow(quakes) * M

	for (t in c(5, 20, 40, 90)) {
		RcppParallel::setThreadOptions(numThreads = 1)
		neighbors <- randomProjectionTreeSearch(t(quakes),
																						K = M,
																						n_trees = t,
																						max_iter = 0,
																						verbose = FALSE,
																						seed = 1974)
		score <- sum(neighbors != bests, na.rm = TRUE)
		expect_lte(score, oldscore, label = paste("n_trees=", t))
		oldscore <- score
	}

	oldscore <- nrow(quakes) * M

	for (t in c(1, 5, 10, 20)) {
		RcppParallel::setThreadOptions(numThreads = 1)
		neighbors <- randomProjectionTreeSearch(t(quakes),
																						K = M,
																						n_trees = 10,
																						max_iter = t,
																						verbose = FALSE,
																						seed = 1974)
		score <- max(0, sum(neighbors != bests, na.rm = TRUE))
		expect_lte(score, oldscore, label = paste("iters=", t))
		oldscore <- score
	}

})

