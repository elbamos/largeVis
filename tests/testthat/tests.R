context("neighbors")

test_that("Trees does not error", {
  data (iris)
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  dat <- t(dat)
  expect_silent(neighbors <- randomProjectionTreeSearch(dat,
                                                        K = 5,
                                                        n_trees = 10,
                                                        tree_threshold = 20,
                                                        max_iter = 0,
                                                        verbose = FALSE))

})

test_that("Trees does not error if neighbors are explored once", {
  data (iris)
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  dat <- t(dat)

  expect_silent(neighbors <- randomProjectionTreeSearch(dat,
                                                        K = 5,
                                                        n_trees = 50,
                                                        tree_threshold = 20,
                                                        max_iter = 1,
                                                        verbose = FALSE))

})

test_that("Trees does not error if neighbors are explored more than once", {
  data (iris)
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  dat <- t(dat)

  expect_silent(neighbors <- randomProjectionTreeSearch(dat,
                                                        K = 5,
                                                        n_trees = 50,
                                                        tree_threshold = 20,
                                                        max_iter = 2,
                                                        verbose = FALSE))
})

test_that("Can determine iris neighbors", {
  data (iris)
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  dat <- t(dat)
  neighbors <- randomProjectionTreeSearch(dat,
                                          K = 5,
                                          n_trees = 20,
                                          tree_threshold = 30,
                                          max_iter = 10,
                                          verbose = FALSE)
  expect_equal(nrow(neighbors), 5)
  expect_equal(ncol(neighbors), ncol(dat))
  expect_lt(sum(neighbors == -1), 20)
  expect_equal(sum(neighbors[, 1:40] > 50), 0)
})

test_that("max threshold is sufficient to find all neighbors", {
  M <- 5
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  data (iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  d_matrix <- as.matrix(dist(dat, method = "euclidean"))
  bests <- apply(d_matrix, MARGIN=1, FUN = function(x) order(x)[1:(M + 1)])
  bests <- bests[-1,] - 1
  dat <- t(dat)

  neighbors <- randomProjectionTreeSearch(dat,
                                          K = M,
                                          n_trees = 1,
                                          tree_threshold = ncol(dat),
                                          max_iter = 0,
                                          verbose = FALSE)
  scores <- lapply(1:ncol(dat), FUN = function(x) sum(neighbors[,x] %in% bests[,x]))
  score <- sum(as.numeric(scores))
  expect_gte(score, M * ncol(dat) - 1) # Two neighbors are equidistanct
})

test_that("exploring after max threshold does not reduce accuracy", {
  M <- 5
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  data (iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  d_matrix <- as.matrix(dist(dat, method = "euclidean"))
  bests <- apply(d_matrix, MARGIN = 1, FUN = function(x) order(x)[1:(M + 1)])
  bests <- bests[-1, ] - 1
  dat <- t(dat)

  neighbors <- randomProjectionTreeSearch(dat,
                                          K = M,
                                          n_trees = 1,
                                          tree_threshold = ncol(dat),
                                          max_iter = 1,
                                          verbose = FALSE)
  scores <- lapply(1:ncol(dat), FUN = function(x) sum(neighbors[, x] %in% bests[, x]))
  score <- sum(as.numeric(scores))
  expect_gte(score, (M * ncol(dat)) - 1)
  oldscore <- score

  neighbors <- randomProjectionTreeSearch(dat,
                                          K = M,
                                          n_trees = 1,
                                          tree_threshold = ncol(dat),
                                          max_iter = 5,
                                          verbose = FALSE)
  scores <- lapply(1:ncol(dat), FUN = function(x) sum(neighbors[, x] %in% bests[, x]))
  score <- sum(as.numeric(scores))
  expect_gte(score, oldscore)
})

test_that("Can determine iris neighbors accurately", {
  M <- 5
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  data (iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  d_matrix <- as.matrix(dist(dat, method = "euclidean"))
  bests <- apply(d_matrix, MARGIN=1, FUN = function(x) order(x)[1:(M + 1)])
  bests <- bests[-1, ] - 1
  dat <- t(dat)

  neighbors <- randomProjectionTreeSearch(dat,
                                          K = M,
                                          n_trees = 10,
                                          tree_threshold = 10,
                                          max_iter = 10,
                                          verbose = FALSE)
  scores <- lapply(1:ncol(dat),
                   FUN = function(x) sum(neighbors[, x] %in% bests[, x]))
  score <- sum(as.numeric(scores))
  expect_gt(score, (ncol(dat) * M) - 15)
})

# test_that("Knows how to converge", {
#   M <- 5
#   set.seed(1974)
#   RcppArmadillo::armadillo_set_seed(1974)
#   data (iris)
#   dat <- as.matrix(iris[, 1:4])
#   dat <- scale(dat)
#   dupes <- which(duplicated(dat))
#   dat <- dat[-dupes, ]
#   d_matrix = as.matrix(dist(dat, method = "euclidean"))
#   bests <- apply(d_matrix, MARGIN=1, FUN = function(x) order(x)[1:(M + 1)])
#   bests <- bests[-1,] - 1
#   dat <- t(dat)
#
#   neighbors <- randomProjectionTreeSearch(dat,
#                                           K = M,
#                                           n_trees = 10,
#                                           tree_threshold = 10,
#                                           max_iter = 100000,
#                                           verbose = FALSE)
#   scores <- lapply(1:ncol(dat), FUN = function(x) sum(neighbors[,x] %in% bests[,x]))
#   score <- sum(as.numeric(scores))
#   expect_gt(score, (ncol(dat) * M) - 15)
# })



test_that("With a bigger dataset, increasing threshold improves result", {
  M <- 10
  data (quakes)
  dat <- as.matrix(quakes)
  dat <- scale(dat)
  d_matrix = as.matrix(dist(dat, method = "euclidean"))
  bests <- apply(d_matrix, MARGIN = 1, FUN = function(x) order(x)[1:(M + 1)])
  bests <- bests[-1, ] - 1
  dat <- t(dat)

  oldscore <- 0

  for (t in c(10, 30, 60, 90)) {
    set.seed(1974)
    RcppArmadillo::armadillo_set_seed(1974)
    neighbors <- randomProjectionTreeSearch(dat,
                                            K = M,
                                            n_trees = 10,
                                            tree_threshold = t,
                                            max_iter = 0,
                                            verbose = FALSE)
    scores <- lapply(1:ncol(dat),
                     FUN = function(x) sum(neighbors[, x] %in% bests[, x]))
    score <- sum(as.numeric(scores)) / (M * ncol(dat))
    expect_gte(score, oldscore * 0.99)  # Allow some gap here to account for randomness
    if (score == 1) break;
    oldscore <- score
  }
})

test_that("With a bigger dataset, increasing n_trees improves result", {
  M <- 10
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  data (quakes)
  dat <- as.matrix(quakes)
  dat <- scale(dat)
  d_matrix = as.matrix(dist(dat, method = "euclidean"))
  bests <- apply(d_matrix, MARGIN=1, FUN = function(x) order(x)[1:(M + 1)])
  bests <- bests[-1,] - 1
  dat <- t(dat)

  oldscore <- 0

  for (t in c(10, 30, 60, 90)) {
    neighbors <- randomProjectionTreeSearch(dat,
                                            K = M,
                                            n_trees = t,
                                            tree_threshold = 10,
                                            max_iter = 0,
                                            verbose = FALSE)
    scores <- lapply(1:ncol(dat),
                     FUN = function(x) sum(neighbors[, x] %in% bests[, x]))
    score <- sum(as.numeric(scores)) / (M * ncol(dat))
    expect_gte(score, oldscore * 0.99)
    if (score == 1) break;
    oldscore <- score
  }
})

test_that("With a bigger dataset, increasing iters improves result", {
  M <- 10
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  data (quakes)
  dat <- as.matrix(quakes)
  dat <- scale(dat)
  d_matrix = as.matrix(dist(dat, method = "euclidean"))
  bests <- apply(d_matrix, MARGIN = 1, FUN = function(x) order(x)[1:(M + 1)])
  bests <- bests[ - 1,] - 1
  dat <- t(dat)

  oldscore <- 0

  for (t in c(0, 1, 5, 10)) {
    neighbors <- randomProjectionTreeSearch(dat,
                                            K = M,
                                            n_trees = 10,
                                            tree_threshold = 10,
                                            max_iter = t,
                                            verbose = FALSE)
    scores <- lapply(1:ncol(dat),
                     FUN = function(x) sum(neighbors[, x] %in% bests[, x]))
    score <- sum(as.numeric(scores)) / (M * ncol(dat))
    expect_gte(score, oldscore * 0.99)
    if (score == 1) break;
    oldscore <- score
  }
})
