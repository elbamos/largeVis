context("largeVis")

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
                                          n_trees = 10,
                                          tree_threshold = 20,
                                          max_iter = 10,
                                          verbose = FALSE)
  expect_equal(nrow(neighbors), 5)
  expect_equal(ncol(neighbors), ncol(dat))
  expect_equal(sum(neighbors == -1), 0)
  expect_equal(sum(neighbors[, 1:40] > 50), 0)
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
  d_matrix = as.matrix(dist(dat, method = 'euclidean'))
  bests <- apply(d_matrix, MARGIN=1, FUN = function(x) order(x)[1:(M + 1)])
  bests <- bests[-1,] - 1
  dat <- t(dat)
  neighbors <- randomProjectionTreeSearch(dat,
                                          K = M,
                                          n_trees = 10,
                                          tree_threshold = 20,
                                          max_iter = 2,
                                          verbose = FALSE)
  scores <- lapply(1:ncol(dat), FUN = function(x) sum(neighbors[,x] %in% bests[,x]))
  score <- sum(as.numeric(scores))
  expect_gt(score, .99 * ncol(dat) * M)
})

test_that("largeVis works", {
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  data(iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  dat <- t(dat)
  visObject <- vis(dat, max_iter = 20, sgd_batches = 1000,
                        K = 10,  gamma = 0.5, verbose = FALSE)
  expect_equal(sum(any(is.na(visObject$coords)) +
                     any(is.nan(visObject$coords)) +
                     any(is.infinite(visObject$coords))),
               0)
})

test_that("largeVis works without weights", {
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  data(iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  dat <- t(dat)
  visObject <- vis(dat,
                   max_iter = 20,
                   sgd_batches = 1000,
                   weight_pos_samples = FALSE,
                   K = 10,
                   verbose = FALSE)
  expect_equal(sum(any(is.na(visObject$coords)) +
                     any(is.nan(visObject$coords)) +
                     any(is.infinite(visObject$coords))),
               0)
})

test_that("largeVis works with cosine", {
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  data(iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  dat <- t(dat)
  visObject <- vis(dat, max_iter = 20, sgd_batches = 1000,
                   K = 10, verbose = FALSE, distance_method="Cosine")
  expect_equal(sum(any(is.na(visObject$coords)) +
                     any(is.nan(visObject$coords)) +
                     any(is.infinite(visObject$coords))),
               0)
})

test_that("Euclidean distances are correct", {
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  test_matrix <- matrix(rnorm(100), nrow = 10)
  distances <- as.matrix(dist(test_matrix, method = "euclidean"))
  index_matrix <- matrix(c(rep(0:9, each = 10), rep(0:9, 10)),
                         ncol = 2, byrow = FALSE)
  test_matrix <- t(test_matrix)
  new_distances <- distance(as.vector(index_matrix[,2]),
                             as.vector(index_matrix[,1]),
                             x = test_matrix,
                             "Euclidean",
                              verbose = FALSE)
  diffs <- as.vector(distances) - new_distances
  expect_lt(sum(diffs), 1e-10)
})

test_that("Cosine distances are correct", {
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  cos.sim <- function(x, i, j) {
    A = x[,i]
    B = x[,j]
    return( 1 - (sum(A*B)/sqrt(sum(A^2)*sum(B^2)) ))
  }
  test_matrix <- matrix(rnorm(100), nrow = 10)
  index_matrix <- matrix(c(rep(0:9, each = 10), rep(0:9, 10)),
                         ncol = 2, byrow = FALSE)
  distances <- apply(index_matrix + 1,
                     MARGIN=1,
                     FUN = function(x) cos.sim(test_matrix, x[1], x[2]))
  new_distances <- distance(as.vector(index_matrix[,2]),
                                       as.vector(index_matrix[,1]),
                                       x = test_matrix,
                                       "Cosine",
                                       verbose = FALSE)
  diffs <- as.vector(distances) - new_distances
  expect_lt(sum(diffs), 1e-10)
})

test_that("buildEdgeMatrix are the same", {
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  dat <- t(dat)
  neighbors <- randomProjectionTreeSearch(dat,
                                          K = 5,
                                          n_trees = 10,
                                          tree_threshold = 20,
                                          max_iter = 10,
                                          verbose = FALSE)
  is <- rep(0:(ncol(dat) - 1), each = 5)
  js <- as.vector(neighbors)
  is <- is[! js == -1]
  js <- js[! js == -1]
  dupes <- duplicated(data.frame(is, js))
  is <- is[! dupes]
  js <- js[! dupes]
  ord <- order(is)
  is <- is[ord]
  js <- js[ord]
  distances <- as.matrix(dist(t(dat)))
  distances <- as.numeric(lapply(1:length(is), FUN = function(x) {
    distances[is[x] + 1, js[x] + 1]
  }))

  ps <- i2p(is)
  sigwij <- buildEdgeMatrix(i = is,
                            j = js,
                            p = ps,
                            d = distances,
                            verbose = F)
  mat <- Matrix::sparseMatrix(i = js,
                           p = ps,
                           x = distances,
                           dims = c(ncol(dat), ncol(dat)),
                           giveCsparse = TRUE,
                           index1=FALSE)
  sigwij2 <- buildEdgeMatrix(mat, verbose = F)

  score <- sum(sigwij$wij@x != sigwij2$wij@x)
  expect_lt(score, 450)
  tmat <- Matrix::sparseMatrix(i = is,
                               j = js,
                               x = distances,
                               dims = c(ncol(dat), ncol(dat)),
                               giveCsparse = FALSE,
                               index1 = FALSE)
  sigwij3 <- buildEdgeMatrix(tmat, verbose = F)
  score <- sum(sigwij$wij@x != sigwij3$wij@x)
  expect_lt(score, 450)
})

test_that("largeVis works when alpha == 0", {
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  data(iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  dat <- t(dat)
  visObject <- vis(dat,
                   max_iter = 20,
                   sgd_batches = 10000,
                   K = 10,
                   alpha = 0,
                   verbose = FALSE,
                   weight_pos_samples = FALSE)
  expect_equal(sum(any(is.na(visObject$coords)) +
                     any(is.nan(visObject$coords)) +
                     any(is.infinite(visObject$coords))),
               0)
})

test_that("sparseDistances", {
  M <- 5
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  data (iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  mat <- Matrix::sparseMatrix(i = rep(1:nrow(dat), ncol(dat)),
                              j = rep(1:ncol(dat), each = nrow(dat)),
                              x = as.vector(dat))
  d = as.matrix(dist(mat, method = 'euclidean'))
  index_matrix <- matrix(c(
    rep(0:(nrow(dat) - 1), nrow(dat)),
    rep(0:(nrow(dat) - 1), each = nrow(dat))
  ), ncol = 2, byrow = FALSE)
  mat <- Matrix::t(mat)
  new_distances <- distance(mat,
                            as.vector(index_matrix[,2]),
                            as.vector(index_matrix[,1]),
                            "Euclidean",
                            verbose = FALSE)
  diffs <- as.vector(d) - new_distances
  expect_lt(sum(diffs), 1e-10)
})

test_that("Can determine sparse iris neighbors accurately", {
  M <- 5
  set.seed(1974)
  RcppArmadillo::armadillo_set_seed(1974)
  data (iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  mat <- Matrix::sparseMatrix(i = rep(1:nrow(dat), ncol(dat)),
                              j = rep(1:ncol(dat), each = nrow(dat)),
                              x = as.vector(dat))
  d_matrix = as.matrix(dist(mat, method = 'euclidean'))
  bests <- apply(d_matrix, MARGIN=1, FUN = function(x) order(x)[1:(M + 1)])
  bests <- bests[-1,] - 1
  mat <- Matrix::t(mat)
  neighbors <- randomProjectionTreeSearch(mat,
                                          K = M,
                                          n_trees = 10,
                                          tree_threshold = 20,
                                          max_iter = 2,
                                          verbose = FALSE)
  scores <- lapply(1:nrow(dat), FUN = function(x) sum(neighbors[,x] %in% bests[,x]))
  score <- sum(as.numeric(scores))
  expect_gt(score, .99 * ncol(dat) * M)
})
