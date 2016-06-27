context("sparse")

test_that("buildEdgeMatrix are the same", {
  set.seed(1974)
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
                              index1 = FALSE)
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

test_that("sparseDistances", {
  M <- 5
  set.seed(1974)
  data (iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  mat <- Matrix::sparseMatrix(i = rep(1:nrow(dat), ncol(dat)),
                              j = rep(1:ncol(dat), each = nrow(dat)),
                              x = as.vector(dat))
  d = as.matrix(dist(mat, method = "euclidean"))
  index_matrix <- matrix(c(
    rep(0:(nrow(dat) - 1), nrow(dat)),
    rep(0:(nrow(dat) - 1), each = nrow(dat))
  ), ncol = 2, byrow = FALSE)
  mat <- Matrix::t(mat)
  new_distances <- distance(mat,
                            as.vector(index_matrix[, 2]),
                            as.vector(index_matrix[, 1]),
                            "Euclidean",
                            verbose = FALSE)
  diffs <- as.vector(d) - new_distances
  expect_lt(sum(diffs), 1e-10)
})

test_that("Can determine sparse iris neighbors accurately", {
  M <- 5
  set.seed(1974)
  data (iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  mat <- Matrix::sparseMatrix(i = rep(1:nrow(dat), ncol(dat)),
                              j = rep(1:ncol(dat), each = nrow(dat)),
                              x = as.vector(dat))
  d_matrix <- as.matrix(dist(mat, method = "euclidean"))
  bests <- apply(d_matrix, MARGIN = 1, FUN = function(x) order(x)[1:(M + 1)])
  bests <- bests[-1,] - 1
  mat <- Matrix::t(mat)
  neighbors <- randomProjectionTreeSearch(mat,
                                          K = M,
                                          n_trees = 10,
                                          tree_threshold = 20,
                                          max_iter = 2,
                                          verbose = FALSE)
  scores <- lapply(1:nrow(dat), FUN = function(x) sum(neighbors[, x] %in% bests[, x]))
  score <- sum(as.numeric(scores))
  expect_gt(score, .99 * ncol(dat) * M)
})
