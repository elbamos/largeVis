context("sparse")
set.seed(1974)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)

test_that("buildEdgeMatrix are the same, Euclidean", {
  neighbors <- randomProjectionTreeSearch(dat,
                                          K = 5,
                                          n_trees = 10,
                                          tree_threshold = 20,
                                          max_iter = 10,
                                          verbose = FALSE)
  edges1 <- buildEdgeMatrix(data = dat, neighbors = neighbors, verbose = FALSE)
  edges2 <- buildEdgeMatrix(data = Matrix::Matrix(dat, sparse = TRUE), neighbors = neighbors, verbose = FALSE)
  score <- sum(edges1@x - edges2@x)
  expect_lt(score, 1)
})

test_that("buildEdgeMatrix are the same, Cosine", {
	neighbors <- randomProjectionTreeSearch(dat,
																					K = 5,
																					n_trees = 10,
																					tree_threshold = 20,
																					max_iter = 10,
																					verbose = FALSE)
	edges1 <- buildEdgeMatrix(data = dat, neighbors = neighbors, verbose = FALSE, distance_method = "Cosine")
	edges2 <- buildEdgeMatrix(data = Matrix::Matrix(dat, sparse = TRUE), neighbors = neighbors, verbose = FALSE, distance_method = "Cosine")
	score <- sum(edges1@x - edges2@x)
	expect_lt(score, 1)
})

M <- 5
mat <- Matrix::sparseMatrix(i = rep(1:nrow(dat), ncol(dat)),
														j = rep(1:ncol(dat), each = nrow(dat)),
														x = as.vector(dat))
d = as.matrix(dist(t(as.matrix(mat)), method = "euclidean"))

test_that("sparseDistances", {
  index_matrix <- matrix(c(
    rep(0:(ncol(dat) - 1), ncol(dat)),
    rep(0:(ncol(dat) - 1), each = ncol(dat))
  ), ncol = 2, byrow = FALSE)
  new_distances <- distance(mat,
                            as.vector(index_matrix[, 2]),
                            as.vector(index_matrix[, 1]),
                            "Euclidean",
                            verbose = FALSE)
  diffs <- as.vector(d) - new_distances
  expect_lt(sum(diffs), 1e-10)
})

test_that("Can determine sparse iris neighbors accurately", {
  d <- as.matrix(dist(t(as.matrix(mat)), method = "euclidean"))
  bests <- apply(d, MARGIN = 1, FUN = function(x) order(x)[1:(M + 1)])
  bests <- bests[-1,] - 1
  neighbors <- randomProjectionTreeSearch(mat,
                                          K = M,
                                          n_trees = 10,
                                          max_iter = 2,
                                          tree_threshold = 20,
                                          verbose = FALSE)
  expect_lte(sum(neighbors - bests, na.rm = TRUE), 5)
})
