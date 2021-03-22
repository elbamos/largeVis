context("sparse")


set.seed(1974)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
M <- 5
mat <- Matrix::sparseMatrix(i = rep(1:nrow(dat), ncol(dat)),
														j = rep(1:ncol(dat), each = nrow(dat)),
														x = as.vector(dat))
d = as.matrix(dist(t(as.matrix(mat)), method = "euclidean"))

RcppParallel::setThreadOptions(numThreads = 2)
test_that("Can determine sparse iris neighbors accurately", {
  bests <- apply(d, MARGIN = 1, FUN = function(x) order(x)[1:(M + 1)])
  bests <- bests[-1,] - 1
  neighbors <- randomProjectionTreeSearch(mat,
                                          K = M,
                                          n_trees = 20,
                                          max_iter = 2,
                                          verbose = FALSE)
  expect_lte(sum(neighbors$neighbors != bests, na.rm = TRUE), 6)
})
