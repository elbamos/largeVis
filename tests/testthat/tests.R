context("largeVis")

test_that("Can determine iris neighbors", {
  data (iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes = which(duplicated(dat))
  dat <- dat[-dupes,]
  neighbors <- randomProjectionTreeSearch(dat, K = 5, n_trees = 10, tree_threshold = 20, max_iter = 10,
                                          verbose = F)
  expect_equal(nrow(neighbors), 5)
  expect_equal(ncol(neighbors), nrow(dat))
  expect_equal(sum(neighbors == -1),0)
  expect_equal(sum(neighbors[1:40] > 50), 0)
})

test_that("largeVis works", {
  data(iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes = which(duplicated(dat))
  dat <- dat[-dupes,]
  visObject <- vis(dat, max_iter = 20, sgd_batches = 1000,
                        K = 10,  gamma = 2, rho = 1, M = 40, alpha = 20,verbose=F)
  expect_equal(sum(any(is.na(visObject$coords)) + any(is.nan(visObject$coords)) + any(is.infinite(visObject$coords))), 0)
})
