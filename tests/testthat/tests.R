context("largeVis")

test_that("Can determine iris neighbors", {
  data (iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes = which(duplicated(dat))
  dat <- dat[-dupes,]
  dat <- t(dat)
  neighbors <- randomProjectionTreeSearch(dat, K = 5, n_trees = 10, tree_threshold = 20, max_iter = 10,
                                          verbose = FALSE)
  expect_equal(nrow(neighbors), 5)
  expect_equal(ncol(neighbors), ncol(dat))
  expect_equal(sum(neighbors == -1),0)
  expect_equal(sum(neighbors[,1:40] > 50), 0)
})

test_that("largeVis works", {
  data(iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes = which(duplicated(dat))
  dat <- dat[-dupes,]
  dat <- t(dat)
  visObject <- vis(dat, max_iter = 20, sgd_batches = 1000,
                        K = 10,  gamma = 0.5, verbose=FALSE)
  expect_equal(sum(any(is.na(visObject$coords)) + any(is.nan(visObject$coords)) + any(is.infinite(visObject$coords))), 0)
  expect_equal(sum(any(visObject$coords > 10) + any(visObject$coords < -10)), 0)
})

test_that("largeVis works without weights", {
  data(iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes = which(duplicated(dat))
  dat <- dat[-dupes,]
  dat <- t(dat)
  visObject <- vis(dat, max_iter = 20, sgd_batches = 1000, weight_pos_samples = FALSE,
                   K = 10, verbose=FALSE)
  expect_equal(sum(any(is.na(visObject$coords)) + any(is.nan(visObject$coords)) + any(is.infinite(visObject$coords))), 0)
  expect_equal(sum(any(visObject$coords > 10) + any(visObject$coords < -10)), 0)
})

test_that("largeVis works when alpha = 0", {
  data(iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes = which(duplicated(dat))
  dat <- dat[-dupes,]
  dat <- t(dat)
  visObject <- vis(dat, alpha = 0,verbose=FALSE, sgd_batches=1000)
  expect_equal(sum(any(is.na(visObject$coords)) + any(is.nan(visObject$coords)) + any(is.infinite(visObject$coords))), 0)
})


