context("all together and variations")

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
                   K = 10, verbose = FALSE, distance_method = "Cosine")
  expect_equal(sum(any(is.na(visObject$coords)) +
                     any(is.nan(visObject$coords)) +
                     any(is.infinite(visObject$coords))),
               0)
})
