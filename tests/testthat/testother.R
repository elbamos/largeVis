context("all together and variations")

test_that("largeVis works", {
  set.seed(1974)
  data(iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  dat <- t(dat)

  visObject <- vis(dat, max_iter = 20, n_trees = 100, tree_threshold = 50, sgd_batches = 1000,
                   K = 20,  gamma = 0.5, verbose = FALSE)
  expect_false(any(is.na(visObject$coords)))
  expect_false(any(is.nan(visObject$coords)))
  expect_false(any(is.infinite(visObject$coords)))
})

test_that("largeVis does not NaN on iris", {
  set.seed(1974)
  data(iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  dat <- t(dat)

  visObject <- vis(dat, max_iter = 20, coords = matrix(rnorm(ncol(dat) * 2), nrow = 2),
                   K = 20,  gamma = 0.5, verbose = FALSE, sgd_batches = 20000 * 150)
  expect_false(any(is.na(visObject$coords)))
  expect_false(any(is.nan(visObject$coords)))
  expect_false(any(is.infinite(visObject$coords)))
})


test_that("largeVis works when alpha == 0", {
  set.seed(1974)
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
                   verbose = FALSE)
  expect_false(any(is.na(visObject$coords)))
  expect_false(any(is.nan(visObject$coords)))
  expect_false(any(is.infinite(visObject$coords)))
})

test_that("sigmas are well behaved", {
  set.seed(1974)
  data(iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  dat <- t(dat)

  neighbors <- randomProjectionTreeSearch(dat, verbose = FALSE)
  neighborIndices <- neighborsToVectors(neighbors)
  distances <- distance(x = dat,
                        i = neighborIndices$i,
                        j = neighborIndices$j, verbose = FALSE)
  wij <- buildEdgeMatrix(i = neighborIndices$i,
                         j = neighborIndices$j,
                         d = distances, verbose = FALSE)
  expect_false(any(is.infinite(wij$sigmas)))
  expect_false(any(is.nan(wij$sigmas)))
  expect_false(any(is.na(wij$sigmas)))
  expect_lt(max(wij$sigmas), 10)
})

test_that("largeVis works with cosine", {
  set.seed(1974)
  data(iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  dat <- t(dat)

  visObject <- vis(dat, max_iter = 20,
                   sgd_batches = 1000,
                   K = 10, verbose = FALSE,
                   distance_method = "Cosine")
  expect_false(any(is.na(visObject$coords)))
  expect_false(any(is.nan(visObject$coords)))
  expect_false(any(is.infinite(visObject$coords)))
})

test_that("largeVis continues to work as it scales up", {
  set.seed(1974)
  data(iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  dat <- t(dat)

  visObject <- vis(dat, max_iter = 20, sgd_batches = 1000,
                   K = 10,  gamma = 0.5, verbose = FALSE)
  expect_false(any(is.na(visObject$coords)))
  expect_false(any(is.nan(visObject$coords)))
  expect_false(any(is.infinite(visObject$coords)))
  for (i in c(10000, 100000, 1000000, 20000 * length(visObject$wij@x))) {
    coords <- projectKNNs(visObject$wij, sgd_batches = i,
                          verbose = FALSE)
    expect_false(any(is.na(coords)))
    expect_false(any(is.nan(coords)))
    expect_false(any(is.infinite(coords)))
  }
})
