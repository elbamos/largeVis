context("largeVis")

test_that("Can determine iris neighbors", {
  data (iris)
  dat <- as.matrix(iris[,1:4])
  dat <- scale(dat)
  dupes = which(duplicated(dat))
  dat <- dat[-dupes,]
  neighbors <- randomProjectionTreeSearch(dat, K = 5, n.trees = 2, tree.threshold = 20, max.iter = 10,
                                          verbose = F)
  expect_equal(nrow(neighbors), 5)
  expect_equal(ncol(neighbors), nrow(dat))
  expect_equal(sum(neighbors == -1),0)
  expect_equal(sum(neighbors[1:40] > 50), 0)
})

test_that("largeVis works", {
  data(iris)
  dat <- as.matrix(iris[,1:4])
  dat <- scale(dat)
  dupes = which(duplicated(dat))
  dat <- dat[-dupes,] # duplicated data potentially can cause the algorithm to fail
  visObject <- largeVis(dat, pca.first = F,
                        max.iter = 20, sgd.batches = 800000,
                        K = 10,  gamma = 2, rho = 1, M = 40, alpha = 20,verbose=F)
  expect_equal(sum(any(is.na(visObject$coords)) + any(is.nan(visObject$coords)) + any(is.infinite(visObject$coords))), 0)
})
