context("wij")

test_that("wij handles small K", {
  set.seed(1974)
  data(iris)
  dat <- as.matrix(iris[, 1:4])
  dat <- scale(dat)
  dupes <- which(duplicated(dat))
  dat <- dat[-dupes, ]
  dat <- t(dat)
  neighbors <- randomProjectionTreeSearch(dat, K = 5, verbose = FALSE)
  edges <- buildEdgeMatrix(dat, neighbors, verbose = FALSE)
  expect_silent(wij <- buildWijMatrix(edges))
})