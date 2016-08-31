skip_old_windows <- function() {
#	testthat::skip_if_not(R.Version()$arch != "i386", "largeVis does not run on 32-bit Windows.")
}
context("distance")

set.seed(1974)
test_matrix <- matrix(rnorm(100), nrow = 10)
index_matrix <- matrix(c(rep(0:9, each = 10), rep(0:9, 10)),
                       ncol = 2, byrow = FALSE)
test_matrix <- t(test_matrix)

test_that("Euclidean distances are correct", {
	skip_old_windows()
  distances <- as.matrix(dist(test_matrix, method = "euclidean"))
  new_distances <- distance(as.vector(index_matrix[, 2]),
                            as.vector(index_matrix[, 1]),
                            x = test_matrix,
                            "Euclidean",
  													threads = 2,
                            verbose = FALSE)
  diffs <- as.vector(distances) - new_distances
  expect_lt(sum(diffs), 1e-10)
})

test_that("Cosine distances are correct", {
	skip_old_windows()
  set.seed(1974)
  cos.sim <- function(x, i, j) {
    A <- x[, i]
    B <- x[, j]
    return( 1 - (sum(A * B) / sqrt(sum(A ^ 2) * sum(B ^ 2)) ))
  }
  distances <- apply(index_matrix + 1,
                     MARGIN=1,
                     FUN = function(x) cos.sim(test_matrix, x[1], x[2]))
  new_distances <- distance(as.vector(index_matrix[, 2]),
                            as.vector(index_matrix[, 1]),
                            x = test_matrix,
                            "Cosine",
  													threads = 2,
                            verbose = FALSE)
  diffs <- as.vector(distances) - new_distances
  expect_lt(sum(diffs), 1e-10)
})