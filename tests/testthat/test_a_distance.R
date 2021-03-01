context("distance")

set.seed(1974)
test_matrix <- matrix(rnorm(100), nrow = 10)
index_matrix <- matrix(c(rep(0:9, each = 10), rep(0:9, 10)),
											 ncol = 2, byrow = FALSE)
test_matrix <- t(test_matrix)

test_that("Euclidean distances are correct", {
  distances <- as.matrix(dist(test_matrix, method = "euclidean"))
  new_distances <- distance(as.vector(index_matrix[, 2]),
                            as.vector(index_matrix[, 1]),
                            x = test_matrix,
                            "Euclidean",
                            verbose = FALSE)
  diffs <- as.vector(distances) - new_distances
  expect_lt(sum(diffs), 1e-10)
  expect_equal(attr(new_distances, "method"), "euclidean")
})

test_that("Cosine distances are correct", {
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
                            verbose = FALSE)
  diffs <- as.vector(distances) - new_distances
  expect_lt(sum(diffs), 1e-10)
  expect_equal(attr(new_distances, "method"), "cosine")
})