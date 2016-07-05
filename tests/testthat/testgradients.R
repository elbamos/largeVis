context("Gradients")

dGrad <- deriv(~ sqrt((x_i - y_i)^2 + (x_j - y_j)^2),
               namevec = c("x_i", "y_i", "x_j", "y_j"),
               function.arg = c("x_i", "x_j", "y_i", "y_j"))

posGrad <- deriv(~ log(1 / (1 + (alpha * d^f))),
               namevec = c("d"),
               function.arg = c("d", "f", "alpha"))

negGrad <- deriv(~ g * log(1 - (1 / (1 + (alpha * d^f)))),
               namevec = c("d"),
               function.arg = c("d", "f", "alpha", "g"))

posGradE <- deriv(~ log(1 / (1 + exp(d^f))),
               namevec = c("d"),
               function.arg = c("d", "f"))

negGradE <- deriv(~ g * log(1 - (1 / (1 + exp(d^f)))),
               namevec = c("d"),
               function.arg = c("d", "f", "g"))

grad <- function(x, y, gfunc, f, ...) {
  d <- dGrad(x[1], x[2], y[1], y[2])
  g <- gfunc(d, f, ...)
  return(as.numeric(attr(g, 'gradient')) * as.numeric(attr(d, 'gradient')))
}


ntests <- 2

test_that("Positive Gradient f = 2", {
  set.seed(1972)
  for (i in 1:ntests) {
    f <- 2
    alpha <- runif(1, 1e-5, 10)
    x_i <- rnorm(2)
    y_i <- rnorm(2)
    suppressWarnings(grads <- grad(x_i, y_i, posGrad, f , alpha = alpha))
    cgrads <- testPositiveGradient(x_i, y_i, alpha, f) * 2
    expect_equal(grads[1], cgrads[1], info = paste("f=", f, "a=", alpha))
    expect_equal(grads[3], cgrads[2], info = paste("f=", f, "a=", alpha))
    expect_equal(grads[2], - cgrads[1], info = paste("f=", f, "a=", alpha))
    expect_equal(grads[4], - cgrads[2], info = paste("f=", f, "a=", alpha))
  }
})

test_that("Negative Gradient f = 2", {
  set.seed(1972)
  f <- 2
  for (i in 1:ntests) {
    alpha <- runif(1, 1e-5, 10)
    g <- runif(1, 1, 10)
    x_i <- rnorm(2) * 100
    y_i <- rnorm(2) * 100
    suppressWarnings(grads <- grad(x_i, y_i, gfunc = negGrad, f , alpha = alpha, g = g))
    cgrads <- testNegativeGradient(x_i, y_i, alpha, g, f) * 2
    expect_equal(grads[1], cgrads[1], info = paste("f=", f, "a=", alpha))
    expect_equal(grads[3], cgrads[2], info = paste("f=", f, "a=", alpha))
    expect_equal(grads[2], - cgrads[1], info = paste("f=", f, "a=", alpha))
    expect_equal(grads[4], - cgrads[2], info = paste("f=", f, "a=", alpha))
  }
})

test_that("Positive Gradient E f = 2", {
  set.seed(1972)
  skip("Positive E gradient test disabled")
  f <- 2
  for (i in 1:ntests) {
    alpha <- runif(1, 1e-5, 10)
    x_i <- rnorm(2)
    y_i <- rnorm(2)
    suppressWarnings(grads <- grad(x_i, y_i, posGradE, f ))
    cgrads <- testPositiveGradient(x_i, y_i, 0, f) * 2
    expect_equal(grads[1], cgrads[1], info = paste("f=", f))
    expect_equal(grads[3], cgrads[2], info = paste("f=", f))
    expect_equal(grads[2], - cgrads[1], info = paste("f=", f))
    expect_equal(grads[4], - cgrads[2], info = paste("f=", f))
  }
})

test_that("Negative Gradient E f = 2", {
  set.seed(1972)
  f <- 2
  for (i in 1:ntests)  {
    g <- runif(1, 1, 10)
    x_i <- rnorm(2) * 5
    y_i <- rnorm(2) * 5
    suppressWarnings(grads <- grad(x_i, y_i, gfunc = negGradE, f, g = g))
    grads <- pmax(pmin(grads, 1), -1) * 2
    cgrads <- testNegativeGradient(x_i, y_i, alpha = 0, g, f)
    expect_equal(grads[1], cgrads[1], info = paste("f=", f))
    expect_equal(grads[3], cgrads[2], info = paste("f=", f))
    expect_equal(grads[2], - cgrads[1], info = paste("f=", f))
    expect_equal(grads[4], - cgrads[2], info = paste("f=", f))
  }
})
