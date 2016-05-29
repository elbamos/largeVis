
#' Build an edge-weight matrix for the LargeVis algorithm.
#'
#' @param x A sparseMatrix, either a \code{\link[Matrix]{CsparseMatrix-class}} or \code{\link[Matrix]{TsparseMatrix-class}}
#' @param i Indices of one node of the nearest-neighbor graph.
#' @param j Indices of the other node.
#' @param p Integer vector of pointers to the initial index of elements of each column. See \code{\link[Matrix]{CsparseMatrix-class}}.
#' @param d The distances between the nodes identified in parameters \code{i} and \code{j}.
#' @param perplexity See the paper for discussion.
#' @param verbose Verbosity
#'
#' @details Implements the portion of the LargeVis algorithm that converts distances between nearest neighbors to an
#' edge-weight graph.
#'
#' @importFrom stats optimize
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @return A list containing: \describe{
#' \item{"sigmas"}{A vector of \eqn{2 \dot \sigma^2} calculated for each node.}
#' \item{"wij"}{A symmetric, sparse matrix of the weights for each edge between nearest neighbors.}
#' }
#' @importClassesFrom Matrix CsparseMatrix
#' @importClassesFrom Matrix TsparseMatrix
#' @export
buildEdgeMatrix <- function(x,
                            i,
                            j,
                            p,
                            d,
                            perplexity,
                            verbose) UseMethod("buildEdgeMatrix")


#' @export
#' @rdname buildEdgeMatrix
buildEdgeMatrix.default <- function(x = NULL,
                                    i,
                                    j,
                                    p = NULL,
                                    d,
                                    perplexity = 50,
                                    verbose = TRUE) {
  if (is.null(p)) p <- i2p(i)
  N <- max(max(i), max(j)) + 1

  if (verbose) {
    progress <- txtProgressBar(max = N, title = "sigmas")
    cat("Estimating sigmas\n")
  }

  perplexity <- log2(perplexity)
  sigmas <- parallel::mclapply(1:N, FUN = function(idx) { # nocov start
    if (verbose) setTxtProgressBar(progress, idx)
    x_i <- d[(p[idx] + 1):(p[idx + 1])]
    ret <- optimize(f = sigFunc,
                    x = x_i,
                    perplexity = perplexity,
                    interval = c(0, 10000))
  }) # nocov end
  sigmas <- sapply(sigmas, `[[`, 1)

  if (verbose) close(progress)

  if (any(is.na(sigmas)) +
      any(is.infinite(sigmas)) +
      any(is.nan(sigmas)) +
      any( (sigmas == 0)) > 0)
    stop("An error has propogated into the sigma vector.")

  if (length(sigmas) != N) stop("Wrong sigma count")

  if (! requireNamespace("Matrix", quietly = T))
    stop("The Matrix package must be available.")

  if (verbose) cat("Calculating w_{ij}.\n")

  wij <- distMatrixTowij(i, j, d, sigmas, N, verbose)

  if (any(is.na(wij@x)) +
      any(is.infinite(wij@x)) +
      any(is.nan(wij@x)) +
      any( (wij@x == 0)) > 0)
    stop(paste("An error has propogated into the w_{ij} vector.",
              "This probably means the input data wasn't scaled."))

  return(list(sigmas = sigmas, wij = wij))
}

#' @export
#' @rdname buildEdgeMatrix
buildEdgeMatrix.CsparseMatrix <- function(x,
                                          i = NULL,
                                          j = NULL,
                                          p = NULL,
                                          d = NULL,
                                          perplexity = 50,
                                          verbose = TRUE) {
  # Will have x@i, which is quickly varying, and x@p, and x@x
  is <- rep(0:(nrow(x) - 1), diff(x@p))
  js <- x@i
  ps <- x@p
  ds <- x@x
  NextMethod("buildEdgeMatrix",
             i = is,
             j = js,
             p = ps,
             d = ds,
             perplexity = perplexity,
             verbose = verbose)
}

#' @inheritParams buildEdgeMatrix.CsparseMatrix
#' @export
#' @rdname buildEdgeMatrix

buildEdgeMatrix.TsparseMatrix <- function(x,
                                      i = NULL,
                                      j = NULL,
                                      p = NULL,
                                      d = NULL,
                                      perplexity = 50,
                                      verbose = TRUE) {
  ps <- i2p(x@i)
  is <- x@i
  js <- x@j
  ds <- x@x
  NextMethod("buildEdgeMatrix",
             i = is,
             j = js,
             p = ps,
             d = ds,
             perplexity = perplexity,
             verbose = verbose)
}

i2p <- function(is) {
  N <- max(is)
  ps <- rep(NA, N + 1)
  diffs <- diff(is)
  ps[is[which(diffs > 0)] + 2] <- which(diffs > 0) + 1
  good <- cumsum(!is.na(ps))
  ps <- ps[good + 1]
  ps[1] <- 1
  ps <- ps - 1
  ps[length(ps) + 1] <- length(is)
  return(ps)
}
