#' Calculate pairwise Euclidean or angular distances efficiently
#'
#' This function is a wrapper around a C++ function that calculates pairwise distances in a memory- and CPU-efficient manner.
#'
#' The Euclidean or angular distances between columns in `x` identified by parameters `i` and `j` are calculated and returned.
#'
#' @param i 0-indexed vector of column indices.
#' @param j 0-indexed vector of column indices.
#' @param x A (potentially sparse) matrix, where examples are columns and features are rows.
#' @param distance_method One of "Euclidean" or "Cosine."
#' @param threads The maximum number of threads to spawn. Determined automatically if \code{NULL} (the default).
#' @param verbose Verbosity.
#'
#' @return A vector of the distances between the columns in `x` indexed by `i` and `j`, with attribute \code{method} giving the \code{distance_method}.
#' @family lowmem
#' @export
distance <- function(x,
                     i,
                     j,
                     distance_method,
										 threads = NULL,
                     verbose) UseMethod("distance")

#' @export
#' @rdname distance
distance.matrix <- function(x,
                     i,
                     j,
                     distance_method = "Euclidean",
										 threads = NULL,
                     verbose = getOption("verbose", TRUE)) {
  ret = fastDistance(i,
                       j,
                       x,
                       distance_method,
  										 threads,
                       verbose)
  attr(ret, "method") <- tolower(distance_method)
  ret
}

#' @export
#' @rdname distance
distance.CsparseMatrix <- function(x,
                                   i,
                                   j,
                                   distance_method = "Euclidean",
																	 threads = NULL,
                                   verbose = getOption("verbose", TRUE)) {
  ret <- fastCDistance(i,
                       j,
                       x@i,
                       x@p,
                       x@x,
                       distance_method,
  										 threads,
                       verbose)
  attr(ret, "method") <- tolower(distance_method)
  ret
}

#' @export
#' @rdname distance
distance.TsparseMatrix <- function(
                                  x,
                                  i,
                                  j,
                                  distance_method="Euclidean",
                                  threads = NULL,
                                  verbose = getOption("verbose", TRUE)) {
  ret <- fastSDistance(i,
                       j,
                       x@i,
                       x@j,
                       x@p,
                       distance_method,
  										 threads,
                       verbose)
  attr(ret, "method") <- tolower(distance_method)
  ret
}

#' A utility function to convert a k-NN graph to a pair of 0-indexed vectors of indices.
#'
#' In the returned list, the nodes indexed by `j` are the identified nearest neighbors of the nodes indexed by `i`.
#' In other words, if `i = c(0,0,0,1,1,1)` and `j = c(1,2,3,2,3,4)`, then nodes `1, 2 & 3` are nearest neighbors of node 0,
#' but node 0 is not a nearest neighbor of node 1.
#'
#' @param x A `[K,N]` matrix of indices of the nearest neighbors of each vertex. 0-indexed.
#'
#' @return A list with fields:
#' \describe{
#' \item{i}{The slowly-varying indices of x}
#' \item{j}{The quickly-varying indices of x}
#' }
#' @export
neighborsToVectors <- function(x) {
  K <- nrow(x)
  N <- ncol(x)
  is <- rep(0:(N - 1), each = K)
  js <- as.vector(x)
  is <- is[! js == -1]
  js <- js[! js == -1]
  return (list(i = is, j = js))
}
