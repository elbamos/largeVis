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
#' @param verbose Verbosity.
#'
#' @return A vector of the distances between the columns in `x` indexed by `i` and `j`, with attribute \code{method} giving the \code{distance_method}.
#' @export
distance <- function(x,
                     i,
                     j,
                     distance_method,
                     verbose) UseMethod("distance")

#' @export
#' @rdname distance
distance.matrix <- function(x,
                     i,
                     j,
                     distance_method = "Euclidean",
                     verbose = getOption("verbose", TRUE)) {
  ret = fastDistance(i,
                       j,
                       x,
                       distance_method,
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
                                   verbose = getOption("verbose", TRUE)) {
  ret <- fastCDistance(i,
                       j,
                       x@i,
                       x@p,
                       x@x,
                       distance_method,
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
                                  verbose = getOption("verbose", TRUE)) {
  ret <- fastSDistance(i,
                       j,
                       x@i,
                       x@j,
                       x@p,
                       distance_method,
                       verbose)
  attr(ret, "method") <- tolower(distance_method)
  ret
}
