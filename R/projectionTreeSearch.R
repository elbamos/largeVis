#' Find approximate k-Nearest Neighbors using random projection tree search.
#'
#' A fast and accurate algorithm for finding approximate k-nearest neighbors.
#'
#' Note that the algorithm does not guarantee that it will find K neighbors for each node. A
#' warning will be issued if it finds fewer neighbors than requested. If the input data contains
#' distinct partitionable clusters, try increasing the \code{tree_threshold} to increase the number
#' of returned neighbors.
#'
#' @param x A (potentially sparse) matrix, where examples are columnns and features are rows.
#' @param K How many nearest neighbors to seek for each node.
#' @param n_trees The number of trees to build.
#' @param tree_threshold The threshold for creating a new branch.  The paper authors suggest
#' using a value equivalent to the number of features in the input set.
#' @param max_iter Number of iterations in the neighborhood exploration phase.
#' @param max_depth The maximum level of recursion.
#' @param distance_method One of "Euclidean" or "Cosine."
#' @param verbose Whether to print verbose logging using the \code{progress} package.
#'
#' @return A [K, N] matrix of the approximate K nearest neighbors for each vertex.
#' @export
randomProjectionTreeSearch <- function(x,
                                       K = 150,
                                       n_trees = 50,
                                       tree_threshold =  max(10, nrow(x)),
                                       max_iter = 1,
                                       max_depth = 32,
                                       distance_method = "Euclidean",
                                       verbose= TRUE)
  UseMethod("randomProjectionTreeSearch")

#' @export
#' @rdname randomProjectionTreeSearch
randomProjectionTreeSearch.matrix <- function(x,
                                       K = 150,
                                       n_trees = 50,
                                       tree_threshold =  max(10, nrow(x)),
                                       max_iter = 1,
                                       max_depth = 32,
                                       distance_method = "Euclidean",
                                       verbose= TRUE) {
  if (verbose) cat("Searching for neighbors.\n")

  if (distance_method == "Cosine") x <- x / rowSums(x)
  
  knns <- searchTrees(threshold = tree_threshold,
                      n_trees = n_trees,
                      K = K,
                      max_recursion_degree = max_depth,
                      maxIter = max_iter,
                      data = x,
                      distance_method,
                      verbose = verbose)

  if (sum(colSums(knns != -1) == 0) > 0)
    stop ("After neighbor search, no candidates for some nodes.")
  if (sum(is.na(knns)) + sum(is.nan(knns)) > 0)
    stop ("NAs or nans in neighbor graph.")
  if (verbose[1] && sum(knns == -1) > 0)
    warning ("Wanted to find", nrow(knns) * ncol(knns),
             " neighbors, but only found",
                  ( (nrow(knns) * ncol(knns)) - sum(knns == -1)))

  return(knns)
}

#' @export
#' @rdname randomProjectionTreeSearch
randomProjectionTreeSearch.CsparseMatrix <- function(x,
                                              K = 150,
                                              n_trees = 50,
                                              tree_threshold =  max(10, nrow(x)),
                                              max_iter = 1,
                                              max_depth = 32,
                                              distance_method = "Euclidean",
                                              verbose= TRUE) {
  if (verbose) cat("Searching for neighbors.\n")

  knns <- searchTreesCSparse(threshold = tree_threshold,
                      n_trees = n_trees,
                      K = K, max_recursion_degree = max_depth,
                      maxIter = max_iter,
                      i = x@i,
                      p = x@p,
                      x = x@x,
                      distance_method,
                      verbose = verbose)

  if (sum(colSums(knns != -1) == 0) > 0)
    stop ("After neighbor search, no candidates for some nodes.")
  if (sum(is.na(knns)) + sum(is.nan(knns)) > 0)
    stop ("NAs or nans in neighbor graph.")
  if (verbose[1] && sum(knns == -1) > 0)
    warning ("Wanted to find", nrow(knns) * ncol(knns),
             " neighbors, but only found",
             ( (nrow(knns) * ncol(knns)) - sum(knns == -1)))

  return(knns)
}

#' @export
#' @rdname randomProjectionTreeSearch
randomProjectionTreeSearch.TsparseMatrix <- function(x,
                                                     K = 150,
                                                     n_trees = 50,
                                                     tree_threshold =
                                                       max(10, nrow(x)),
                                                     max_iter = 1,
                                                     max_depth = 32,
                                                     distance_method =
                                                       "Euclidean",
                                                     verbose= TRUE) {
  if (verbose) cat("Searching for neighbors.\n")

  knns <- searchTreesTSparse(threshold = tree_threshold,
                             n_trees = n_trees,
                             K = K, max_recursion_degree = max_depth,
                             maxIter = max_iter,
                             i = x@i,
                             j = x@j,
                             x = x@x,
                             distance_method,
                             verbose = verbose)

  if (sum(colSums(knns != -1) == 0) > 0)
    stop ("After neighbor search, no candidates for some nodes.")
  if (sum(is.na(knns)) + sum(is.nan(knns)) > 0)
    stop ("NAs or nans in neighbor graph.")
  if (verbose[1] && sum(knns == -1) > 0)
    warning ("Wanted to find", nrow(knns) * ncol(knns),
             " neighbors, but only found",
             ( (nrow(knns) * ncol(knns) ) - sum(knns == -1)))

  return(knns)
}
