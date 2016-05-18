#' Find approximate k-Nearest Neighbors using random projection tree search.
#'
#' This is a very fast and accurate algorithm for finding k-nearest neighbors.
#'
#' @param x A matrix
#' @param K How many nearest neighbors to seek for each node
#' @param n.trees The number of trees to build
#' @param tree.threshold The threshold for creating a new branch.  A value of 1.5 * K guarantees that
#' the algorithm will return K candidate nearest neighbors for each vertex. The paper authors suggested
#' using a value equivalent to the number of features in the input set.  However this can lead to very lengthy
#' calculation times.
#' @param max.iter Number of iterations in the neighborhood exploration phase
#' @param verbose Whether to print verbose logging using the \code{progress} package
#'
#' @return A [N, K] integer matrix showing the estimated K nearest neighbors for each vertex.
#' @export
#'
#' @examples
#' \dontrun{
#' library(kknn) # Shamelessly borrowing sample data from another package
#' data(miete)
#' miete <- model.matrix(~ ., miete)
#' system.time(neighbors <- randomProjectionTreeSearch(miete))
#' }
#'
randomProjectionTreeSearch <- function(x,
                                       K = 5, #
                                       n.trees = 2, # how many trees to build
                                       tree.threshold =  max(10, ncol(x)), # the maximum number of nodes per leaf
                                       max.iter = 2, # in the neighborhood exploration phase, the number of iterations
                                       verbose= TRUE) {
  N <- nrow(x)
  if (tree.threshold < 10) stop("The tree threshold must be at least 10.")
  if (any(is.nan(x))) stop("NaNs detected in x.")
  if (any(is.na(x))) stop("NAs detected in x.")
  if (any(is.infinite(x))) stop("Infs detected in x.")



  # random projection trees
  if (verbose[1]) ptick <- progress::progress_bar$new(total = (n.trees * N) + (N * (max.iter + 1)),
                                                      format = 'Random projection trees [:bar] :percent/:elapsed eta: :eta', clear = FALSE)$tick
  else ptick <- function(ticks) {}
  ptick(0)
  knns <- searchTrees(tree.threshold, n.trees, K, 32, max.iter, x, callback = ptick)

  if (sum(colSums(knns != -1) == 0) + sum(is.na(knns)) + sum(is.nan(knns)) > 0)
    stop("After neighbor search, no candidates for some nodes.")
  if (verbose[1] && sum(knns == -1) > 0)
    warning("Wanted to find", nrow(knns) * ncol(knns), " neighbors, but only found",
                  ((nrow(knns) * ncol(knns)) - sum(knns == -1)))

  return(knns)
}
