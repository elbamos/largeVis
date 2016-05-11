#' Find approximate k-Nearest Neighbors using random projection tree search.
#'
#' This is a very fast and accurate algorithm for finding k-nearest neighbors. It is implemented here in C++ and parallel execution
#' is enabled using `parallel::mclapply`.
#'
#' @param x A matrix
#' @param K How many nearest neighbors to seek for each node
#' @param n.trees The number of trees to build
#' @param tree.threshold The threshold for creating a new branch
#' @param max.iter Number of iterations in the neighborhood exploration phase
#' @param verbose Whether to print verbose logging using the \code{progress} package
#'
#' @return A [nrow(x), K] integer matrix showing the estimated K nearest neighbors for each vertex.
#' @export
#'
#' @examples
#'
randomProjectionTreeSearch <- function(x,
                                       K = 5, #
                                       n.trees = 2, # how many trees to build
                                       tree.threshold = max(10, K * 2), # the maximum number of nodes per leaf
                                       max.iter = 2, # in the neighborhood exploration phase, the number of iterations
                                       verbose= TRUE) {
  N <- nrow(x)
  if (tree.threshold < 10) stop("The tree threshold must be at least 7.")
  # random projection trees
  tree_assignments <- list()
  if (verbose[1]) ptick <- progress::progress_bar$new(total = n.trees * N, format = 'Exploring random projection trees [:bar] :percent/:elapsed eta: :eta', clear = FALSE)$tick
  else ptick <- function(ticks) {}
  ptick(0)

  knns <- parallel::mclapply(1:n.trees, FUN=function(T) {
    someknns <- matrix(0, ncol = N, nrow = tree.threshold)
    searchTree(tree.threshold, 1:N, x, someknns, ptick)
    someknns
  })
  knns <- do.call(rbind, knns)

  if (sum(colSums(knns) == 0) > 0) stop("Random projection trees found no candidates for some nodes.")

  if (verbose[1]) ptick <- progress::progress_bar$new(total = max.iter * N, format = 'Exploring Neighbors [:bar] :percent/:elapsed eta: :eta', clear = FALSE)$tick
  else ptick <- function(ticks) {}

  outputKnns <- matrix(0, nrow = K, ncol = N)
  neighbors_inner(max.iter, knns, x, outputKnns, ptick)
  return(outputKnns)
}
