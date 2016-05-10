#' Find approximate k-Nearest Neighbors using random projection tree search
#'
#' @param x A matrix
#' @param K How many nearest neighbors to seek for each node
#' @param n.trees The number of trees to build
#' @param tree.threshold The threshold for creating a new branch
#' @param distance.method Method used to calculate distance between vertices, passed to \code{proxy::dist()}.
#' @param max.iter Number of iterations in the neighborhood exploration phase
#' @param search.harder Whether to check 3d-degree neighbors if \code{< K} local neighbors were found
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
                                       tree.threshold = K * 2, # the maximum number of nodes per leaf
                                       distance.method = 'euclidean',
                                       max.iter = 2, # in the neighborhood exploration phase, the number of iterations
                                       search.harder = TRUE, # whether to check 3d-degree neighbors if insufficient neighbors were found
                                       verbose= TRUE) {
  N <- nrow(x)
  # random projection trees
  tree_assignments <- list()
  if (verbose) cat("Creating random projection trees...")

  tree_assignments <- parallel::mclapply(1:n.trees, FUN = function(T) {
    tree <- partition(indices = 1:N, .threshold = tree.threshold, .data = x)
    knns <- matrix(0, nrow = tree.threshold, ncol = N)
    for (leaf in tree) {
      knns[1:length(leaf),leaf] <- rep(leaf, length(leaf))
    }
    knns
  })
  new_knns <- do.call(rbind, tree_assignments)
  if (verbose[1]) cat("...done\n")

  if (verbose[1]) ptick <- progress::progress_bar$new(total = max.iter * N, format = 'Exploring Neighbors [:bar] :percent eta: :eta', clear = FALSE)$tick
  else ptick <- function(ticks) {}

  for (T in 1:max.iter) {
    old_knns <- new_knns
    new_knns <- neighbors_inner(K, old_knns, x, ptick)
  }

  if (verbose[1]) print('...done!')
  return(new_knns)
}
