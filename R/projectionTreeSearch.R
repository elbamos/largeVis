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
                                       tree.threshold = k * 2, # the maximum number of nodes per leaf
                                       distance.method = 'euclidean',
                                       max.iter = 2, # in the neighborhood exploration phase, the number of iterations
                                       search.harder = TRUE, # whether to check 3d-degree neighbors if insufficient neighbors were found
                                       verbose= TRUE) {
  N <- nrow(x)
  # random projection trees
  tree_assignments <- list()
  if (verbose) cat("Creating random projection trees...")

  tree_assignments <- lapply(1:n.trees, FUN = function(T) {
    tree <- partition(indices = 1:N, .threshold = tree.threshold, .data = x)
    knns <- matrix(0, nrow = tree.threshold, ncol = N)
    for (leaf in tree) {
      knns[1:length(leaf),leaf] <- rep(leaf, length(leaf))
    }
    knns
  })
  new_knns <- do.call(rbind, tree_assignments)
  if (verbose) cat("...done\n")

  if (verbose) progress <- progress::progress_bar$new(total = max.iter, format = 'Exploring Neighbors [:bar] :percent eta: :eta')#,
                           #       utils::txtProgressBar(min = 0, max = sgd.batches, style = 3))
  for (T in 1:max.iter) {
    if (verbose) progress$tick()
    knns <- new_knns
    new_knns <- foreach(i = 1:N, .combine = cbind, .multicombine = TRUE) %do% {
      candidates <- c(knns[,i], knns[,knns[,i]])
      candidates <- candidates[candidates != i]
      candidates <- candidates[! candidates == 0]
      candidates <- unique(candidates)
      if (length(candidates) <= K) {
        if (search.harder) {
          candidates <- c(candidates, knns[,knns[,knns[,i]]])
          candidates <- candidates[candidates != i]
          candidates <- candidates[! candidates == 0]
          candidates <- unique(candidates)
        } else {
          stop(paste("Attempted to visit neighbors but the number of knn candidates was smaller than K."))
        }
      }
      distances <- proxy::dist(x = x[i,,drop=FALSE],
                               y = x[candidates,],
                               method = distance.method)
      candidates <- candidates[order(distances,decreasing=FALSE)[1:min(K, length(candidates))]]
      if (length(candidates) < K) candidates <- c(candidates,(rep(0, K - length(candidates))))
      candidates
    }
  }
  return(new_knns)
}
