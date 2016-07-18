#' Apply the LargeVis algorithm for visualizing large high-dimensional datasets.
#'
#' Implements the \code{vis}
#'
#' Note that this implementation expects the data to be free of \code{NaN}'s, \code{NA}'s, \code{Inf}'s, and duplicate rows.
#' If any of these assumptions are violated, the algorithm will fail. It is also usually a good idea to scale the input data
#' to have unit norm and mean 0. If there are large values in the input matrix, some computations may oveflow.
#'
#' @param x A matrix, where the features are rows and the examples are columns.
#' @param dim The number of dimensions in the output
#' @param K The number of nearest-neighbors to use in computing the kNN graph
#' @param n_trees See \code{\link{randomProjectionTreeSearch}}.  The default is set at 50, which is the number
#' used in the examples in the original paper.
#' @param tree_threshold See \code{\link{randomProjectionTreeSearch}}.  By default, this is the number of features
#' in the input set, which is the setting used in the examples in the original paper.  Note the time and memory requirements:
#' the first pass through the neighborhood exploration phases will involve up to \eqn{N * nTrees * threshold} comparisons.
#' @param max_depth See \code{\link{randomProjectionTreeSearch}}
#' @param max_iter See \code{\link{randomProjectionTreeSearch}}.
#' @param distance_method One of "Euclidean" or "Cosine."  See \code{\link{randomProjectionTreeSearch}}.
#' @param perplexity See paper
#' @param sgd_batches See \code{\link{projectKNNs}}.
#' @param M See \code{\link{projectKNNs}}.
#' @param alpha See \code{\link{projectKNNs}}.
#' @param gamma See \code{\link{projectKNNs}}.
#' @param rho See \code{\link{projectKNNs}}.
#' @param min_rho \code{\link{projectKNNs}}.
#' @param save_neighbors Whether to include in the output the adjacency matrix of nearest neighbors.
#' @param save_sigmas Whether to include in the output the esimates values of sigma.
#' @param coords A [N,K] matrix of coordinates to use as a starting point -- useful for refining an embedding in stages.
#' @param verbose Verbosity
#' @param ... See paper
#'
#' @return A `largeVis` object with the following slots:
#'  \describe{
#'    \item{'knns'}{An [N,K] 0-indexed integer matrix, which is an adjacency list of each vertex' identified nearest neighbors.
#'    If the algorithm failed to find \code{K} neighbors, the matrix is padded with \code{NA}'s.}
#'    \item{'wij'}{A sparse [N,N] matrix where each cell represents \eqn{w_{ij}}.}
#'    \item{'call'}{The call.}
#'    \item{'coords'}{A [N,D] matrix of the embedding of the dataset in the low-dimensional space.}
#'    \item{'sigmas'}{A [N] vector of the values of sigma estimated for each vertex. Primarily useful for debugging
#'    purposes and therefore not returned by default.}
#'  }
#'
#' @export
#' @references Jian Tang, Jingzhou Liu, Ming Zhang, Qiaozhu Mei. \href{https://arxiv.org/abs/1602.00370}{Visualizing Large-scale and High-dimensional Data.}
#'
#' @examples
#' # iris
#' data(iris)
#' dat <- as.matrix(iris[,1:4])
#' dat <- scale(dat)
#' dupes = which(duplicated(dat))
#' dat <- dat[-dupes,] # duplicated data potentially can cause the algorithm to fail
#' dat <- t(dat)
#' visObject <- vis(dat, max_iter = 20, sgd_batches = 800000,
#'                      K = 10,  gamma = 2, rho = 1, M = 40, alpha = 20,verbose=FALSE)
#'\dontrun{
#' # mnist
#' load("./mnist.Rda")
#' dat <- mnist$images
#' dim(dat) <- c(42000, 28 * 28)
#' dat <- (dat / 255) - 0.5
#' dat <- t(dat)
#' coords <- vis(dat, check=FALSE,
#'              n_tree = 50, tree_th = 200,
#'              K = 50, alpha = 2, max.iter = 4)
#' }
#'
vis <- function(x,
                     dim = 2,
                     K = 40,

                     n_trees = 50,
                     tree_threshold = max(10, ncol(x)),
                     max_iter = 1,
                     max_depth = 32,
                     distance_method = "Euclidean",

                     perplexity = 50,

                     sgd_batches = NULL,
                     M = 5,
                     alpha = 1,
                     gamma = 7,
                     rho = 1,
                     min_rho = 0,

                     coords = NULL,

                     save_neighbors = TRUE,
                     save_sigmas = FALSE,

                     verbose = TRUE,
                    ...) {

  #############################################
  # Search for kNearestNeighbors
  #############################################
  knns <- randomProjectionTreeSearch(x,
                                     n_trees = n_trees,
                                     tree_threshold = tree_threshold,
                                     K = K,
                                     max_iter = max_iter,
                                     max_depth = max_depth,
                                     distance_method = distance_method,
                                     verbose = verbose)
  #############################################
  # Clean knns
  #############################################
  if (verbose[1]) cat("Calculating edge weights...")
  neighbor_indices <- neighborsToVectors(knns)
  if (! save_neighbors) rm(knns)
  gc()

  #######################################################
  # Calculate edge weights for candidate neighbors
  #######################################################
  if (verbose) cat("Calculating neighbor distances.\n")

  xs <- distance(x = x,
                 neighbor_indices$i,
                 neighbor_indices$j,
                 distance_method,
                 verbose)[, 1]
  if (distance_method == "Euclidean") xs <- xs^2

  if (verbose) cat("\n")

  if ( (any(is.na(xs)) +
        any(is.infinite(xs)) +
        any(is.nan(xs)) +
        any(xs == 0)) > 0)
  stop("An error leaked into the distance calculation - check for duplicates")
  if (any(xs > 27)) { # nocov start
    warning(paste(
    "The Distances between some neighbors are large enough to cause the calculation of p_{j|i} to overflow.",
    "Scaling the distance vector."))
    xs <- scale(xs, center = FALSE)
  } # nocov end

  #######################################################
  # Get w_{ij}
  #######################################################

  # sigwij <- buildEdgeMatrix(i = neighbor_indices$i,
  #                        j = neighbor_indices$j,
  #                        d = xs,
  #                        perplexity = perplexity,
  #                        verbose = verbose)

  wij <- referenceWij(i = neighbor_indices$i,
                        j = neighbor_indices$j,
                        d = xs,
                        perplexity = perplexity)
  sigwij <- list(wij = wij, sigmas = rep(0, ncol(x)))
  rm(neighbor_indices, xs)
  if (! save_sigmas) sigwij$sigmas <- NULL
  gc()
  #######################################################
  # Estimate embeddings
  #######################################################
  coords <- projectKNNs(wij = sigwij$wij,
                        dim = dim,
                        sgd_batches = sgd_batches,
                        M = M,
                        gamma = gamma,
                        verbose = verbose,
                        alpha = alpha,
                        coords = coords,
                        rho = rho,
                        min_rho = min_rho,
                        ...)

  #######################################################
  # Cleanup
  #######################################################

  returnvalue <- list(
    knns = t(knns),
    wij = sigwij$wij,
    call = sys.call(),
    coords = coords
  )

  if (save_neighbors) {
    knns[knns == -1] <- NA
    returnvalue$knns <- t(knns)
  }
  if (save_sigmas) {
    sigmas <- sqrt(sigwij$sigmas / 2)
    returnvalue$sigmas <- sigmas
  }

  class(returnvalue) <- "largeVis"
  return(returnvalue)
}

##########################################
# Some helper functions useful for debugging
##########################################

pji <- function(x_i, sigma)
  exp(- (x_i ^ 2) / sigma) / sum(exp( - (x_i ^ 2) / sigma))
perp <- function(x_i, sigma)
  - sum(log2(pji(x_i, sigma))) / length(x_i)
pdiff <- function(x_i, sigma, perplexity)
  (perplexity - perp(x_i, sigma))^2
