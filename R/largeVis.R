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
#' @param check.assumptions Whether to check the input matrix for duplicates, \code{NA}`s, etc.
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
#' @param weight_pos_samples See \code{\link{projectKNNs}}.
#' @param alpha See \code{\link{projectKNNs}}.
#' @param gamma See \code{\link{projectKNNs}}.
#' @param rho See \code{\link{projectKNNs}}.
#' @param min_rho \code{\link{projectKNNs}}.
#' @param coords A [N,K] matrix of coordinates to use as a starting point -- useful for refining an embedding in stages.
#' @param verbose Verbosity
#' @param ... See paper
#'
#' @return A `largeVis` object with the following slots:
#'  \describe{
#'    \item{'knns'}{An [N,K] integer matrix, which is an adjacency list of each vertex' identified nearest neighbors.
#'    If the algorithm failed to find \code{K} neighbors, the matrix is padded with \code{NA}'s.}
#'    \item{'wij'}{A sparse [N,N] matrix where each cell represents \eqn{w_{ij}}.}
#'    \item{'call'}{The call.}
#'    \item{'coords'}{A [N,D] matrix of the embedding of the dataset in the low-dimensional space.}
#'  }
#'
#'
#' @export
#' @importFrom stats optimize princomp
#' @importFrom utils setTxtProgressBar txtProgressBar
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

                     check.assumptions = TRUE,

                     n_trees = 50,
                     tree_threshold = max(10, nrow(x)),
                     max_iter = 3,
                     max_depth = 32,
                     distance_method = "Euclidean",

                     perplexity = 50,

                     sgd_batches = ncol(x) * 20000,
                     M = 5,
                     weight_pos_samples = TRUE,
                     alpha = 1,
                     gamma = 7,
                     rho = 1,
                     min_rho = 0,

                     coords = NULL,

                     verbose = TRUE,
                    ...) {
  N <- ncol(x)

  if (check.assumptions)   {
    if ( (any(is.na(x)) +
         any(is.infinite(x)) +
         any(is.nan(x))) > 0)
      stop("Missing values present in input matrix.")
  }

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
  # These vectors are analogous to the components of a sparse matrix,
  # but both triple and C-compressed forms are created.
  # The i and j vectors are 0-indexed while p is 1-indexed.
  is <- rep(0:(N - 1), each = K)
  js <- as.vector(knns)
  is <- is[! js == -1]
  js <- js[! js == -1]
  dupes <- duplicated(data.frame(is, js))
  is <- is[! dupes]
  js <- js[! dupes]
  ord <- order(is)
  is <- is[ord]
  js <- js[ord]

  #######################################################
  # Calculate edge weights for candidate neighbors
  #######################################################
  if (verbose) cat("Calculating neighbor distances.\n")

  xs <- distance(is, js, x, distance_method,verbose)[, 1]

  if (verbose) cat("\n")

  if ( (any(is.na(xs)) + any(is.infinite(xs)) + any(is.nan(xs)) + any(xs == 0)) > 0)
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

  ps <- i2p(is)
  sigwij <- buildEdgeMatrix(i = is,
                         j = js,
                         p = ps,
                         d = xs,
                         perplexity = perplexity,
                         verbose = verbose)


  #######################################################
  # Estimate embeddings
  #######################################################
  coords <- projectKNNs(wij = sigwij$wij,
                        dim = dim,
                        sgd_batches = sgd_batches,
                        M = M,
                        weight_pos_samples = weight_pos_samples,
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
  knns[knns == -1] <- NA

  returnvalue <- list(
    knns = t(knns),
    wij = sigwij$wij,
    call = sys.call(),
    coords = coords,
    sigmas = sqrt(sigwij$sigmas / 2)
  )

  class(returnvalue) <- "largeVis"
  return(returnvalue)
}

##########################################
# Some helper functions useful for debugging
##########################################

pji <- function(x_i, sigma)  exp(- (x_i^2) / sigma) / sum(exp(- (x_i^2) / sigma))
perp <- function(x_i, sigma) - sum(log2(pji(x_i, sigma))) / length(x_i)
pdiff <- function(x_i, sigma, perplexity) (perplexity - perp(x_i, sigma))^2
