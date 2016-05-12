#' Apply the LargeVis algorithm for visualizing large high-dimensional datasets.
#'
#' Implements the \code{largeVis} algorithm by Tang et al.
#'
#' \code{largeVis} estimates a low-dimensional embedding for high-dimensional data, where the distance between vertices
#' in the low-dimensional space is proportional to the distance between them in the high-dimensional space. The algorithm
#' works in 4 phases:
#'
#' \itemize{
#' \item  Estimate candidate nearest-neighbors for each vertex by building \code{n.trees} random projection trees.
#' \item  Estimate \code{K} nearest-neighbors for each vertex by visiting each vertex' 2d-degree neighbors (its neighbors' neighbors).
#' This is repeated \code{max.iter} times.  Note that the original paper suggested a \code{max.iter} of 1, however a larger number
#' may be appropriate for some datasets if the algorithm has trouble finding K neighbors for every vertex.
#' \item Estimate \eqn{p_{j|i}}, the conditional probability that each edge found in the previous step is actually to a
#' nearest neighbor of each of its nodes.
#' \item Using stochastic gradient descent, estimate an embedding for each vertex in the low-dimensional space.
#' }
#'
#' Note that this implementation expects the data to be free of \code{NaN}'s, \code{NA}'s, \code{Inf}'s, and duplicate rows.
#' If any of these assumptions are violated, the algorithm will fail. It is also usually a good idea to sale the input data
#' to have unit norm and mean 0. If there are large values in the input matrix, some computations may oveflow.
#'
#' @param x A matrix. Ideally, the columns should be scaled and normalized to avoid the risk of errors caused by overflow.
#' @param dim The number of dimensions in the output
#' @param K The number of nearest-neighbors to use in computing the kNN graph
#' @param check.assumptions Whether to check the input matrix for duplicates, \code{NA}`s, etc.
#' @param pca.first Whether to apply pca first (can speed-up distance calculations)
#' @param pca.dims How many pca dimensions to use
#' @param n.trees See \code{\link{randomProjectionTreeSearch}}.
#' @param tree.threshold See \code{\link{randomProjectionTreeSearch}}.
#' @param max.iter See \code{\link{randomProjectionTreeSearch}}.
#' @param perplexity See paper
#' @param sgd.batches See \code{\link{projectKNNs}}.
#' @param M See \code{\link{projectKNNs}}.
#' @param weight.pos.samples See \code{\link{projectKNNs}}.
#' @param alpha See \code{\link{projectKNNs}}.
#' @param gamma See \code{\link{projectKNNs}}.
#' @param coords A [N,K] matrix of coordinates to use as a starting point -- useful for refining an embedding in stages.
#' @param verbose Verbosity
#' @param ... See paper
#'
#' @return A `largeVis` object with the following slots:
#'  \itemize{
#'    \item{'knns'}{An [N,K] integer matrix, which is an adjacency list of each vertex' identified nearest neighbors.
#'    If the algorithm failed to find \code{K} neighbors, the matrix is padded with \code{NA}'s.}
#'    \item{'wij'}{A sparse [N,N] adjacency matrix where each cell represents \eqn{w_{ij}}.}
#'    \item{'call'}
#'    \item{'coords'}{A [N,D] matrix of the embedding of the dataset in the low-dimensional space.}
#'  }
#'
#'
#' @export
#' @references Jian Tang, Jingzhou Liu, Ming Zhang, Qiaozhu Mei. \href{https://arxiv.org/abs/1602.00370}{Visualizing Large-scale and High-dimensional Data.}
#' @
#' @examples
#'
#' @useDynLib largeVis
#' @importFrom Rcpp sourceCpp
#'
largeVis <- function(x,
                     dim = 2,
                     K = 40, # number of knn edges per vertex

                     check.assumptions = TRUE,
                     pca.first = TRUE, # whether to apply dimensional reduction first
                     pca.dims = 50,

                     n.trees = 2, # in the random projection tree phase, how many trees to build
                     tree.threshold = K * 2, #the maximum number of nodes per leaf
                     max.iter = 2, # in the neighborhood exploration phase, the number of iterations

                     perplexity = max(K, 50), # hyperparameter for calculating p(j|i)

                     sgd.batches = nrow(x) * 10000,
                     M = 5,
                     weight.pos.samples = TRUE,
                     alpha = 1,
                     gamma = 7,
                     rho = 1,
                     min.rho = 0,

                     coords = NULL,

                     verbose = TRUE,
                    ...) {
  N = nrow(x)
  # Handle pca.first
  shrunken.x <- x
  if (pca.first) {
    if (pca.dims >= ncol(x)) stop("Called for dimensional reduction from ", ncol(x), " to ", pca.dims, " using pca.")
    if (verbose[1]) cat("PCA...")
    shrunken.x <- princomp(x, scores = TRUE)$scores[,1:pca.dims]
    shrunken.x <- scale(shrunken.x)
    if(any(is.nan(shrunken.x))) stop("NaNs in the pca-reduced matrix imply features with 0 variance.")
    if (verbose[1]) cat("done\n")
  }

  if (check.assumptions)   {
    if (any(duplicated(shrunken.x))) stop("Duplicates found.")
    if ((any(is.na(shrunken.x)) + any(is.infinite(shrunken.x)) + any(is.nan(shrunken.x)) + any(shrunken.x == 0)) > 0)
      stop("Missing values present in input matrix.")
  }

  #############################################
  # Search for kNearestNeighbors
  #############################################
  knns <- randomProjectionTreeSearch(shrunken.x,
                                     n.trees = n.trees,
                                     tree.threshold = tree.threshold,
                                     K = K,
                                     max.iter = max.iter,
                                     verbose = verbose)

  #############################################
  # Clean knns and calculate edge distances
  #############################################
  if (verbose[1]) cat("Calculating edge weights...")
  # calculate distance for knns
  is <- rep(1:N, each = K)
  js <- as.vector(knns)
  is <- is[! js == 0]
  js <- js[! js == 0]
  # eliminate symmetrical duplicates
  wrongtri <- is < js
  ts <- is[wrongtri]
  is[wrongtri] <- js[wrongtri]
  js[wrongtri] <- ts
  dupes <- duplicated(data.frame(is, js))
  is <- is[! dupes]
  js <- js[! dupes]
  rm(ts, wrongtri, dupes)
  if (verbose[1]) cat("Calculating neighbor distances...")
  xs <- rep(0, length(is)) # pre-allocate
  distance(is, js, xs, shrunken.x)
  if (verbose) cat("done!\n")

  if ((any(is.na(xs)) + any(is.infinite(xs)) + any(is.nan(xs)) + any(xs == 0)) > 0)
    stop("An error leaked into the distance calculation - check for duplicates")

  # assemble edges into graph, as symmetric sparse matrix
  dist_matrix <- Matrix::sparseMatrix(i = is,
                                      j = js,
                                      x = xs,
                                      symmetric = TRUE)

  ########################################################
  # Estimate sigmas
  ########################################################
  # select denominators sigma to match distributions to given perplexity
  if (verbose[1]) ptick <- progress::progress_bar$new(total = N, format = 'Calculate sigmas [:bar] :percent/:elapsed eta: :eta', clear=FALSE)$tick
  else ptick <- function(tick) {}

  p <- dist_matrix@p
  xvec <- dist_matrix@x
  perplexity = log2(perplexity)
  sigmas <- parallel::mclapply(1:N, FUN = function(idx) {
    ptick(1)
    optimize(f = sigFunc,
             idx = idx,
             p = p,
             x = xvec,
             perplexity = perplexity,
             interval = c(0,100))
  })
  sigmas <- sapply(sigmas, `[[`, 1)
  # note that sigmas are actually (2 * (sigmas^2))
  if (any(is.na(sigmas)) + any(is.infinite(sigmas)) + any(is.nan(sigmas)) + any((sigmas == 0)) > 0) {
    stop("Bad sigma")
  }

  #######################################################
  # Calculate w_{ij}
  #######################################################
  wijVector = rep(0, length(xs * 2))
  if (verbose[1]) progress <- progress::progress_bar$new(total = N, format = 'Calculate p_{j|i} and w_{ij} [:bar] :percent/:elapsed eta: :eta', clear=FALSE)$tick
  else progress <- function(tick) {}
  distMatrixTowij(is, js, xs, sigmas, wijVector, N, progress)
  if (any(is.na(xs)) + any(is.infinite(xs)) + any(is.nan(xs)) + any((xs == 0)) > 0) {
    stop("xs")
  }

  #######################################################
  # Estimate embeddings
  #######################################################
  coords <- projectKNNs(i = is, j = js, x = wijVector,
                        dim = dim,
                        sgd.batches = sgd.batches,
                        M = M,
                        weight.pos.samples = weight.pos.samples,
                        gamma = gamma,
                        verbose = verbose,
                        alpha = alpha,
                        .coords = coords,
                        rho = rho,
                        min.rho = min.rho,
                        ...)

  #######################################################
  # Cleanup
  #######################################################
  knns[knns == 0] <- NA
  dist_matrix@x <- xs

  returnvalue <- list(
    knns = t(knns),
    wij = dist_matrix,
    call = sys.call(),
    coords = coords,
    sigmas = sqrt(sigmas / 2)
  )

  class(returnvalue) <- 'largeVis'
  return(returnvalue)
}
