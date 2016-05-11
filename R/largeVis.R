#' Apply the LargeVis algorithm for visualizing large high-dimensional datasets.
#'
#' Implements the \code{largeVis} algorithm by Tang et al.
#'
#' \code{largeVis} estimates a low-dimensional embedding for high-dimensional data, where the distance between vertices
#' in the low-dimensional space is proportional to the distance between them in the high-dimensional space. The algorithm
#' works in 4 phases:
#'
#' \itemize {
#' \item  Estimate candidate nearest-neighbors for each vertex by building \code{n.trees} random projection trees.
#' \item  Estimate \code{K} nearest-neighbors for each vertex by visiting each vertex' 2d-degree neighbors (its neighbors' neighbors).
#' This is repeated \code{max.iter} times.  Note that the original paper suggested a \code{max.iter} of 1, however a larger number
#' may be appropriate for some datasets if the algorithm has trouble finding K neighbors for every vertex.
#' \item Estimate \eqn{p_{j|i}}, the conditional probability that each edge found in the previous step is actually to a
#' nearest neighbor of each of its nodes.
#' \item Using stochastic gradient descent, estimate an embedding for each vertex in the low-dimensional space.
#' }
#'
#' @param x A matrix
#' @param dim The number of dimensions in the output
#' @param K The number of nearest-neighbors to use in computing the graph
#' @param pca.first Whether to apply pca first (can speed-up distance calculations)
#' @param pca.dims How many pca dimensions to use
#' @param n.trees See \code{\link{randomProjectionTreeSearch}}
#' @param tree.threshold See \code{\link{randomProjectionTreeSearch}}
#' @param max.iter See \code{\link{randomProjectionTreeSearch}}
#' @param perplexity See paper
#' @param sgd.batches See \code{\link{projectKNNs}}
#' @param M See \code{\link{projectKNNs}}
#' @param weight.pos.samples See \code{\link{projectKNNs}}
#' @param alpha See \code{\link{projectKNNs}}
#' @param gamma See `projectKNNs`
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
#'
#' @examples
#'
#' @useDynLib largeVis
#' @importFrom Rcpp sourceCpp
#'
largeVis <- function(x,
                     dim = 2,
                     K = 40, # number of knn edges per vertex

                     pca.first = TRUE, # whether to apply dimensional reduction first
                     pca.dims = 50,

                     n.trees = 2, # in the random projection tree phase, how many trees to build
                     tree.threshold = K * 2, #the maximum number of nodes per leaf
                     max.iter = 2, # in the neighborhood exploration phase, the number of iterations

                     perplexity = K, # hyperparameter for calculating p(j|i)

                     sgd.batches = nrow(x) * 10000,
                     M = 5,
                     weight.pos.samples = TRUE,
                     alpha = 1,
                     gamma = 7,
                     rho = 1,
                     min.rho = 0,

                     verbose = TRUE,
                    ...) {
  N = nrow(x)
  # Handle pca.first
  shrunken.x <- x
  if (pca.first) {
    if (pca.dims >= ncol(x)) stop("Called for dimensional reduction from ", ncol(x), " to ", pca.dims, " using pca.")
    if (verbose[1]) cat("PCA...")
    shrunken.x <- princomp(x, scores = TRUE)$scores[,1:pca.dims]
    if (verbose[1]) cat("done\n")
  }
  knns <- randomProjectionTreeSearch(shrunken.x,
                                     n.trees = n.trees,
                                     tree.threshold = tree.threshold,
                                     K = K,
                                     max.iter = max.iter,
                                     verbose = verbose)
  if (sum(colSums(knns) == 0) > 0) stop("Found no neighbors for some nodes.")
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

  # calculate distances
  if (verbose[1]) cat("neighbor distances...")
  xs <- rep(0, length(is)) # pre-allocate

  if (verbose[1]) cat("calculating...")
  distance(is, js, xs, shrunken.x)
  if (verbose) cat("done!\n")

  # assemble edges into graph, as symmetric sparse matrix

  dist_matrix <- Matrix::sparseMatrix(i = is,
                                      j = js,
                                      x = xs,
                                      symmetric = TRUE)
  # select denominators sigma to match distributions to given perplexity
  if (verbose[1]) cat("Estimating sigmas...")
  suppressMessages(sigmas <- parallel::mclapply(1:N, FUN = function(idx) {
    x_i <- dist_matrix[idx,,drop=FALSE]
    x_i <- x_i[x_i > 0]
    optimize(f = function(sigma, xi) {
      lxs <- exp(-(xi^2)/(sigma))
      softxs <- lxs / sum(lxs)
      p <- -sum(log2(softxs))/length(xi)
      (log2(perplexity) - p)^2
    },
    x = x_i,
    interval = c(0,100))$minimum # note that sigmas are actually (2 * (sigmas^2))
  }))
  sigmas <- unlist(sigmas)
  if (verbose[1]) cat("done!\n")
  wijVector = rep(0, length(xs * 2))
  if (verbose[1]) progress <- #utils::txtProgressBar(min = 0, max = sgd.batches, style = 3)
    progress::progress_bar$new(total = N, format = 'Calculate p_{j|i} and w_{ij} [:bar] :percent eta: :eta', clear=FALSE)$tick
  else progress <- function(tick) {}

  distMatrixTowij(is, js, xs, sigmas, wijVector, N, progress)

  # SGD PHASE
  coords <- projectKNNs(i = is, j = js, x = wijVector,
                        dim = dim,
                        sgd.batches = sgd.batches,
                        M = M,
                        weight.pos.samples = weight.pos.samples,
                        gamma = gamma,
                        verbose = verbose,
                        alpha = alpha,
                        ...)

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
