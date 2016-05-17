#' Apply the LargeVis algorithm for visualizing large high-dimensional datasets.
#'
#' Implements the \code{vis}
#'
#' Note that this implementation expects the data to be free of \code{NaN}'s, \code{NA}'s, \code{Inf}'s, and duplicate rows.
#' If any of these assumptions are violated, the algorithm will fail. It is also usually a good idea to scale the input data
#' to have unit norm and mean 0. If there are large values in the input matrix, some computations may oveflow.
#'
#' @param x A matrix. Ideally, the columns should be scaled and normalized to avoid the risk of errors caused by overflow.
#' @param dim The number of dimensions in the output
#' @param K The number of nearest-neighbors to use in computing the kNN graph
#' @param check.assumptions Whether to check the input matrix for duplicates, \code{NA}`s, etc.
#' @param pca.first Whether to apply pca first (can speed-up distance calculations)
#' @param pca.dims How many pca dimensions to use
#' @param n.trees See \code{\link{randomProjectionTreeSearch}}.  The default is set at 50, which is the number
#' used in the examples in the original paper.
#' @param tree.threshold See \code{\link{randomProjectionTreeSearch}}.  By default, this is the number of features
#' in the input set, which is the setting used in the examples in the original paper.  Note the time and memory requirements:
#' the first pass through the neighborhood exploration phases will involve up to \eqn{N * nTrees * threshold} comparisons.
#' @param max.iter See \code{\link{randomProjectionTreeSearch}}.
#' @param perplexity See paper
#' @param sgd.batches See \code{\link{projectKNNs}}.
#' @param M See \code{\link{projectKNNs}}.
#' @param weight.pos.samples See \code{\link{projectKNNs}}.
#' @param alpha See \code{\link{projectKNNs}}.
#' @param gamma See \code{\link{projectKNNs}}.
#' @param rho See \code{\link{projectKNNs}}.
#' @param min.rho \code{\link{projectKNNs}}.
#' @param coords A [N,K] matrix of coordinates to use as a starting point -- useful for refining an embedding in stages.
#' @param verbose Verbosity
#' @param ... See paper
#'
#' @return A `largeVis` object with the following slots:
#'  \itemize{
#'    \item{'knns'} {An [N,K] integer matrix, which is an adjacency list of each vertex' identified nearest neighbors.
#'    If the algorithm failed to find \code{K} neighbors, the matrix is padded with \code{NA}'s.}
#'    \item{'wij'} {A sparse [N,N] matrix where each cell represents \eqn{w_{ij}}.}
#'    \item{'call'}
#'    \item{'coords'} {A [N,D] matrix of the embedding of the dataset in the low-dimensional space.}
#'  }
#'
#'
#' @export
#' @importFrom stats optimize princomp
#' @references Jian Tang, Jingzhou Liu, Ming Zhang, Qiaozhu Mei. \href{https://arxiv.org/abs/1602.00370}{Visualizing Large-scale and High-dimensional Data.}
#'
#' @examples
#' # iris
#' data(iris)
#' dat <- as.matrix(iris[,1:4])
#' dat <- scale(dat)
#' dupes = which(duplicated(dat))
#' dat <- dat[-dupes,] # duplicated data potentially can cause the algorithm to fail
#' visObject <- vis(dat, pca.first = FALSE,
#'                      max.iter = 20, sgd.batches = 800000,
#'                      K = 10,  gamma = 2, rho = 1, M = 40, alpha = 20,verbose=FALSE)
#'
#'  # mnist
#'  load("./mnist.Rda")
#'  dat <- mnist$images
#'  dim(dat) <- c(42000, 28 * 28)
#'  dat <- (dat / 255) - 0.5
#'  coords <- vis(dat, pca.first = FALSE,
#'                   n.tree = 10, tree.threshold = 40,
#'                   K = 40, sgd = 20000 * 42000, alpha = 1, max.iter = 10)
#'
#'
vis <- function(x,
                     dim = 2,
                     K = 40,

                     check.assumptions = TRUE,
                     pca.first = TRUE,
                     pca.dims = 50,

                     n.trees = 50,
                     tree.threshold = max(10, if (pca.first) {pca.dims} else {ncol(x)}),
                     max.iter = 3,

                     perplexity = 50,

                     sgd.batches = nrow(x) * 20000,
                     M = 5,
                     weight.pos.samples = TRUE,
                     alpha = 2,
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
    if ((any(is.na(shrunken.x)) + any(is.infinite(shrunken.x)) + any(is.nan(shrunken.x))) > 0)
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
  # Clean knns
  #############################################
  if (verbose[1]) cat("Calculating edge weights...")
  # These vectors are analogous to the components of a sparse matrix,
  # but both triple and C-compressed forms are created.
  # The i and j vectors are 0-indexed while p is 1-indexed.
  is <- rep(0:(N-1), each = K)
  js <- as.vector(knns)
  is <- is[! js == -1]
  js <- js[! js == -1]
  dupes <- duplicated(data.frame(is, js))
  is <- is[! dupes]
  js <- js[! dupes]
  ord <- order(is)
  is <- is[ord]
  js <- js[ord]
  ps <- rep(NA,N + 1)
  diffs <- diff(is)
  ps[is[which(diffs > 0)] + 2] <- which(diffs > 0) + 1
  good <- cumsum(!is.na(ps))
  ps <- ps[good + 1]
  ps[1] <- 1
  ps[length(ps) + 1] <- length(is)

  #######################################################
  # Calculate edge weights for candidate neighbors
  #######################################################
  if (verbose[1]) ptick <- progress::progress_bar$new(total = N, format = 'Calculate neighbor distances [:bar] :percent/:elapsed eta: :eta', clear=FALSE)$tick
  else ptick <- function(tick) {}

  xs <- rep(0, length(is)) # pre-allocate
  distance(is, js, xs, shrunken.x, ptick)

  if (verbose) cat("\n")
  if ((any(is.na(xs)) + any(is.infinite(xs)) + any(is.nan(xs)) + any(xs == 0)) > 0)
    stop("An error leaked into the distance calculation - check for duplicates")

  ########################################################
  # Estimate sigmas
  ########################################################
  if (verbose[1]) ptick <- progress::progress_bar$new(total = N, format = 'Calculate sigmas [:bar] :percent/:elapsed eta: :eta', clear=FALSE)$tick
  else ptick <- function(tick) {}
  # TODO: MAKE SURE THAT THE C++ CODE HANDLES THE SITUATION WHERE A COLUMN IS EMPTY (HAS NO PRESENCE IN P)

  perplexity = log2(perplexity)
  sigmas <- parallel::mclapply(1:N, FUN = function(idx) {
    ptick(1)
    x_i <- xs[ps[idx]:(ps[idx + 1] - 1)]
    ret <- optimize(f = sigFunc,
             x = x_i,
             perplexity = perplexity,
             interval = c(0,10000))
  })
  sigmas <- sapply(sigmas, `[[`, 1)

  if (any(is.na(sigmas)) + any(is.infinite(sigmas)) + any(is.nan(sigmas)) + any((sigmas == 0)) > 0)
    stop("An error has propogated into the sigma vector.")


  #######################################################
  # Calculate w_{ij}
  #######################################################
  if (! requireNamespace(Matrix,quietly=T)) stop("The Matrix package must be available.")
  if (verbose[1]) progress <- progress::progress_bar$new(total = length(xs) * 2, format = 'Calculate p_{j|i} and w_{ij} [:bar] :percent/:elapsed eta: :eta', clear=FALSE)$tick
  else progress <- function(tick) {}

  wij <- distMatrixTowij(is, js, xs, sigmas, N, progress)

  if (any(is.na(wij@x)) || any(is.infinite(wij@x)) || any(is.nan(wij@x)) || any((wij@x == 0)) > 0)
    stop("An error has propogated into the w_{ij} vector.  This probably means the input data wasn't scaled.")

  # Symmetricize
  wij <- wij + Matrix::t(wij)

  #######################################################
  # Estimate embeddings
  #######################################################
  coords <- projectKNNs(wij = wij,
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
  knns[knns == -1] <- NA

  returnvalue <- list(
    knns = t(knns),
    wij = wij,
    call = sys.call(),
    coords = t(coords),
    sigmas = sqrt(sigmas / 2)
  )

  class(returnvalue) <- 'largeVis'
  return(returnvalue)
}

##########################################
# Some helper functions useful for debugging
##########################################

pji <- function(x_i, sigma)  exp(- (x_i^2) / sigma) / sum(exp(- (x_i^2) / sigma))
perp <- function(x_i, sigma) - sum(log2(pji(x_i, sigma))) / length(x_i)
pdiff <- function(x_i, sigma, perplexity) (perplexity - perp(x_i, sigma))^2
