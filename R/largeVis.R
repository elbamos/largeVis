#' Apply the LargeVis algorithm for visualizing large high-dimensional datasets.
#'
#' See \href{https://arxiv.org/abs/1602.00370}{Original paper}
#'
#' @param x A matrix
#' @param dim The number of dimensions in the output
#' @param K The number of nearest-neighbors to use in computing the graph
#' @param pca.first Whether to apply pca first (can speed-up distance calculations)
#' @param pca.dims How many pca dimensions to use
#' @param n.trees See \link{randomProjectionTreeSearch}
#' @param tree.threshold See \link{randomProjectionTreeSearch}
#' @param max.iter See \link{randomProjectionTreeSearch}
#' @param perplexity See paper
#' @param sgd.batches See \link{projectKNNs}
#' @param M See \link{projectKNNs}
#' @param weight.pos.samples See \link{projectKNNs}
#' @param alpha See \link{projectKNNs}
#' @param gamma See `projectKNNs`
#' @param verbose Verbosity
#' @param ... See paper
#'
#' @return A `largeVis` object with the following slots:
#'
#' @export
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
    if (verbose) cat("PCA...")
    shrunken.x <- princomp(x, scores = TRUE)$scores[,1:pca.dims]
    if (verbose) cat("done\n")
  }
  knns <- randomProjectionTreeSearch(shrunken.x,
                                     n.trees = n.trees,
                                     tree.threshold = tree.threshold,
                                     K = K,
                                     max.iter = max.iter,
                                     verbose = verbose)
  if (verbose) cat("Calculating edge weights...")
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
  rm(ts, wrongtri)
  # calculate distances
  if (verbose) cat("neighbor distances...")
  xs <- rep(0, length(is)) # pre-allocate

  if (verbose) cat("calculating...")
  distance(is, js, xs, shrunken.x)
  if (verbose) cat("done!\n")

  # assemble edges into graph, as symmetric sparse matrix

  dist_matrix <- Matrix::sparseMatrix(i = is,
                                      j = js,
                                      x = as.vector(xs),
                                      symmetric = TRUE)
  # select denominators sigma to match distributions to given perplexity
  if (verbose[1]) progress <- #utils::txtProgressBar(min = 0, max = sgd.batches, style = 3)
    progress::progress_bar$new(total = N, format = 'Estimate sigma [:bar] :percent eta: :eta', clear=FALSE)

  sigmas <- parallel::mclapply(1:N, FUN = function(idx) {
    x_i <- dist_matrix[idx,,drop=FALSE]
    x_i <- x_i[x_i > 0]
    optimize(f = function(sigma, x) {
      xs <- exp(-(x^2)/(sigma))
      softxs <- xs / sum(xs)
      p <- -sum(log2(softxs))/length(x)
      (log2(perplexity) - p)^2
    },
    x = x_i,
    interval = c(0,100))$minimum # note that sigmas are actually (2 * (sigmas^2))
  })
  sigmas <- unlist(sigmas)

  if (verbose) cat("p(j|i)...")
  pji <- Matrix::tril(dist_matrix) + Matrix::triu(dist_matrix)
  stopifnot(length(pji@i) == 2 * length(dist_matrix@i)) # confirm that matrix is no longer symmetric
  pji@x <- exp(-(pji@x^2) / sigmas[pji@i + 1])
  pji <- pji / Matrix::rowSums(pji)

  if (verbose) cat("w_{ji}...")
  # calculate w_{ij} by symmetrizing p(i|j) + p(j|i)/2N
  wij <- Matrix::forceSymmetric((pji + Matrix::t(pji)) / (2 * N))
  wij <- as(wij, "dgTMatrix")
  if (verbose) cat("Done!\n")
  # SGD PHASE
  coords <- projectKNNs(x = wij,
                        dim = dim,
                        sgd.batches = sgd.batches,
                        M = M,
                        weight.pos.samples = weight.pos.samples,
                        gamma = gamma,
                        verbose = verbose,
                        alpha = alpha,
                        ...)

  returnvalue <- list(
    knns = t(knns),
    pji = pji,
    call = sys.call(),
    coords = coords,
    sigmas = sqrt(sigmas / 2),
    wij = wij
  )

  class(returnvalue) <- 'largeVis'
  return(returnvalue)
}
