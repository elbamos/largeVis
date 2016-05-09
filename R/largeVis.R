#' Apply the LargeVis algorithm for visualizing large high-dimensional datasets.
#'
#' See \href{https://arxiv.org/abs/1602.00370}{Original paper}
#'
#' @param x A matrix
#' @param dim The number of dimensions in the output
#' @param K The number of nearest-neighbors to use in computing the graph
#' @param pca.first Whether to apply pca first (can speed-up distance calculations)
#' @param pca.dims How many pca dimensions to use
#' @param n.trees See \code\link{randomProjectionTreeSearch}
#' @param tree.threshold See \code\link{randomProjectionTreeSearch}
#' @param max.iter See \code\link{randomProjectionTreeSearch}
#' @param distance.method See \code\link{randomProjectionTreeSearch}
#' @param perplexity See paper
#' @param sgd.batches See \code\link{projectKNNs}
#' @param M See \code\link{projectKNNs}
#' @param weight.pos.samples See \code\link{projectKNNs}
#' @param distance.function See \code\link{projectKNNs}
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
largeVis <- function(x,
                     dim = 2,
                     K = 40, # number of knn edges per vertex

                     pca.first = TRUE, # whether to apply dimensional reduction first
                     pca.dims = 50,

                     n.trees = 2, # in the random projection tree phase, how many trees to build
                     tree.threshold = K * 2, #the maximum number of nodes per leaf
                     max.iter = 2, # in the neighborhood exploration phase, the number of iterations

                     distance.method = 'euclidean', # a distance method under by proxy::dist() for calculating the distance between neighbors

                     perplexity = K, # hyperparameter for calculating p(j|i)

                     sgd.batches = nrow(x) * 10000,
                     M = 5,
                     weight.pos.samples = TRUE,
                     distance.function = #'1 / (1 + (X)^2)',
                     #  '1 / (1 + (X)^2)',
                       '1 / (1 + (2 * (X)^2))',
                    #   '1 / (1 + exp(- (X)^2))',
                     gamma = 7,

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
                                     verbose = verbose,
                                     distance.method = distance.method,
                                     ...)

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
  xs <- proxy::dist(x = shrunken.x[is,],
                    y = shrunken.x[js,],
                    method='euclidean',
                    pairwise = TRUE)
  # assemble edges into graph, as symmetric sparse matrix
  dist_matrix <- Matrix::sparseMatrix(i = is,
                                      j = js,
                                      x = as.vector(xs),
                                      symmetric = TRUE)
  calcpji <- function(s, e, x) {
    xs <- (x^2)/s
    softxs <- exp(-x) / sum(exp(-x))
    p <- -sum(softxs * log2(softxs))
    (e - p)^2
  }

  # select denominators sigma to match distributions to given perplexity
  sigmas <- lapply(1:N,
                     FUN = function(idx) optimize(f = calcpji,
                                                  e = perplexity,
                                                  x = dist_matrix[idx,,drop=FALSE],
                                                  interval = c(0,2)))
  sigmas <- unlist(sigmas) # note that sigmas are actually (2 * (sigmas^2))
  pji <- tril(dist_matrix) + triu(dist_matrix)
  stopifnot(length(pji@i) == 2 * length(dist_matrix@i)) # confirm that matrix is no longer symmetric
  pji@x <- exp(-(pji@x^2) / sigmas[pji@i + 1])
  pji <- pji / rowSums(pji)

  if (verbose) cat("Calculated p(j|i)...")

  # calculate w_{ij} by symmetrizing p(i|j) + p(j|i)/2N
  wij <-forceSymmetric((pji + t(pji)) / (2 * N))
  wij <- as(wij, "dgTMatrix")

  if (verbose) cat("w_{ij}\n")

  # SGD PHASE
  coords <- projectKNNs(x = wij,
                        dim = dim,
                        sgd.batches = sgd.batches,
                        M = M,
                        weight.pos.samples = weight.pos.samples,
                        gamma = gamma,
                        verbose = verbose,
                        distance.function = distance.function)

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

