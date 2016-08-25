

#' Project a distance matrix into a lower-dimensional space.
#'
#' Takes as input a sparse matrix of the edge weights connecting each node to its nearest neighbors, and outputs
#' a matrix of coordinates embedding the inputs in a lower-dimensional space.
#'
#' The algorithm attempts to estimate a \code{dim}-dimensional embedding using stochastic gradient descent and
#' negative sampling.
#'
#' The objective function is: \deqn{ O = \sum_{(i,j)\in E} w_{ij} (\log f(||p(e_{ij} = 1||) + \sum_{k=1}^{M} E_{jk~P_{n}(j)} \gamma \log(1 - f(||p(e_{ij_k} - 1||)))}
#' where \eqn{f()} is a probabilistic function relating the distance between two points in the low-dimensional projection space,
#' and the probability that they are nearest neighbors.
#'
#' The default probabilistic function is \eqn{1 / (1 + \alpha \dot ||x||^2)}. If \eqn{\alpha} is set to zero,
#' an alternative probabilistic function, \eqn{1 / (1 + \exp(x^2))} will be used instead.
#'
#' Note that the input matrix should be symmetric.  If any columns in the matrix are empty, the function will fail.
#'
#' @param wij A symmetric sparse matrix of edge weights, in C-compressed format, as created with the \code{Matrix} package.
#' @param dim The number of dimensions for the projection space.
#' @param sgd_batches The number of edges to process during SGD.
#' @param M The number of negative edges to sample for each positive edge.
#' @param gamma The strength of the force pushing non-neighbor nodes apart.
#' @param alpha Hyperparameter used in the default distance function, \eqn{1 / (1 + \alpha \dot ||y_i - y_j||^2)}.  The function relates the distance
#' between points in the low-dimensional projection to the likelihood that the two points are nearest neighbors. Increasing \eqn{\alpha} tends
#' to push nodes and their neighbors closer together; decreasing \eqn{\alpha} produces a broader distribution. Setting \eqn{\alpha} to zero
#' enables the alternative distance function. \eqn{\alpha} below zero is meaningless.
#' @param rho Initial learning rate.
#' @param coords An initialized coordinate matrix.
#' @param seed Random seed to be passed to the C++ functions; sampled from hardware entropy pool if \code{NULL} (the default).
#' Note that if the seed is not \code{NULL} (the default), the maximum number of threads will be set to 1 in phases of the algorithm
#' that would otherwise be non-deterministic.
#' @param threads The maximum number of threads to spawn. Determined automatically if \code{NULL} (the default).
#' @param verbose Verbosity
#'
#' @note If specified, \code{seed} is passed to the C++ and used to initialize the random number generator. This will not, however, be
#' sufficient to ensure reproducible results, because the initial coordinate matrix is generated using the \code{R} random number generator.
#' To ensure reproducibility, call \code{\link[base]{set.seed}} before calling this function, or pass it a pre-allocated coordinate matrix.
#'
#' @return A dense [N,D] matrix of the coordinates projecting the w_ij matrix into the lower-dimensional space.
#' @export
#' @importFrom stats runif
#'
projectKNNs <- function(wij, # symmetric sparse matrix
                        dim = 2, # dimension of the projection space
                        sgd_batches = NULL,
                        M = 5,
                        gamma = 7,
                        alpha = 1,
                        rho = 1,
                        coords = NULL,
												seed = NULL,
												threads = NULL,
                        verbose = getOption("verbose", TRUE)) {

  if (alpha < 0) stop("alpha < 0 is meaningless")
  N <-  (length(wij@p) - 1)
  js <- rep(0:(N - 1), diff(wij@p))
  if (any(is.na(js))) stop("NAs in the index vector.")
  is <- wij@i

  ##############################################
  # Initialize coordinate matrix
  ##############################################
  if (is.null(coords)) #coords <- matrix(rnorm(N * dim), nrow = dim)
                        coords <- matrix((runif(N * dim) - 0.5) / dim * 0.0001, nrow = dim)
  if (is.null(sgd_batches)) {
    if (N < 10000) {
      sgd_batches <- 20000 * length(wij@x)
    } else if (N < 1000000) {
      sgd_batches <- (N - 10000) * 9000 / (1000000 - 10000) + 1000
      sgd_batches <- sgd_batches * 1000000
    } else {
      sgd_batches <- 10000 * N
    }
  }

  #################################################
  # SGD
  #################################################
  if (verbose) cat("Estimating embeddings.\n")
  coords <- sgd(coords,
                targets_i = is,
                sources_j = js,
                ps = wij@p,
                weights = wij@x,
                alpha = alpha, gamma = gamma, M = M,
                rho = rho,
                n_samples = sgd_batches,
  							seed = seed,
  							threads = threads,
                verbose = verbose)

  return(coords)
}
