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
#' @param sgd_batches The number of edges to process during SGD. Defaults to a value set based on the size of the dataset. If the parameter given is
#' between \code{0} and \code{1}, the default value will be multiplied by the parameter.
#' @param M The number of negative edges to sample for each positive edge.
#' @param gamma The strength of the force pushing non-neighbor nodes apart.
#' @param alpha Hyperparameter used in the default distance function, \eqn{1 / (1 + \alpha \dot ||y_i - y_j||^2)}.  The function relates the distance
#' between points in the low-dimensional projection to the likelihood that the two points are nearest neighbors. Increasing \eqn{\alpha} tends
#' to push nodes and their neighbors closer together; decreasing \eqn{\alpha} produces a broader distribution. Setting \eqn{\alpha} to zero
#' enables the alternative distance function. \eqn{\alpha} below zero is meaningless.
#' @param rho Initial learning rate.
#' @param coords An initialized coordinate matrix.
#' @param useDegree Whether to use vertex degree to determine weights in negative sampling (if \code{TRUE}), or the sum of the vertex's edges (the default). See Notes.
#' @param momentum If not \code{NULL} (the default), SGD with momentum is used, with this multiplier, which must be between 0 and 1. Note that
#' momentum can drastically speed-up training time, at the cost of additional memory consumed.
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
#' @note The original paper called for weights in negative sampling to be calculated according to the degree of each vertex, the number of edges
#' connecting to the vertex. The reference implementation, however, uses the sum of the weights of the edges to each vertex. In experiments, the
#' difference was imperceptible with small (MNIST-size) datasets, but the results seems aesthetically preferrable using degree. The default
#' is to use the edge weights, consistent with the reference implementation.
#'
#' @return A dense [N,D] matrix of the coordinates projecting the w_ij matrix into the lower-dimensional space.
#' @export
#' @importFrom stats runif
#' @examples
#' data(CO2)
#' CO2$Plant <- as.integer(CO2$Plant)
#' CO2$Type <- as.integer(CO2$Type)
#' CO2$Treatment <- as.integer(CO2$Treatment)
#' co <- scale(as.matrix(CO2))
#' # Very small datasets often produce a warning regarding the alias table.  This is safely ignored.
#' suppressWarnings(vis <- largeVis(t(co), K = 20, sgd_batches = 1, threads = 2))
#' suppressWarnings(coords <- projectKNNs(vis$wij, threads = 2))
#' plot(t(coords))
projectKNNs <- function(wij, # symmetric sparse matrix
                        dim = 2, # dimension of the projection space
                        sgd_batches = NULL,
                        M = 5,
                        gamma = 7,
                        alpha = 1,
                        rho = 1,
                        coords = NULL,
												useDegree = FALSE,
												momentum = NULL,
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
  if (is.null(coords)) coords <- matrix((runif(N * dim) - 0.5) / dim * 0.0001, nrow = dim)

  if (! is.null(sgd_batches) && sgd_batches < 0) stop("sgd batches must be > 0")
  if (! is.null(sgd_batches) && sgd_batches < 1) {
  	multiplier <- sgd_batches
  	sgd_batches <- NULL
  } else multiplier <- 1

  if (is.null(sgd_batches)) sgd_batches <- sgdBatches(N, length(wij@x / 2))
  sgd_batches <- sgd_batches * multiplier

  if (!is.null(threads)) threads <- as.integer(threads)
  if (!is.null(momentum)) momentum <- as.double(momentum)

  #################################################
  # SGD
  #################################################
  if (verbose) cat("Estimating embeddings.\n")
  coords <- sgd(coords,
                targets_i = is,
                sources_j = js,
                ps = wij@p,
                weights = wij@x,
                alpha = as.double(alpha), gamma = as.double(gamma), M = as.integer(M),
                rho = as.double(rho),
                n_samples = sgd_batches,
  							momentum = momentum,
  							seed = seed,
  							threads = threads,
                verbose = as.logical(verbose))

  return(coords)
}

#' sgdBatches
#'
#' Calculate the default number of batches for a given number of vertices and edges.
#'
#' The formula used is the one used by the \code{LargeVis} reference implementation.  This is substantially less than the recommendation \eqn{E * 10000} in the original paper.
#'
#' @param N Number of vertices.
#' @param E Number of edges.
#'
#' @return The recommended number of sgd batches.
#' @export
#'
#' @examples
#' # Observe that increasing K has no effect on processing time
#' N <- 70000 # MNIST
#' K <- 10:250
#' plot(K, sgdBatches(rep(N, length(K)), N * K / 2))
#'
#' # Observe that processing time scales linarly with N
#' N <- c(seq(from = 1, to = 10000, by = 100), seq(from = 10000, to = 10000000, by = 1000))
#' plot(N, sgdBatches(N))
sgdBatches <- function(N, E = 150 * N / 2) {
	ifelse(N < 10000, 2000 * E, ifelse(N < 1000000, 1000000 * (9000 * (N - 10000) / (1000000 - 10000) + 1000), N * 10000))
}
