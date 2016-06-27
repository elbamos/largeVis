

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
#' The \code{weight_pos_samples} parameter controls how to handle edge-weights.  The paper authors recommend using a weighted
#' sampling approach to select edges, and treating edge-weight as binary in calculating the objective. This is the default.
#'
#' Note that the input matrix should be symmetric.  If any columns in the matrix are empty, the function will fail.
#'
#' @param wij A symmetric sparse matrix of edge weights, in C-compressed format, as created with the \code{Matrix} package.
#' @param dim The number of dimensions for the projection space.
#' @param sgd_batches The number of edges to process during SGD; defaults to 20000 * the number of rows in x, as recommended
#' by the paper authors.
#' @param M The number of negative edges to sample for each positive edge.
#' @param gamma The strength of the force pushing non-neighbor nodes apart.
#' @param alpha Hyperparameter used in the default distance function, \eqn{1 / (1 + \alpha \dot ||y_i - y_j||^2)}.  The function relates the distance
#' between points in the low-dimensional projection to the likelihood that the two points are nearest neighbors. Increasing \eqn{\alpha} tends
#' to push nodes and their neighbors closer together; decreasing \eqn{\alpha} produces a broader distribution. Setting \eqn{\alpha} to zero
#' enables the alternative distance function. \eqn{\alpha} below zero is meaningless.
#' @param weight_pos_samples Whether to sample positive edges according to their edge weights (the default) or take the
#' weights into account when calculating gradient.  See also the Details section.
#' @param rho Initial learning rate.
#' @param min_rho Final learning rate.
#' @param coords An initialized coordinate matrix.
#' @param verbose Verbosity
#'
#' @return A dense [N,D] matrix of the coordinates projecting the w_ij matrix into the lower-dimensional space.
#' @export
#' @examples
#' \dontrun{
#' data(wiki)
#' coords <- projectKNNs(wiki)
#' coords <- scale(coords)
#' plot(coords, xlim = c(-1.5,1.5), ylim = c(-1.5,1.5))
#' }
#' @importFrom stats rnorm
#'

projectKNNs <- function(wij, # symmetric sparse matrix
                        dim = 2, # dimension of the projection space
                        sgd_batches = (length(wij@p) - 1) * 20000,
                        M = 5,
                        weight_pos_samples = TRUE,
                        gamma = 7,
                        alpha = 1,
                        rho = 1,
                        coords = NULL,
                        min_rho = 0,
                        verbose = TRUE) {

  if (alpha < 0) stop("alpha < 0 is meaningless")
  if (alpha == 0) warning ("alpha == 0 may lead to numeric instability")
  N <-  (length(wij@p) - 1)
  js <- rep(0:(N - 1), diff(wij@p))
  if (any(is.na(js))) stop("NAs in the index vector.")
  is <- wij@i

  ##############################################
  # Initialize coordinate matrix
  ##############################################
  if (is.null(coords)) coords <- matrix(rnorm(N * dim), nrow = dim)

  #################################################
  # SGD
  #################################################
  if (verbose) cat("Estimating embeddings.\n")
  coords <- sgd(coords,
              is = is,
              js = js,
              ps = wij@p,
              ws = wij@x,
              gamma = gamma, rho = rho, minRho = min_rho,
              useWeights = ! weight_pos_samples, nBatches = sgd_batches,
              M = M, alpha = alpha, verbose = verbose)

  return(coords)
}
