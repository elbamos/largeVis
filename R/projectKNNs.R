

#' Project a distance matrix into a lower-dimensional space.
#'
#' The input is a sparse triplet matrix showing the weights to give the edges, which are presumably estimated
#' k-nearest-neighbors.
#'
#' The algorithm attempts to estimate a \code{dim}-dimensional embedding using stochastic gradient descent and
#' negative sampling.
#'
#' The objective function is: \deqn{ O = \sum_{(i,j)\in E} w_{ij} (\log f(||p(e_{ij} = 1||) + \sum_{k=1}^{M} E_{jk~P_{n}(j)} \gamma \log(1 - f(||p(e_{ij_k} - 1||)))}
#' where \eqn{f()} is a probabilistic function relating the distance between two points in the low-dimensional projection space,
#' and the probability that they are nearest neighbors.  See the discussion of the alpha parameter below.

#' @param wij A sparse matrix of edge weights.
#' @param dim The number of dimensions for the projection space.
#' @param sgd.batches The number of edges to process during SGD; defaults to 20000 * the number of rows in x, as recommended
#' by the paper authors.
#' @param M The number of negative edges to sample for each positive edge.
#' @param alpha Hyperparameter used in the default distance function, \eqn{1 / (1 + \alpha \dot ||y_i - y_j||^2)}.  If \code{alpha} is 0, the alternative distance
#' function \eqn{1 / 1 + exp(||y_i - y_j||^2)} is used instead.  These functions relate the distance between points in the low-dimensional projection to the likelihood
#' that they two points are nearest neighbors. Note: the alternative probabilistic distance function is not yet implemented.
#' @param gamma Hyperparameter analogous to the strength of the force operating to push-away negative examples.
#' @param weight.pos.samples Whether to sample positive edges according to their edge weights (the default) or take the
#' weights into account when calculating gradient.  Note:  Applying weights to the gradients is not yet implemented.
#' @param rho Initial learning rate.
#' @param min.rho Final learning rate.
#' @param coords An initialized coordinate matrix.
#' @param verbose Verbosity
#'
#' @return A dense [nrow(x),dim] matrix of the coordinates projecting x into the lower-dimensional space.
#' @export
#' @importFrom stats rnorm
#'

projectKNNs <- function(wij, # sparse matrix
                        dim = 2, # dimension of the projection space
                        sgd.batches = nrow(N) * 20000,
                        M = 5,
                        weight.pos.samples = TRUE,
                        gamma = 7,
                        alpha = 2,
                        rho = 1,
                        coords = NULL,
                        min.rho = 0.1,
                        verbose = TRUE) {
  N <- length(wij@p) - 1
  nnzs <- diff(wij@p)
  js = rep(0:(N-1), diff(wij@p))
  is = wij@i


  ##############################################
  # Prepare vector of positive samples
  ##############################################
  pos.edges <- NULL

  if (weight.pos.samples) pos.edges <- sample(length(wij@x), sgd.batches, replace = T, prob = wij@x) - 1
  else pos.edges <- sample(length(wij@x), sgd.batches, replace = T) - 1

  ##############################################
  # Initialize coordinate matrix
  ##############################################
  if (is.null(coords)) coords <- matrix(rnorm(N * dim), nrow = dim)

  #################################################
  # SGD
  #################################################
  callback <- function(tick) {}
  progress <- #utils::txtProgressBar(min = 0, max = sgd.batches, style = 3)
    progress::progress_bar$new(total = sgd.batches, format = 'SGD [:bar] :percent/:elapsed eta: :eta', clear=FALSE)
  if (verbose[1]) callback <- progress$tick
  callback(0)
  coords <- sgd(coords,
              pos.edges,
              is = is,
              js = js,
              ps = wij@p,
              ws = wij@x,
              gamma = gamma, rho = rho, minRho = min.rho,
              useWeights = ! weight.pos.samples, M = M,
              alpha = alpha, callback = callback)

  return(coords)
}
