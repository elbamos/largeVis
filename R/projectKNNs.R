

#' Project a distance matrix into a lower-dimensional space
#'
#' The input is a sparse triplet matrix showing the distances between vertices which are presumably estimated k-nearest-neighbors.
#' The algorithm attempts to estimate a \code{dim}-dimensional embedding using stochastic gradient descent.
#'
#' The objective function is: \deqn{ O = \sum_{(i,j)\in E} w_{ij} (\log p(e_{ij} = 1) + \sum_{k=1}^{M} E_{jk~P_{n}(j)} \gamma \log(1 - p(e_{ij_k} - 1)))  }

#' @param i The i-vector component of a sparse triplet matrix
#' @param j the j-vector component of a sparse triplet matrix
#' @param x the x-vector component of a sparse triplet matrix. See \link{Matrix::sparseMatrix}
#' @param dim The number of dimensions for the projection space
#' @param sgd.batches The number of edges to process during SGD; defaults to 10000 * the number of rows in x
#' @param M The number of negative edges to sample for each positive edge
#' @param alpha Hyperparameter used in the default distance function, \eqn{1 / (1 + \alpha \dot ||y_i - y_j||^2)}.  If \code{alpha} is 0, the alternative distance
#' function \eqn{1 / 1 + exp(||y_i - y_j||^2)} is used instead.  These functions relate the distance between points in the low-dimensional projection to the likelihood
#' that they two points are nearest neighbors.
#' @param gamma Hyperparameter analogous to the strength of the force operating to push-away negative examples.
#' @param weight.pos.samples Whether to sample positive edges according to their edge weights (the default) or multiply the edge-loss by the edge-weight in the objective function.
#' @param rho Initial learning rate.
#' @param min.rho Final learning rate. The learning rate declines non-linearly.  \eqn{\rho_t = \rho_{t-1} - ((\rho_{t-1} - \rho_{min}) / sgd.batches)}
#' @param .coords An initialized coordinate matrix.
#' @param verbose Verbosity
#'
#' @return A dense [nrow(x),dim] matrix of the coordinates projecting x into the lower-dimensional space.
#' @export
#'

projectKNNs <- function(wij, # sparse matrix
                        dim, # dimension of the projection space
                        sgd.batches = nrow(N) * 10000,
                        M = 5,
                        weight.pos.samples = TRUE,
                        gamma = 7,
                        alpha = 2,
                        rho = 1,
                        .coords = NULL,
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
  # Prepare SGD
  ##############################################
  if (is.null(.coords)) .coords <- matrix(rnorm(N * dim), ncol = dim)

  plotcounter <- 0

  callback <- function(tick) {}
  progress <- #utils::txtProgressBar(min = 0, max = sgd.batches, style = 3)
    progress::progress_bar$new(total = sgd.batches, format = 'SGD [:bar] :percent/:elapsed eta: :eta', clear=FALSE)

  if (verbose[1] && is.factor(verbose)) callback <- function(tick) {
    progress$tick(tick)
    if (is.factor(verbose) && tick %% 50000 == 1) {
      plotcounter <- plotcounter + 1
      if (plotcounter %% 20 == 1) {
        x11()
        par(mfrow=c(5, 4))
      }
      plot(coords[as.integer(verbose) == 1,],
           col = rainbow(nlevels(verbose))[1],
           main = paste("Batch", counter))
      for (i in 2:nlevels(verbose)) {
        points(.coords[as.integer(verbose) == i,],
               col = rainbow(nlevels(verbose))[i])
      }
    }
  }
  else if (verbose[1]) callback <- progress$tick

  #################################################
  # SGD
  #################################################
  callback(0)
  sgd(.coords,
      pos.edges,
      is = is,
      js = js,
      ps = wij@p,
      ws = wij@x,
      gamma = gamma, rho = rho, minRho = min.rho,
      useWeights = ! weight.pos.samples, M = M,
      alpha = alpha, callback = callback)

  return(.coords)
}
