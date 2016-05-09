

#' Project a distance matrix into a lower-dimensional space
#'
#' The input is a sparse triplet matrix showing the distances between vertices which are presumably estimated k-nearest-neighbors.
#' The algorithm attempts to estimate a \code{dim}-dimensional embedding using stochastic gradient descent.
#'
#' The objective function is: \deqn{ O = \sum_{(i,j)\in E} w_{ij} (\log p(e_{ij} = 1) + \sum_{k=1}^{M} E_{jk~P_{n}(j)} \gamma \log(1 - p(e_{ij_k} - 1)))  }

#' @param x A sparse matrix in triplet form
#' @param dim The number of dimensions for the projection space
#' @param sgd.batches The number of edges to process during SGD; defaults to 10000 * the number of rows in x
#' @param M The number of negative edges to sample for each positive edge
#' @param distance.function A function mapping the distance between two vertices in the lower-dimensional space to the probability that they are kNN's of each other.  (\eqn{f(||y_1 - y_2||) = P(e_{ij} = 1)}) The choices are a, \eqn{1 / (1 + \alpha * x^2)},
#' and b, \eqn{1 / (1 + exp(-x^2))}.
#' @param alpha Hyperparameter used in distance function a
#' @param gamma Hyperparameter controlling the weight given to each negative sample.
#' @param weight.pos.samples Whether to sample positive edges according to their edge weights (the default) or multiply the edge-loss by the edge-weight in the objective function.
#' @param rho Initial learning rate.
#' @param min.rho Final learning rate. The learning rate declines non-linearly.  \eqn{\rho_t = \rho_{t-1} - ((\rho_{t-1} - \rho_{min}) / sgd.batches)}
#' @param verbose Verbosity
#'
#' @return A dense [nrow(x),dim] matrix of the coordinates projecting x into the lower-dimensional space.
#' @export
#'

projectKNNs <- function(x, # a sparse distance matrix in triplet form
                        dim, # dimension of the projection space
                        sgd.batches = nrow(x) * 10000,
                        M = 5,
                        weight.pos.samples = TRUE,
                        distance.function =  "a",
                        gamma = 7,
                        alpha = 2,
                        rho = 1,
                        min.rho = 0.1,
                        verbose = TRUE) {

  ##################################
  # SETUP SGD
  ##################################
  neg.sample.weights <- colSums(x > 0)^0.75
  # Select positive samples
  pos.edges <- NULL
  if (weight.pos.samples) pos.edges <- sample(length(x@x), sgd.batches, replace = T, prob = x@x)
  else pos.edges <- sample(length(x@x), sgd.batches, replace = T)

  # initialize coordinate matrix
  coords <- matrix(rnorm(nrow(x) * dim), ncol = dim)
  .loss <- 0

  rho <- 1 # learning rate

  if (verbose) progress <- #utils::txtProgressBar(min = 0, max = sgd.batches, style = 3)
    progress::progress_bar$new(total = sgd.batches, format = 'SGD [:bar] :percent eta: :eta', clear=FALSE)

  #################################
  # SGD
  #################################
  counter <- 0
  plotcounter <- 0
  for (eij in pos.edges) {
    if (verbose) progress$tick()

    .args <- list()
    if (! weight.pos.samples) .args$wij <- wij@x[eij]

    i <- x@i[eij] + 1
    j <- x@j[eij] + 1

    # Half the time, transpose the edges
    if (rnorm(1) < 0) {
      t <- j
      j <- i
      i <- t
    }

    # Positive edge
    if (distance.function == "a") grads <- afuncpos(coords[i,], coords[j,], ifelse(weight.pos.samples, 1, wij@x[eij]), alpha)
    else grads <- bfuncpos(coords[i,], coords[j,],  ifelse(weight.pos.samples, 1, wij@x[eij]))
    coords[i,] <- coords[i,] + (grads * rho)
    coords[j,] <- coords[j,] - (grads * rho)

    availjs <- which(x[i,] == 0)
    availjs <- availjs[availjs != i]
    js <- sample(availjs, min(M, length(availjs)), replace = F, prob = neg.sample.weights[availjs])
    if (length(js) > 0) {
      for (j in js) {
        if (distance.function == "a") grads <- afuncneg(coords[i,], coords[j,], gamma, alpha)
        else grads <- bfuncneg(coords[i,], coords[j,],  gamma)
        coords[i,] <- coords[i,] + (grads * rho)
        coords[j,] <- coords[j,] - (grads * rho)
      }
    }
    rho <- rho - ((rho - min.rho) / (sgd.batches + 1))

    counter <- counter + 1

    if (is.factor(verbose) && counter %% 50000 == 1) {
      print(rho)
      plotcounter <- plotcounter + 1
      if (plotcounter %% 9 == 1) {
        x11()
        par(mfrow=c(3, 3))
      }
      plot(coords[as.integer(verbose) == 1,],
           col = rainbow(nlevels(verbose))[1],
           main = paste("Batch", counter))
      for (i in 2:nlevels(verbose)) {
        points(coords[as.integer(verbose) == i,],
               col = rainbow(nlevels(verbose))[i])
      }
    }
  }
  return(coords)
}

afuncpos <- function(i, j, wij, alpha) - (2 * wij * alpha * (i - j)) / (alpha * (sum((i-j)^2)) + 1)
afuncneg <- function(i,j,gamma, alpha) (2 * alpha * gamma * (i - j) ) / ((1 - (1 / (1 + (alpha * (sum((i-j)^2)))))) * (1 + (alpha * (sum((i-j)^2))))^2)
bfuncpos <- function(i,j,wij)  - (4 * wij * (i - j) * (sum((i^2)) + sum(j^2) - (2 * i * j)) ) / (1 + exp(sum((i - j)^2)^2))
bfuncneg <- function(i,j,gamma)  (4 * gamma * (i - j) * (sum((i^2)) + sum(j^2) - (2 * i * j)) ) * exp(sum((i - j)^2)^2) / (1 + exp(sum((i - j)^2)^2))

