

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
#' @param distance.function A function mapping the distance between two vertices in the lower-dimensional space to the probability that they are kNN's of each other.
#'  \eqn{f(||y_1 - y_2||) = P(e_{ij} = 1)}
#' @param gamma Hyperparameter used in the objective function.  See the papers for details
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
                        distance.function =  # the probability that two nodes will be knns, given the distance (X) between them
                          # "1 / (1 + (X)^2)",
                          # "1 / (1 + (X)^2)",
                          "1 / (1 + (2 * (X)^2))",
                        # "1 / (1 + exp(- (X)^2))",
                        gamma = 7,
                        rho = 1,
                        min.rho = 0.1,
                        verbose = TRUE) {
  ########################################
  # GET FUNCTION FOR AUTODIFFERENTIATION
  ########################################
  nansaver <- NULL # " + 1e-10"
  peij <- function(kindex = 0) {
    distance <- list()
    for (d in 1:dim) {
      distance <- c(distance,
                    paste(sep = "",
                          "( yi", d, " - yj", d, nansaver, ")^2")
      )
    }
    distance <- paste("sqrt(", distance, ")", collapse = "+")
    sub("X", distance, distance.function, fixed = TRUE)
  }

  posobjective <- paste("~ log(", peij(0), ")")
  negobjective <- paste("~ gamma * log(1 - ", peij(0), ")")
  # TODO: ADD WIJ WEIGHTING

  # get autodifferentiation calls
  namevec <- c(paste(sep="", "yi", 1:dim),
              paste(sep="", "yj", 1:dim))
  # TODO: add wij to namevec
  print(posobjective)
  print(negobjective)
  print(namevec)

  posfunc <- deriv(as.formula(posobjective), namevec = namevec, function.arg = c(namevec))
  negfunc <- deriv(as.formula(negobjective), namevec = namevec, function.arg = c("gamma", namevec))

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
    if (verbose) if (require(progress)) {progress$tick()} else {utils::setTxtProgressBar(progress, utils::getTxtProgressBar(progress) + 1)}

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
    .args[paste(sep = '', "y", unlist(outer(c("i", "j"), 1:2, FUN=paste, sep="")))] <- as.vector(coords[c(i,j),])
    # Get gradients
    .val <-  do.call(posfunc, .args)
    grads <- attr(.val, 'gradient')
    .loss <- .loss + .val # Track training loss

    # Negative edges
    .args <- list(gamma = gamma)
    .args[paste(sep='', "yi", 1:dim)] <- as.vector(coords[i,])
    availjs <- which(x[i,] == 0)
    availjs <- availjs[availjs != i]
    js <- sample(availjs, min(M, length(availjs)), replace = F, prob = neg.sample.weights[availjs])
    if (length(js) > 0) {

      for (j in js) {
        .args[paste(sep='', "yj", 1:dim)] <- as.vector(coords[j,])
        .val <-  do.call(negfunc, .args)
        grads <- attr(.val, 'gradient')
        # Update parameters, maximizing objective function
        coords[i,] <- coords[i,] + (grads[1:dim] * rho / M) # update i
        coords[j,] <- coords[j,] + (grads[(dim +1):(dim * 2)] * rho / M) # update j
      }
      .loss <- .loss + .val # Track training loss
    }
    rho <- rho - ((rho - min.rho) / sgd.batches)

    counter <- counter + 1

    if (is.factor(verbose) && counter %% 50000 == 1) {
      plotcounter <- plotcounter + 1
      if (plotcounter %% 9 == 1) {
        x11()
        par(mfrow=c(3, 3))
      }
      plot(c(min(coords[,1]), max(coords[,1])),
           c(min(coords[,2]), max(coords[,2])),
           main = paste("Batch", counter))
      for (i in 1:nlevels(verbose)) {
        points(coords[as.integer(verbose) == i,],
               col = rainbow(nlevels(verbose))[i])
      }
    }
  }
  return(coords)
}
