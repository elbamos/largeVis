#'  Project a distance matrix (presumably of k-nearest neighbors) into a lower-dimensional
#'  space using a proportionality function and stochastic gradient descent.
#'
#'  @param x A sparse matrix in triplet form
#'  @param dim The number of dimensions for the projection space
#'  @param sgd.batches The number of edges to process during SGD; defaults to 10000 * the number of rows in x
#'  @param neg.samples The number of negative edges to sample for each positive edge
#'  @param distance.function A function mapping the distance between two vertices in the lower-dimensional space to the probability that they are kNN's of each other.
#'  @param gamma Hyperparameter used in the objective function.  See the papers for details
#'  @param weight.pos.samples Whether to sample positive edges according to their edge weights (the default) or multiply the edge-loss by the edge-weight in the objective function.
#'  @return A dense [nrow(x),dim] matrix of the coordinates projecting x into the lower-dimensional space.

projectKNNs <- function(x, # a sparse distance matrix in triplet form
                               dim, # dimension of the projection space
                               sgd.batches = nrow(x) * 10000,
                               neg.samples = 5,
                               weight.pos.samples = TRUE,
                               distance.function =  # the probability that two nodes will be knns, given the distance (X) between them
                                  #'1 / (1 + (X)^2)',
                                 # '1 / (1 + (X)^2)',
                                 '1 / (1 + (2 * (X)^2))',
                                  #'1 / (1 + exp(- (X)^2))',
                               gamma = 7,
                               verbose = TRUE) {
    ########################################
    # GET FUNCTION FOR AUTODIFFERENTIATION
    ########################################
  nansaver <- NULL # " + 1e-10"
  peij <- function(kindex = 0) {
    distance <- list()
    for (d in 1:dim) {
      j <- ifelse(kindex == 0, "j", paste(sep="", "k", kindex))
      distance <- c(distance,
                    paste(sep = "",
                     "( yi", d, " - y", j, d, nansaver, ")^2")
      )
    }
    distance <- paste("sqrt(", distance, ")", collapse = "+")
    sub("X", distance, distance.function, fixed = TRUE)
  }

  objective <- paste("log(", peij(0), ") + gamma * (")
  for (k in 1:neg.samples) {
    objective <- paste(objective, "+ log(1 - ", peij(k), ")")
  }
  objective <- paste(objective, ")")
  if (weight.pos.samples) objective <- paste("~", objective)
  else objective <- paste("~ wij * (", objective, ")")

  # get autodifferentiation calls
  namevec <- c(paste(sep="", "yi", 1:dim),
               paste(sep="", "yj", 1:dim),
               paste('yk', sep = '', outer(1:neg.samples,1:dim, FUN='paste', sep = ''))
  ) # TODO: add wij to namevec
  print(objective)
  objfunc <- deriv(as.formula(objective), namevec = namevec, function.arg = c("gamma", namevec))
  print(namevec)

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

  klabels <- paste(sep = '', "yk", unlist(outer(1:neg.samples,1:dim, FUN=paste,sep="")))

  if (verbose) progress <- #utils::txtProgressBar(min = 0, max = sgd.batches, style = 3)
    progress::progress_bar$new(total = sgd.batches, format = 'SGD [:bar] :percent eta: :eta')

  #################################
  # SGD
  #################################

  for (eij in pos.edges) {
    if (verbose) if (require(progress)) {progress$tick()} else {utils::setTxtProgressBar(progress, utils::getTxtProgressBar(progress) + 1)}

    .args <- list(gamma = gamma)
    if (! weight.pos.samples) .args$wij <- wij@x[eij]

    i <- x@i[eij] + 1
    j <- x@j[eij] + 1

    # Half the time, transpose the edges
    if (rnorm(1) < 0) {
      t <- j
      j <- i
      i <- t
    }
    # Select negative samples
    availks <- which(x[i,] == 0)
    availks <- availks[availks != i] #     When a negative node is the identity node
    # EDGE CASES:
    #     When a node has no negative edges (i.e., is fully connected)
    ks <- sample(availks, neg.samples, replace = T, prob = neg.sample.weights[availks])

    .args[paste(sep = '', "y", unlist(outer(c("i", "j"), 1:2, FUN=paste, sep="")))] <- as.vector(coords[c(i,j),])
    .args[klabels] <- unlist((coords[ks,1:dim]))

    # Get gradients
    .val <-  do.call(objfunc, .args)
    grads <- attr(.val, 'gradient')

    if (any(is.nan(grads)) || is.nan(.val) || any(is.na(grads)) || is.na(.val) ||
        any(is.infinite(grads)) || is.infinite(.val)) {
      storage <<- list(
        val = .val,
        objective = objective,
        args = .args,
        coords = coords,
        i = i,
        j = j,
        ks = ks,
        objfunc = objfunc,
        namevec = namevec
      )
      stop()
    }

    # Update parameters, maximizing objective function
    coords[i,] <- coords[i,] + (grads[1:dim] * rho) # update i
    coords[j,] <- coords[j,] + (grads[(dim +1):(dim * 2)] * rho) # update j

    kgrads <- matrix(grads[(1 + (dim * 2)):length(grads)], ncol = dim)
    if (nrow(kgrads) > 0) {
      coords[ks,1:dim] <- coords[ks,1:dim] + (kgrads * rho)
    }
    .loss <- .loss + .val # Track training loss
    rho <- rho - (1 / sgd.batches)
  }
  return(coords)
}
