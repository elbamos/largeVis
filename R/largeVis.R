
require(proxy)
require(foreach)

calcpji <- function(s, e, xsquareds) {
  pi <- exp(- xsquareds/(2 * s^2))
  pis <- pi / sum(pi)
  h <- -sum(pis * log2(pis))
  abs(e - h)
}

largeVis <- function(x,
                     dim = 2,
                     k = 150, # number of knn edges per vertex

                     pca.first = TRUE, # whether to apply dimensional reduction first
                     pca.dims = 50,

                     n.trees = 2, # in the random projection tree phase, how many trees to build
                     tree.threshold = k * 2, # the maximum number of nodes per leaf

                     max.iter = 1, # in the neighborhood exploration phase, the number of iterations

                     perplexity = k, # hyperparameter for calculating p(j|i)

                     sgd.max.iter = 10,
                     neg.samples = 5,
                     pos.samples = nrow(x) * 10,
                     distance.function = '1 / (1 + (X)^2)',
                     gamma = 7,

                     verbose = TRUE) {
  N = nrow(x)
  K = k
  # Handle pca.first
  shrunken.x <- x
  if (pca.first) {
    if (pca.dims >= ncol(x)) {
      warning("Called for dimensional reduction from ", ncol(x), " to ", pca.dims, " using pca.",
              "This is probably not what you intended.  Setting pca.dims = ", ncol(x))
      pca.dims <- ncol(x)
    }
    if (verbose) cat("PCA...")
    shrunken.x <- princomp(x, scores = TRUE)$scores[,1:pca.dims]
    if (verbose) cat("done\n")
  }

    # random projection trees
  tree_assignments <- list()
  if (verbose) cat("Creating random projection trees...")

  tree_assignments <- lapply(1:n.trees, FUN = function(T) {
    tree <- partition(indices = 1:N, .threshold = tree.threshold, .data = shrunken.x)
    knns <- matrix(0, nrow = tree.threshold, ncol = N)
    for (leaf in tree) {
      knns[1:length(leaf),leaf] <- rep(leaf, length(leaf))
    }
    knns
  })
  new_knns <- do.call(rbind, tree_assignments)
  # convert list of matrices to single [j,i] matrix

  if (verbose) cat("...done\n")

  if (verbose) cat("Exploring neighbors...")
  for (T in 1:max.iter) {
    if (verbose) cat(paste(T, " "))
    knns <- new_knns
    new_knns <- foreach(i = 1:N, .combine = cbind, .multicombine = TRUE) %do% {
      candidates <- c(knns[,i], knns[,knns[,i]])
      candidates <- candidates[candidates != i]
      candidates <- candidates[! candidates == 0]
      candidates <- unique(candidates)
      if (length(candidates) < K) {
        stop(paste("Attempted to visit neighbors but the number of knn candidates was smaller than K.",
                   "This can happen when there is a fully-connected distinct cluster in your dataset.",
                   "Try increasing the tree threshold or pca dimensions, or reducing K."))
      }
      distances <- proxy::dist(x = matrix(shrunken.x[i,], ncol = ncol(shrunken.x)),
                        y = shrunken.x[candidates,])
      candidates[order(distances)[1:K]]
    }
  }
  if (verbose) cat("done\n")

  # copy new_knns into adjacency matrix, where cells are the distance^2
  is <- rep(1:N, each = K)
  js <- as.vector(new_knns)
  # TODO: currently we get some distances twice

  xs <- proxy::dist(x = shrunken.x[is,],
                    y = shrunken.x[js,],
                    pairwise = TRUE)
  xs <- xs^2
  adj_matrix <- forceSymmetric(Matrix::spMatrix(
                                      nrow = N,
                                      ncol = N,
                                      i = c(is, js),
                                      j = c(js, is),
                                      x = rep(as.vector(xs), 2))) # TODO: we may have to check the diagonals

  # select denominators sigma to match distributions to given perplexity
  e <- log2(perplexity)
  i_sigmas <- lapply(1:N, FUN = function(idx) {
    xsquareds <- adj_matrix[idx,]
    xsquareds <- xsquareds[xsquareds > 0]
    sigma <- optimize(f = calcpji, e = e, xsquareds = xsquareds, interval = c(0,1))
    sigma$minimum
  })
  i_sigmas <- unlist(i_sigmas)
  i_sigmas <- 2 * (i_sigmas^2)
  i_sigmas <- i_sigmas[adj_matrix@i + 1]
  adj_matrix@x <- exp( - adj_matrix@x / i_sigmas)

  pji <- adj_matrix + t(adj_matrix) # so we get both directions
  pji <- pji / rowSums(pji) # softmax
  if (verbose) cat("Calculated p(j|i)...")

  # calculate w_{ij} by symmetrizing p(i|j) + p(j|i)/2N
  wij <-forceSymmetric((pji + t(pji)) / (2 * N))
  wij <- as(wij, "dgTMatrix")
  diag(wij) <- 0

  if (verbose) cat("w_{ij}\n")

  # SGD PHASE

  # create formulae for autodifferentiation
  insertion <- NULL
  for (i in 1:dim) {
    xx <- paste(sep="", "y", c("i", "j"))
    insertion <- c(insertion, paste("(", paste(xx, i, collapse = " - ", sep = ""), ")^2"))
  }
  distance <- paste("sqrt(", paste(insertion, collapse = "+"), ")")
  distance <- sub("X", distance, distance.function, fixed = TRUE)
#  posform <- paste("~ wij * log(", distance, ")")
  posform <- paste("~ log(", distance, ")") # treat samples as binary edges
  negform <- paste("~ gamma * log(1 - (", distance, "))")
  # get autodifferentiation calls
  namevec <- c(paste(sep="", "yi", 1:dim), paste(sep="", "yj", 1:dim))
  posfunc <- deriv(as.formula(posform), namevec = namevec)
  negfunc <- deriv(as.formula(negform), namevec = namevec)

  n_edges <- length(wij@x)
  degrees <- colSums(wij > 0)^0.75 # weights for each vertex when negative sampling j's from noisy distribution

  # initialize coordinate matrix
  .coords <- matrix(rnorm(nrow(x) * dim), ncol = dim)

  returnvalue <- list(
    knns = t(new_knns),
    pji = pji,
    call = sys.call(),
    losses <- list(),
    positivefunction <- posform,
    negativefunction <- negform
  )

  if (verbose) cat("Calculating coordinates by stochastic gradient descent...\n")
  for (T in 1:sgd.max.iter) {
    if (verbose) progress <- utils::txtProgressBar(min = 1, max = N + pos.samples, title = paste("Iteration", T))
    print(sum(.coords))
    .loss <- 0
    .dfunc.envir <- new.env()
    .dfunc.envir$gamma <- gamma

    learnRate <- 1 - (T / pos.samples)

    # sample pos.sample edges proportional to weight
    poses <- sample(length(wij@x), pos.samples, replace = T, prob = wij@x)

    for (.idx in poses) {
      if (verbose) utils::setTxtProgressBar(progress, utils::getTxtProgressBar(progress) + 1)

      # Get gradients for positive samples
     # .dfunc.envir$wij <- wij@x[.idx]
      i <- wij@i[.idx] + 1
      j <- wij@j[.idx] + 1

      for (d in 1:dim) {
        .dfunc.envir[[paste("yi", d, sep = "")]] <- .coords[i,d]
        .dfunc.envir[[paste("yj", d, sep = "")]] <- .coords[j,d]
      }
      .val <-  eval(posfunc, envir = .dfunc.envir)
      grads <- attr(.val, 'gradient')
      grads[is.nan(grads)] <- 0
      # .val <<- .val
      # stop()

      # Update coordinates
      .coords[i,] <- .coords[i,] - (grads[1:dim] * learnRate)
      .coords[j,] <- .coords[j,] - (grads[(dim +1):(dim * 2)] * learnRate)

      .loss <- .loss + .val # Track training loss
    }
    .loss <- .loss / length(poses)

    .negloss <- 0
    for (i in  (wij@i[poses] + 1)) {
      if (verbose) utils::setTxtProgressBar(progress, utils::getTxtProgressBar(progress) + 1)

      availjs <- which(wij[i,] == 0)
      js <- availjs
      if (length(availjs) > neg.samples) js <- sample(availjs, neg.samples, replace = F, prob = degrees[availjs])

      for (d in 1:dim) {
        .dfunc.envir[[paste("yi", d, sep = "")]] <- .coords[i,d]
        .dfunc.envir[[paste("yj", d, sep = "")]] <- .coords[js,d]
      }
      .val <- eval(negfunc, envir = .dfunc.envir)
      grads <- attr(.val, 'gradient')
      grads[is.nan(grads)] <- 0

#      grads[is.nan(grads)] <- 0
      # Update coordinate matrix
      .coords[i,] <- .coords[i,] - (colSums(grads)[1:dim] * learnRate)
      .coords[js,] <- .coords[js,] - (grads[,(dim +1):(dim * 2)] * learnRate)

      # Track training loss
      .negloss <- .negloss + sum(.val)/length(.val)
    }

    .negloss <- .negloss / N
    if (verbose) cat(paste("Epoch", T, "positive loss", .loss, "negative loss", .negloss, "\n"))
    returnvalue$losses[[T]] <- .loss
  }
  returnvalue$coords <- .coords
  class(returnvalue) <- 'largeVis'
  return(returnvalue)
}

