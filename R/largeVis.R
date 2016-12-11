#' Apply the LargeVis algorithm for visualizing large high-dimensional datasets.
#'
#' @param x A matrix, where the features are rows and the examples are columns.
#' @param dim The number of dimensions in the output
#' @param K The number of nearest-neighbors to use in computing the kNN graph
#' @param n_trees See \code{\link{randomProjectionTreeSearch}}.  The default is set at 50, which is the number
#' used in the examples in the original paper.
#' @param tree_threshold See \code{\link{randomProjectionTreeSearch}}.  By default, this is the number of features
#' in the input set.
#' @param max_iter See \code{\link{randomProjectionTreeSearch}}.
#' @param distance_method One of "Euclidean" or "Cosine."  See \code{\link{randomProjectionTreeSearch}}.
#' @param perplexity See \code{\link{buildWijMatrix}}.
#' @param save_neighbors Whether to include in the output the adjacency matrix of nearest neighbors.
#' @param save_edges Whether to include in the output the distance matrix of nearest neighbors.
#' @param threads The maximum number of threads to spawn. Determined automatically if \code{NULL} (the default).  It is unlikely that
#' this parameter should ever need to be adjusted.  It is only available to make it possible to abide by the CRAN limitation that no package
#' use more than two cores.
#' @param verbose Verbosity
#' @param ... Additional arguments passed to \code{\link{projectKNNs}}.
#'
#' @return A `largeVis` object with the following slots:
#'  \describe{
#'    \item{'knns'}{If \code{save_neighbors=TRUE}, An [N,K] 0-indexed integer matrix, which is an adjacency list of each vertex' identified nearest neighbors.
#'    If the algorithm failed to find \code{K} neighbors, the matrix is padded with \code{NA}'s. Note that this matrix is not identical to the output
#'    from \code{\link{randomProjectionTreeSearch}}: missing neighbors are \code{NA}'s rather than \code{-1}'s, and the matrix is transposed.}
#'    \item{'edges'}{If \code{save_edges=TRUE}, a [N,N] sparse matrix of distances between nearest neighbors.}
#'    \item{'wij'}{A sparse [N,N] matrix where each cell represents \eqn{w_{ij}}.}
#'    \item{'call'}{The call.}
#'    \item{'coords'}{A [D,N] matrix of the embedding of the dataset in the low-dimensional space.}
#'  }
#'
#' @export
#' @references Jian Tang, Jingzhou Liu, Ming Zhang, Qiaozhu Mei. \href{https://arxiv.org/abs/1602.00370}{Visualizing Large-scale and High-dimensional Data.}
#'
#' @examples
#' \dontrun{
#' # iris
#' data(iris)
#' dat <- as.matrix(iris[,1:4])
#' dat <- scale(dat)
#' dupes = which(duplicated(dat))
#' dat <- dat[-dupes,] # duplicates can cause the algorithm to fail
#' dat <- t(dat)
#' visObject <- largeVis(dat, max_iter = 20, K = 10, sgd_batches = 10000)
#' plot(t(visObject$coords))
#'
#' # mnist
#' # Note: The MNIST dataset may be obtained using the deepnet package.
#' load("./mnist.Rda")
#' dat <- mnist$images
#' dim(dat) <- c(42000, 28 * 28)
#' dat <- (dat / 255) - 0.5
#' dat <- t(dat)
#' visObject <- largeVis(dat, n_trees = 50, tree_th = 200, K = 50)
#' plot(t(visObject$coords))
#' }
#'
largeVis <- function(x,
                     dim = 2,
                     K = 50,

                     n_trees = 50,
                     tree_threshold = max(10, ncol(x)),
                     max_iter = 1,
                     distance_method = "Euclidean",

                     perplexity = max(50, K / 3),

                     save_neighbors = TRUE,
										 save_edges = TRUE,

										 threads = NULL,

                     verbose = getOption("verbose", TRUE),
                    ...) {

	if (!(is.matrix(x) && is.numeric(x)) && !is.data.frame(x)) stop("LargeVis requires a matrix or data.frame")
	if (is.data.frame(x)) x <- t(as.matrix(x[, sapply(x, is.numeric)]))

  #############################################
  # Search for kNearestNeighbors
  #############################################
  knns <- randomProjectionTreeSearch(x,
                                     n_trees = n_trees,
                                     tree_threshold = tree_threshold,
                                     K = K,
                                     max_iter = max_iter,
                                     distance_method = distance_method,
  																	 threads,
                                     verbose = verbose)
  #############################################
  # Clean knns
  #############################################
  if (verbose[1]) cat("Calculating edge weights...\n")
  edges <- buildEdgeMatrix(data = x,
  												 neighbors = knns,
  												 distance_method = distance_method,
  												 verbose = verbose)
  if (!save_neighbors) rm(knns)
  gc()
  if (any(edges$x > 27)) {
  	warning(paste(
  		"The Distances between some neighbors are large enough to cause the calculation of p_{j|i} to overflow.",
  		"Scaling the distance vector."))
  	edges$x <- edges$x / max(edges$x)
  }
  wij <- buildWijMatrix(edges, threads, perplexity)
  if (!save_edges) rm(edges)

  #######################################################
  # Estimate embeddings
  #######################################################
  coords <- projectKNNs(wij = wij,
                        dim = dim,
                        verbose = verbose,
  											threads = threads,
                        ...)

  #######################################################
  # Cleanup
  #######################################################

  returnvalue <- list(
    knns = t(knns),
    wij = wij,
    call = sys.call(),
    coords = coords
  )

  if (save_neighbors) {
    knns[knns == -1] <- NA
    returnvalue$knns <- t(knns)
  }
  if (save_edges) {
  	returnvalue$edges <- edges
  }

  class(returnvalue) <- "largeVis"
  return(returnvalue)
}
