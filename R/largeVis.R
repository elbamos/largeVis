#' Apply the LargeVis algorithm for visualizing large high-dimensional datasets.
#'
#' @param x A matrix, where the features are rows and the examples are columns.
#' @param dim The number of dimensions in the output
#' @param K The number of nearest-neighbors to use in computing the kNN graph
#' @param n_trees See \code{\link{randomProjectionTreeSearch}}.  The default is set at 50, which is the number
#' used in the examples in the original paper.
#' @param max_iter See \code{\link{randomProjectionTreeSearch}}.
#' @param distance_method See \code{\link{randomProjectionTreeSearch}}.
#' @param perplexity See \code{\link{buildWijMatrix}}.
#' @param save_neighbors Whether to include in the output the adjacency matrix of nearest neighbors.
#' @param save_edges Whether to include in the output the distance matrix of nearest neighbors.
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
#' # mnist
#' # Note: The MNIST dataset may be obtained using the deepnet package.
#' load("./mnist.Rda")
#' dat <- mnist$train$x
#' dat <- (dat / 255) - 0.5
#' dat <- t(dat)
#' visObject <- largeVis(dat, n_trees = 50, K = 50)
#' plot(t(visObject$coords))
#' }
#'
largeVis <- function(x,
                     dim = 2,
                     K = 50,

                     n_trees = 50,
                     max_iter = 1,
                     distance_method = "Euclidean",

                     perplexity = max(50, K / 3),

                     save_neighbors = TRUE,
										 save_edges = TRUE,

                     verbose = getOption("verbose", TRUE),
                    ...) {

	if (!(is.matrix(x) && is.numeric(x)) && !is.data.frame(x) && ! inherits(x, "Matrix")) stop("LargeVis requires a matrix or data.frame")
	if (is.data.frame(x)) x <- t(as.matrix(x[, sapply(x, is.numeric)]))

  #############################################
  # Search for kNearestNeighbors
  #############################################
  knns <- randomProjectionTreeSearch(x,
                                     n_trees = n_trees,
                                     K = K,
                                     max_iter = max_iter,
                                     distance_method = distance_method,
                                     verbose = verbose)

  if (!save_neighbors) knns$neighbors <- NULL
  gc()
  if (any(knns$edgematrix$x > 27)) {
  	warning(paste(
  		"The Distances between some neighbors are large enough to cause the calculation of p_{j|i} to overflow.",
  		"Scaling the distance vector."))
  	knns$edgematrix$x <- knns$edgematrix$x / max(knns$edgematrix$x)
  }
  wij <- buildWijMatrix(knns$edgematrix, perplexity)
  if (!save_edges && !save_neighbors) rm(knns)

  #######################################################
  # Estimate embeddings
  #######################################################
  coords <- projectKNNs(wij = wij,
                        dim = dim,
                        verbose = verbose,
                        ...)

  #######################################################
  # Cleanup
  #######################################################

  returnvalue <- structure(list(
    wij = wij,
    call = sys.call(),
    coords = coords
  ), class = "largeVis")

  if (save_neighbors) {
    knns$neighbors[knns$neighbors == -1] <- NA
    returnvalue$knns <- t(knns$neighbors)
  }
  if (save_edges) {
  	returnvalue$edges <- knns$edgematrix
  }

  return(returnvalue)
}
