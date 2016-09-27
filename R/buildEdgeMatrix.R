#' Build an nearest-neighbor graph weighted by distance.
#'
#' @param data A matrix with a number of columns equal to the number of columns in `x`
#' @param neighbors An adjacency matrix of the type produced by \code{\link{randomProjectionTreeSearch}}. If \code{NULL}, \code{\link{randomProjectionTreeSearch}}
#' will be run with parameters given by \code{...}.
#' @param distance_method One of "Euclidean" or "Cosine"
#' @param threads The number of threads to use in calculating distance; set automatically if \code{NULL} (the default).
#' @param verbose Verbosity
#' @param ... Additional parameters passed to \code{\link{randomProjectionTreeSearch}} if \code{neighbors} is \code{NULL}.
#'
#' @return A `sparseMatrix`, with the distance method stored in attribute \code{method} and additional class `edge_matrix.`
#' @importFrom Matrix sparseMatrix
#' @export
buildEdgeMatrix <- function(data,
                            neighbors = NULL,
                            distance_method = "Euclidean",
														threads = NULL,
                            verbose = getOption("verbose", TRUE),
														...) {
	if (is.null(neighbors)) neighbors <- randomProjectionTreeSearch(data, threads = threads, ...)
	indices <- neighborsToVectors(neighbors)
	distances <- distance(i = indices$i, j = indices$j, x = data, distance_method = distance_method, verbose = verbose)
	distances <- pmax(distances, 1e-5)
	mat <- sparseMatrix(
		                  i = indices$i + 1,
											j = indices$j + 1,
											x = as.vector(distances),
											dims = c(ncol(data), ncol(data)))
	structure(.Data = mat,
						method = tolower(distance_method))
}


#' edgeMatrixToDist
#'
#' Convert an edge matrix to a \code{dist} object.
#'
#' @param x An edge matrix.
#'
#' @return A \code{\link[stats]{dist}} object.
#'
#' @note This method converts the otherwise sparse edge matrix into a dense \code{dist} object,
#' where any distances absent from the edge matrix are represented as \code{NA}.
#'
#' @export
#' @rdname buildEdgeMatrix
#' @importFrom Matrix triu tril t
edgeMatrixToDist <- function(x) {
	y <- Matrix::tril(x)
	z <- Matrix::t(Matrix::triu(x))
	zeros <- y == 0
	y[zeros] <- z[zeros]
	y[y == 0] <- NA
	diag(y) <- 0
	structure(as.dist(as.matrix(y)),
						class = "dist",
						Size = ncol(x),
						Diag = FALSE,
						Upper = FALSE,
						method = attr(x, "method"),
						call = sys.call())
}

#' buildWijMatrix
#'
#' Rescale the weights in an edge matrix to match a given perplexity.
#'
#' @param x A sparse matrix
#' @param threads The maximum number of threads to spawn. Determined automatically if \code{NULL} (the default).
#' @param perplexity Given perplexity.
#'
#' @return A \code{list} with the following components: \describe{
#'    \item{'dist'}{An [N,K] matrix of the distances to the nearest neighbors.}
#'    \item{'id'}{An [N,K] matrix of the node indexes of the neartest neighbors.  Note that this matrix is 1-indexed,
#'    unlike most other matrices in this package.}
#'    \item{'k'}{The number of nearest neighbors.}
#'  }
#' @export
buildWijMatrix <- function(x,
													 threads = NULL,
										       perplexity = 50) UseMethod("buildWijMatrix")
#' @export
#' @rdname buildWijMatrix
buildWijMatrix.TsparseMatrix <- function(x,
																				 threads = NULL,
																	 perplexity = 50) {
	wij <- referenceWij(x@j, x@i, x@x^2, as.integer(threads), perplexity);
	return(wij)
}
#' @export
#' @rdname buildWijMatrix
buildWijMatrix.CsparseMatrix <- function(x, threads = NULL, perplexity = 50) {
	is <- rep(0:(ncol(x) - 1), diff(x@p))
  wij <- referenceWij(is, x@i, x@x^2, as.integer(threads), perplexity)
  return(wij)
}