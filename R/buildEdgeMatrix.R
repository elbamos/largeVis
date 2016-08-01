
#' Build an nearest-neighbor graph weighted by distance.
#'
#' @param x An edge matrix, 0-indexed, where each vertex is a column
#' @param data A matrix with a number of columns equal to the number of columns in `x`
#' @param distance_method One of "Euclidean" or "Cosine"
#' @param verbose Verbosity
#'
#' @return A `sparseMatrix`
#' @importFrom Matrix sparseMatrix
#' @export
buildEdgeMatrix <- function(data,
                            neighbors,
                            distance_method = "Euclidean",
                            verbose = TRUE) {
	indices <- neighborsToVectors(neighbors)
	distances <- distance(indices$i, indices$j, x = data, distance_method, verbose)
	mat <- Matrix::sparseMatrix(
		                  i = indices$i + 1,
											j = indices$j + 1,
											x = as.vector(distances),
											dims = c(ncol(data), ncol(data)))
	return(mat)
}

#' buildWijMatrix
#'
#' Rescale the weights in an edge matrix to match a given perplexity.
#'
#' @param x A sparse matrix
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
										       perplexity) UseMethod("buildWijMatrix")
#' @export
#' @rdname buildWijMatrix
buildWijMatrix.TsparseMatrix <- function(x,
																	 perplexity = 50) {
	wij <- referenceWij(x@j, x@i, x@x^2, perplexity)
	return(wij)
}
#' @export
#' @rdname buildWijMatrix
buildWijMatrix.CsparseMatrix <- function(x, perplexity = 50) {
	is <- rep(0:(ncol(x) - 1), diff(x@p))
  wij <- referenceWij(is, x@i, x@x^2,perplexity)
  return(wij)
}
