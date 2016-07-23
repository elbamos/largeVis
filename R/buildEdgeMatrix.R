
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
                            verbose) {
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
#' @return A sparse matrix
#' @export
buildWijMatrix <- function(x,
										       perplexity) UseMethod("buildWijMatrix")
#' @export
#' @rdname buildWijMatrix
buildWijMatrix.TsparseMatrix <- function(x,
																	 perplexity) {
	wij <- referenceWij(x@j, x@i, x@x^2, perplexity)
	return(wij)
}
#' @export
#' @rdname buildWijMatrix
buildWijMatrix.CsparseMatrix <- function(x, perplexity) {
	is <- rep(0:(ncol(x) - 1), diff(x@p))
  wij <- referenceWij(is, x@i, x@x^2,perplexity)
  return(wij)
}
