#' @importFrom Matrix sparseMatrix
toMatrix <- function(x) {
	sparseMatrix(
		i = x$j,
		j = x$i,
		x = x$x,
		dims = attr(x, "dims")
	)
}

#' as.dist.edgematrix
#'
#' Convert an edge matrix to a \code{dist} object.
#'
#' @param m An `edgematrix` object.
#' @param diag logical value indicating whether the diagonal of the distance matrix should be printed by print.dist.
#' @param upper	logical value indicating whether the upper triangle of the distance matrix should be printed by print.dist.
#'
#' @return A \code{\link[stats]{dist}} object.
#'
#' @note This method converts the otherwise sparse edge matrix into a dense \code{dist} object,
#' where any distances absent from the edge matrix are represented as \code{NA}.
#'
#' @export
#' @rdname buildEdgeMatrix
#' @importFrom Matrix triu tril t as.matrix diag sparseMatrix
#' @importFrom stats as.dist
as.dist.edgematrix <- function(m, diag = FALSE, upper = FALSE) {
	x <- toMatrix(m)
	y <- Matrix::tril(x)
	z <- Matrix::t(Matrix::triu(x))
	zeros <- y == 0
	y[zeros] <- z[zeros]
	y[y == 0] <- NA
	Matrix::diag(y) <- 0
	structure(stats::as.dist(Matrix::as.matrix(y)),
						class = c("dist", "dissimilarity"),
						Size = ncol(x),
						Diag = FALSE,
						Upper = FALSE,
						Metric = attr(m, "Metric"),
						NA.message = "Sparse dist matrix",
						Diag = diag,
						Upper = upper,
						call = attr(m, "call"))
}

#' buildWijMatrix
#'
#' Rescale the weights in an edge matrix to match a given perplexity.
#'
#' @param x An edgematrix, either an `edgematrix` object or a sparse matrix.
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
										       perplexity = 50) UseMethod("buildWijMatrix")
#' @export
#' @rdname buildWijMatrix
buildWijMatrix.edgematrix <- function(x,
																		 perplexity = 50) {
	buildWijMatrix(toMatrix(x),  perplexity)
}
#' @export
#' @rdname buildWijMatrix
buildWijMatrix.TsparseMatrix <- function(x,
																	       perplexity = 50) {
	wij <- referenceWij(x@j, x@i, x@x^2, perplexity);
	return(wij)
}
#' @export
#' @rdname buildWijMatrix
buildWijMatrix.CsparseMatrix <- function(x, perplexity = 50) {
	is <- rep(0:(ncol(x) - 1), diff(x@p))
  wij <- referenceWij(is, x@i, x@x^2, perplexity)
  return(wij)
}