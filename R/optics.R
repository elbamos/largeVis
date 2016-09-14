#' optics
#'
#' Experimental implementation of the OPTICS algorithm.
#'
#' @param edges A weighted graph of the type produced by \code{\link{buildEdgeMatrix}}.
#' @param neighbors An adjacency matrix of the type produced by \code{\link{randomProjectionTreeSearch}}
#' @param eps See \code{\link[dbscan]{optics}}.
#' @param minPts See \code{\link[dbscan]{optics}}.
#' @param eps_cl See \code{\link[dbscan]{optics}}.
#' @param xi See \code{\link[dbscan]{optics}}.
#' @param verbose Vebosity level.
#'
#' @details This is a preliminary implementation of a variant of the OPTICS algorithm that attempts
#' to leverage the \code{largeVis} nearest-neighbor search.
#'
#'
#' @note Support for dbscan and optics are preliminary, and not fully tested for
#' correctness.
#'
#' @note This is not the original OPTICS algorithm. In particular, the neighbor-search strategy in
#' OPTICS is not used, in favor of using the pre-calculated neighbor matrix produced incidentally by
#' `largeVis`.
#'
#' @return An \code{\link[dbscan]{optics}} object.
#' @references  Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Jörg Sander (1999). OPTICS: Ordering Points To Identify the Clustering Structure. ACM SIGMOD international conference on Management of data. ACM Press. pp. 49–60.
#' @export
#' @importFrom dbscan optics_cut opticsXi
optics <- function(edges,
									 neighbors,
									 eps = Inf,
									 minPts = nrow(neighbors),
									 eps_cl,
									 xi,
									 verbose = getOption("verbose", TRUE)) {
	if (is.null(edges) || is.null(neighbors)) stop("Both edges and neighbors must be specified.")
	ret <- optics_cpp(edges = edges, neighbors = neighbors, eps = as.double(eps), minPts = as.integer(minPts), verbose = as.logical(verbose))

	ret$minPts <- minPts
	ret$eps <- eps
	ret$eps_cl <- NA
	ret$call <- sys.call()
	class(ret) <- "optics"

	if(!missing(eps_cl)) ret <- dbscan::optics_cut(ret, eps_cl)
	if(!missing(xi)) ret <- dbscan::opticsXi(ret, xi)

	ret
}
