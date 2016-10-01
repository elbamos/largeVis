#' lv_optics
#'
#' Experimental implementation of the OPTICS algorithm.
#'
#' @param edges A weighted graph of the type produced by \code{\link{buildEdgeMatrix}}. Alternatively, a \code{largeVis} object,
#' in which case \code{edges} and \code{neighbors} will be taken from the \code{edges} and \code{knns} parameters, respectively.
#' @param neighbors An adjacency matrix of the type produced by \code{\link{randomProjectionTreeSearch}}
#' @param eps See \code{\link[dbscan]{optics}}.
#' @param minPts See \code{\link[dbscan]{optics}}.
#' @param eps_cl See \code{\link[dbscan]{optics}}.
#' @param xi See \code{\link[dbscan]{optics}}.
#' @param useQueue Whether to process points in order of core distznce.  (See note.)
#' @param verbose Vebosity level.
#'
#' @details This is a preliminary implementation of a variant of the OPTICS algorithm that attempts
#' to leverage the \code{largeVis} nearest-neighbor search.
#'
#'
#' @note Support for dbscan and optics are preliminary, and not fully tested for
#' correctness.
#'
#' @note the \code{useQueue} parameter controls the order in which points that have not yet been visisted are processed. If \code{FALSE},
#' points are processed in order of rows. If \code{TRUE}, they are processed in ascending order of core distance. \code{FALSE} is more
#' compatible with the implementations in the \code{dbscan} package and in the \code{ELKI} Java clustering package. \code{TRUE} may produce
#' preferrable results.
#'
#' @note This is not the original OPTICS algorithm. In particular, the neighbor-search strategy in
#' OPTICS is not used, in favor of using the pre-calculated neighbor matrix produced incidentally by
#' `largeVis`.
#'
#' @return An \code{\link[dbscan]{optics}} object.
#' @references  Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Jorg Sander (1999). OPTICS: Ordering Points To Identify the Clustering Structure. ACM SIGMOD international conference on Management of data. ACM Press. pp. 49–60.
#' @export
#' @importFrom dbscan optics_cut opticsXi
lv_optics <- function(edges,
										 neighbors,
										 eps = Inf,
										 minPts = nrow(neighbors),
										 eps_cl,
										 xi,
										 useQueue = TRUE,
										 verbose = getOption("verbose", TRUE)) {
	if (inherits(edges, "largeVis")) {
		if (missing(neighbors)) neighbors <- edges$knns
		edges <- edges$edges
	}
	if (!is.null(neighbors)) {
		neighbors[is.na(neighbors)] <- -1
		if (ncol(neighbors) != ncol(edges)) neighbors <- t(neighbors)
	}
	if (is.null(edges) || is.null(neighbors)) stop("Both edges and neighbors must be specified.")
	ret <- optics_cpp(edges = edges,
										neighbors = neighbors,
										eps = as.double(eps),
										minPts = as.integer(minPts),
										useQueue = as.logical(useQueue),
										verbose = as.logical(verbose))

	ret$minPts <- minPts
	ret$eps <- eps
	ret$eps_cl <- NA
	ret$call <- sys.call()
	class(ret) <- "optics"

	if(!missing(eps_cl)) ret <- dbscan::optics_cut(ret, eps_cl)
	if(!missing(xi)) ret <- dbscan::opticsXi(ret, xi)

	ret
}
