#' lv_dbscan
#'
#' Implementation of the DBSCAN algorithm using largeVis datastructures.
#'
#' @param edges A weighted graph of the type produced by \code{\link{buildEdgeMatrix}}. Alternatively, a \code{largeVis} object,
#' in which case \code{edges} and \code{neighbors} will be taken from the \code{edges} and \code{knns} parameters, respectively.
#' @param neighbors An adjacency matrix of the type produced by \code{\link{randomProjectionTreeSearch}}
#' @param eps See \code{\link[dbscan]{dbscan}}.
#' @param minPts See \code{\link[dbscan]{dbscan}}.
#' @param verbose Vebosity level.
#'
#' @details The DBSCAN algorithm attempts to find clusters of a minimum density given by \code{eps}. This
#' implementation leverages the nearest neighbor data assembled by largeVis.
#'
#' @return A \code{\link[dbscan]{dbscan}} object.
#' @export
#'
#' @references Martin Ester, Hans-Peter Kriegel, Jorg Sander, Xiaowei Xu (1996). Evangelos Simoudis, Jiawei Han, Usama M. Fayyad, eds. A density-based algorithm for discovering clusters in large spatial databases with noise. Proceedings of the Second International Conference on Knowledge Discovery and Data Mining (KDD-96). AAAI Press. pp. 226–231. ISBN 1-57735-004-9.
lv_dbscan <- function(edges,
									 neighbors,
									 eps = Inf,
									 minPts = nrow(neighbors - 1),
									 verbose = getOption("verbose", TRUE)) {
	if (inherits(edges, "largeVis")) {
		if (missing(neighbors)) neighbors <- edges$knns
		edges <- edges$edges
	}
	if (!is.null(neighbors)) {
		neighbors[is.na(neighbors)] <- -1
		if (ncol(neighbors) != ncol(edges)) neighbors <- t(neighbors)
	}
	if (is.null(edges) || is.null(neighbors)) stop("Both edges and neighbors must be provided.")

	clusters <- dbscan_cpp(edges, neighbors, as.double(eps), as.integer(minPts), as.logical(verbose))

	structure(list(cluster = clusters, eps = eps, minPts = minPts, call = sys.call()),
						class = c("dbscan_fast", "dbscan"))
}

#' @title Local Outlier Factor Score
#'
#' @description Calculate the Local Outlier Factor (LOF) score for each data point given knowledge
#' of k-Nearest Neighbors.
#'
#' @param edges An edge matrix of the type produced by \code{\link{buildEdgeMatrix}}.
#'
#' @references Based on code in the \code{\link[dbscan]{dbscan}} package.
#'
#' @return A vector of LOF values for each data point.
#' @export
lof <- function(edges) {

	id <- apply(edges,MARGIN = 1, FUN = function(x) which(x != 0))
	dist <- apply(edges, MARGIN = 1, FUN = function(x) x[x != 0])
	for (i in 1:ncol(id)) {
		ord <- order(dist[,i])
		id[,i] <- id[,i][ord]
		dist[,i] <- dist[,i][ord]
	}
	K <- nrow(id)
	N <- ncol(id)
	dist <- t(dist)
	id <- t(id)

	lrd <- rep(0, N)

	for (i in 1:N) {
		merged <- cbind(dist[id[i,], K], dist[i, ])
		lrd[i] <- 1/(sum(apply(merged, 1, max))/K)
	}

	ret <- rep(0, N)
	for (i in 1:N) ret[i] <- sum(lrd[id[i,]])/K / lrd[i]

	ret[is.nan(ret)] <- NA

  ret
}

