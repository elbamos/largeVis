#' OPTICS
#'
#' Experimental implementation of the OPTICS algorithm.
#'
#' @param data Input data, where examples are columns.
#' @param neighbors An adjacency matrix of the type produced by \code{\link{randomProjectionTreeSearch}}
#' @param edges A weighted graph of the type produced by \code{\link{buildEdgeMatrix}}.
#' @param eps See \code{\link[dbscan]{optics}}.
#' @param minPts See \code{\link[dbscan]{optics}}.
#' @param eps_cl See \code{\link[dbscan]{optics}}.
#' @param xi See \code{\link[dbscan]{optics}}.
#' @param verbose Vebosity level.
#'
#' @details This is a preliminary implementation of a variant of the OPTICS algorithm that attempts
#' to leverage the \code{largeVis} nearest-neighbor search.
#'
#' One of \code{neighbors} or \code{edges} must be specified. If \code{edges} is missing,
#' \code{data} must also be given. If \code{data} is given along with either \code{edges}
#' or \code{neighbors}, the algorithm will attempt a more thorough search.
#'
#' @note Support for dbscan and optics are preliminary, and not fully tested for
#' correctness.
#'
#' @note This is not the original OPTICS algorithm. In particular, the neighbor-search strategy in
#' OPTICS is not used, in favor of using a pre-calculated neighbor matrix produced incidentally by
#' `largeVis`.
#'
#' @return An \code{\link[dbscan]{optics}} object.
#' @export
#' @importFrom dbscan optics_cut opticsXi
optics <- function(data = NULL,
                   neighbors = NULL,
                   edges = NULL,
                   eps,
                   minPts = nrow(data) + 1,
                   eps_cl,
                   xi,
                   verbose = getOption("verbose", TRUE)) {
  if (! is.null(edges) && is.null(data))
    ret <- optics_e(edges = edges,
                    eps = as.double(eps), minPts = as.integer(minPts),
                    verbose = verbose)
  else if (! is.null(edges))
    ret <- optics_ed(edges = edges, data = data,
                     eps = as.double(eps), minPts = as.integer(minPts),
                     verbose = verbose)
  else
    ret <- optics_nd(neighbors = neighbors, data = data,
                     eps = as.double(eps), minPts = as.integer(minPts),
                     verbose = verbose)

  ret$minPts <- minPts
  ret$eps <- eps
  ret$eps_cl <- NA
  class(ret) <- "optics"

  if(!missing(eps_cl)) ret <- dbscan::optics_cut(ret, eps_cl)
  if(!missing(xi)) ret <- dbscan::opticsXi(ret, xi)

  ret
}

edgeMatrixToKNNS <- function(edges) {
  id = apply(edges,MARGIN = 1, FUN = function(x) which(x != 0))
  dist = apply(edges, MARGIN = 1, FUN = function(x) x[x != 0])
  for (i in 1:ncol(id)) {
    ord <- order(dist[,i])
    id[,i] <- id[,i][ord]
    dist[,i] <- dist[,i][ord]
  }
  k = nrow(id)
  list(dist = t(dist), id = t(id), k = k)
}

# The source code for function lof is based on code that bore this license:
#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2015 Michael Hahsler

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


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
  kNNlist <- edgeMatrixToKNNS(edges)
  N <- nrow(kNNlist$id)
  K <- kNNlist$k

  # lrd <- rep(0, N)
  lrd <- rep(0, N)
  # for(i in 1:N) {
  #   input <- kNNlist$dist[c(i, kNNlist$id[i, ]) ,]
  #   lrd[i] <- 1 / (sum(apply(input, MARGIN = 1, max)) / K)
  # }
  for(i in 1:N) lrd[i] <- 1/(sum(apply(
    cbind(kNNlist$dist[kNNlist$id[i,], K], kNNlist$dist[i,]),
    1, max)) / K)

  ret <- rep(0, N)
  for (i in 1:N) ret[i] <- sum(lrd[kNNlist$id[i,]])/K / lrd[i]

  ret[is.nan(ret)] <- NA

  ret
}

#' hdbscan
#'
#' Implemenation of the hdbscan algorithm.
#'
#' @param edges An edge matrix of the type returned by \code{\link{buildEdgeMatrix}}.
#' @param minPts The minimum number of points in a cluster.
#' @param K The number of points in the core neighborhood. (See details.)
#' @param neighbors An adjacency matrix of the type returned by \code{\link{randomProjectionTreeSearch}}.
#' @param threads The maximum number of threads to spawn. Determined automatically if \code{NULL} (the default).
#' @param verbose Verbosity.
#'
#' @details The hyperparameter \code{K} controls the size of core neighborhoods.
#' The algorithm measures the density around a point as 1 / the distance between
#' that point and its \eqn{k}th nearest neighbor. A low value of \code{K} is similar
#' to clustering nearest neighbors rather than based on density. A high value of
#' \code{K} may cause the algorithm to miss some (usually contrived) clustering
#' patterns, such as where clusters are made up of points arranged in lines to form
#' shapes.
#'
#' If \code{neighbors} is specified, some costly sorts of neighbors in the edge
#' matrix may be avoided.
#'
#' The function must be provided sufficient nearest-neighbor data for whatever
#' is specified for \eqn{k}. If \eqn{k} = 5, for example, the edge matrix (and
#' neighbor matrix, if specified) must contain data on at least 5 neighbors for
#' each point. This should not be problematic in typical use in connection with
#' \code{\link{largeVis}}, which is ordinarily run with a far higher \eqn{k}-value
#' than hdbscan.
#'
#' @return An object of type \code{hdbscan} with the following fields:
#' \describe{
#'    \item{'clusters'}{A vector of the cluster membership for each vertex. Outliers
#'    are given \code{NA}}
#'    \item{'probabilities'}{A vector of the degree of each vertex' membership. This
#'    is calculated as each vertex' \eqn{\lambda_p} over the highest \eqn{\lambda_p}
#'    in the cluster. }
#'    \item{'tree'}{The minimum spanning tree used to generate the clustering.}
#'    \item{'hierarchy'}{A representation of the condensed cluster hierarchy.}
#'    \item{'call'}{The call.}
#'  }
#'
#'  The hierarchy describes the complete post-condensation structure of the tree:
#'  \describe{
#'  \item{'nodemembership'}{The node ID of the vertex's immediate parent, after condensation.}
#'  \item{'lambda'}{\eqn{\lambda_p}}
#'  \item{'parent'}{The node ID of each node's parent.}
#'  \item{'stability'}{The node's stability, taking into account child-node stabilities.}
#'  \item{'selected'}{Whether the node was selected.}
#'  \item{'coredistances'}{The core distance determined for each vertex.}
#'  }
#'
#' @references R. Campello, D. Moulavi, and J. Sander, Density-Based Clustering Based on Hierarchical Density Estimates In: Advances in Knowledge Discovery and Data Mining, Springer, pp 160-172. 2013
#' @seealso \url{https://github.com/lmcinnes/hdbscan}
#' @note This is not precisely the \code{HDBSCAN} algorithm because it relies on the
#' nearest neighbor data generated by the \code{LargeVis} algorithm. In particular,
#' \code{HDBSCAN} assumes that all points can be connected in a single minimal-spanning
#' tree. This implementation uses a minimal-spanning forest, because the graph may not
#' be fully connected depending on the amount of nearest-neighbor data provided.
#' If, for example, the data has distinct clusters where no member of some cluster is a
#' nearest neighbor of a point in any other cluster, which can easily happen, the algorithm will
#' generate distinct trees for each disconnected set of points. In testing, this
#' improved the performance of the algorithm.
#'
#' @note Do not rely on the content of the \code{probabilities} field for outliers. A future version
#' will hopefully provide some measure of outlier factor.
#' @examples
#' \dontrun{
#' library(largeVis)
#' library(clusteringdatasets)
#' data(spiral)
#' dat <- as.matrix(spiral[, 1:2])
#' neighbors <- randomProjectionTreeSearch(t(dat), K = 10, tree_threshold = 100,
#'                                        max_iter = 5)
#' edges <- buildEdgeMatrix(t(dat), neighbors)
#' clusters <- hdbscan(edges, verbose = FALSE)
#'
#' largeHighDimensionalDataset <- matrix(rnorm(50000), ncol = 50)
#' vis <- largeVis(largeHighDimensionalDataset)
#' edges <- buildEdgeMatrix(largeHighDimensionalDataset, vis$knns)
#' clustering <- hdbscan(edges)
#' gplot(clustering, t(vis$coords))
#' }
#' @export
#' @importFrom stats aggregate
hdbscan <- function(edges, minPts = 20, K = 5, neighbors = NULL,
										threads = NULL,
                    verbose = getOption("verbose", TRUE)) {

	if (! is.null(neighbors)) {
		neighbors[is.na(neighbors)] <- -1
		if (ncol(neighbors) != ncol(edges)) neighbors <- t(neighbors)
	}
  clustersout <- hdbscanc(edges, neighbors, K, minPts, threads, verbose)
  clusters <- clustersout$clusters[1, ]
  clusters[clusters == -1] <- NA
  clusters = factor(clusters, exclude = NULL)
  probs <- data.frame(
    probs = clustersout$clusters[2, ]
  )
  mins = stats::aggregate(probs, by = list(clusters), FUN = "min")$probs
  maxes = stats::aggregate(probs, by = list(clusters), FUN = "max")$probs - mins
  probs$probs[!is.na(clusters)] <- (probs$probs[!is.na(clusters)] -
                                    mins[as.integer(clusters)[!is.na(clusters)]]) /
    maxes[as.integer(clusters)[!is.na(clusters)]]

  ret <- list(
    clusters = clusters,
    probabilities = probs$probs,
    tree = clustersout$tree,
    hierarchy = clustersout$hierarchy,
    call = sys.call()
  )
  class(ret) <- "hdbscan"
  return(ret)
}

#' gplot
#'
#' Plot an \code{hdbscan} object, using \code{\link[ggplot2]{ggplot}}. The
#' plot is primarily intended for diagnostic purposes, but can be useful to undertand
#' how clusters were generated.
#'
#' Point color corresponds to clusters, with outliers as the \code{NA} color. Alpha
#' corresponds to the centrality of the node in the cluster (i.e., \eqn{\lambda_p} relative
#' to \eqn{\lambda_{birth}} and \eqn{\lambda_{death}}). The segments on the plot
#' correspond to the connections on the minimum spanning tree. Segment alpha
#' corresponds to \eqn{\lambda_p}.
#'
#' If the parameter \code{text} is set to \code{TRUE} or \code{"parent"}, the nodes will
#' be labelled with the node index number or cluster index number, respectively.
#'
#' @param x An \code{hdbscan} object.
#' @param coords Coordinates for the points clustered in \code{x}.
#' @param text If \code{TRUE}, include on the plot labels for each node's index.
#' If \code{"parent"}, the labels will instead be the index number of the node's
#' cluster.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#' @export
#'
#' @examples
#' \dontrun{
#' library(largeVis)
#' library(clusteringdatasets)
#' data(Aggregation)
#' dat <- as.matrix(Aggregation[, 1:2])
#' neighbors <- randomProjectionTreeSearch(t(dat), K = 10, tree_threshold = 100,
#'                                        max_iter = 5)
#' edges <- buildEdgeMatrix(t(dat), neighbors)
#' clusters <- hdbscan(edges, verbose = FALSE)
#' plot(clusters, dat)
#' }
#' @importFrom ggplot2 ggplot unit geom_label geom_point geom_segment aes_
gplot <- function(x, coords, text = FALSE) {
  dframe <- data.frame(coords)
  colnames(dframe) <- c("x", "y")
  dframe$cluster = x$clusters
  dframe$probabilities = x$probabilities
  dframe$probabilities[is.nan(dframe$probabilities)] <-
    x$hierarchy$lambda[is.nan(dframe$probabilities)]
  tree <- x$tree
  tree[tree == -1] <- NA
  xy <- data.frame(coords[tree + 1, ])
  colnames(xy) <- c("x2", "y2")
  dframe <- cbind(dframe, xy)
  dframe$lambda <- x$hierarchy$lambda / max(x$hierarchy$lambda)
  dframe$label <- 0:(nrow(dframe) - 1)
  dframe$parent <- x$hierarchy$nodemembership
  plt <- ggplot2::ggplot(dframe,
                         ggplot2::aes_(x = quote(x), y = quote(y),
                    xend = quote(x2), yend = quote(y2), color = quote(cluster))) +
    ggplot2::geom_point(aes_(alpha = quote(probabilities)), size = 0.7) +
    ggplot2::geom_segment(size = 0.5, ggplot2::aes_(alpha = quote(lambda), size = quote(lambda)))
  if (text == "parent") {
    plt <- plt + ggplot2::geom_label(ggplot2::aes_(label = quote(parent)), size = 2.5,
                            label.padding = ggplot2::unit(0.1, "lines"),
                            label.size = 0.1, alpha = 0.7)
  } else if (text) {
    plt <- plt + ggplot2::geom_label(ggplot2::aes_(label = quote(label)), size = 2.5,
                            label.padding = ggplot2::unit(0.1, "lines"),
                            label.size = 0.1, alpha = 0.7)
  }
  plt
}