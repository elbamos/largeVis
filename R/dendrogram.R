#' @title as_dendrogram_hdbscan
#' @description Convert an hdbscan object into dendrogram compatible with the \code{stats} package.
#' @note The hdbscan algorithm works by first building a hierarchy based on a minimal spanning tree, then consolidating nodes according to
#' rules in the algorithm. The algorithm then selects some of the consolidated nodes as clusters, deselecting others. For example, if Node A has children
#' B and C, the algorithm might select A, and then all points under A, B, and C would be assigned to the same cluster. Or, it might deselect A,
#' and select B and C instead. In that case, the nodes under B would be assigned to one cluster, the nodes under C to a different cluster, and nodes
#' under A but not B or C would not be assigned to any cluster. This function returns a dendrogram of the middle stage, the hierarchy of consolidated
#' nodes. Whether a node was selected as as cluster is an attribute of each node.
#' @param object An \code{hdbscan} object.
#' @return A \code{dendrogram} object, where nodes have the following attributes:
#' \describe{
#' \item{'leaf'}{As in \code{\link[stats]{dendrogram}}.}
#' \item{'members'}{As in \code{\link[stats]{dendrogram}}.}
#' \item{'height'}{The \eqn{\lamba_{birth}} of the node or leaf.}
#' \item{'probability'}{The probability that the leaf is a true member of its assigned cluster.}
#' \item{'stability'}{The node's determined stability, taking into account child-node stabilities. Missing for leaves.}
#' \item{'selected'}{Whether the node was selected as a cluster.  Missing for leaves.  Note that when a node is selected,
#' all points under child branches are assigned to the same cluster.}
#' \item{'cluster'}{The cluster number, for reference against the \code{hdbscan} object.}
#' }
#' @export
#' @examples
#' data(iris)
#' neighbors <- randomProjectionTreeSearch(t(iris[,1:4]), K = 20)
#' edges <- buildEdgeMatrix(t(iris[,1:4]), neighbors)
#' hdbscanobj <- hdbscan(edges, minPts = 10, K = 5)
#' plot(as_dendrogram_hdbscan(hdbscanobj))
as_dendrogram_hdbscan <- function(object) {
	C <- length(object$hierarchy$parent)
	counts <- tabulate(object$hierarchy$nodemembership, nbins = C) + tabulate(object$hierarchy$parent, nbins = C)

	leafs <- lapply(1:length(object$hierarchy$nodemembership), FUN = function(i)
		structure(as.list(i),
			leaf = TRUE,
			members = 1L,
			label = i,
			cluster = object$clusters[i],
			probability = object$probabilities[i],
			height = object$hierarchy$lambda[i],
			class = "dendrogram"
		))

	clusters <- vector("list", C)
	h <- function(z) attr(z, "height")
	m <- function(z) attr(z, "members")
	for (i in C:1) {
		children <- which(object$hierarchy$nodemembership == i)
		children <- leafs[children]
		cousins <- which(object$hierarchy$parent == i & (1:C != i))
		cousins <- clusters[cousins]
		cousins <- cousins[!sapply(cousins, is.null)]
		newcluster <- structure(
			.Data = c(children, cousins),
			members = sum(unlist(vapply(children, FUN.VALUE = 0L, FUN = m))) +
				        sum(unlist(vapply(cousins, FUN.VALUE = 0L, FUN = m))),
			leaf = FALSE,
			height = 1.1 * max(0, unlist(vapply(children, FUN.VALUE = 0, FUN = h)),
												    unlist(vapply(cousins, FUN.VALUE = 0, FUN = h))),
			selected = object$hierarchy$selected[i],
			cluster = i,
			stability = object$hierarchy$stability[i],
			label = paste("cluster", i, "stability",object$hierarchy$stability[i], ifelse(object$hierarchy$selected[i], "selected", "")),
			class = "dendrogram"
		)
		if (attr(newcluster, "members") > 0) {
			clusters[[i]] <- newcluster
		}
		clusters[(object$hierarchy$parent == i) & (1:C != i)] <- NA
	}

	clusters[which(is.na(clusters))] <- NULL
	clusters <- clusters[!sapply(clusters, is.null)]

	if (length(clusters) == 1) clusters[[1]]
	else {
		structure(clusters,
			members = sum(unlist(vapply(clusters, FUN.VALUE = 0L, FUN = m))),
			leaf = FALSE,
			height = 1.1 * max(unlist(vapply(clusters, FUN.VALUE = 0, FUN = h))),
			selected = 0,
			cluster = 0,
			stability = Inf,
			label = "root",
			class = "dendrogram"
		)
	}
}
