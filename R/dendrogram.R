#' @title as.dendrogram.hdbscan
#' @description Convert an hdbscan object into a dendrogram compatible with the \code{stats} package.
#' @note The hdbscan algorithm works by first building a hierarchy based on a minimal spanning tree, then consolidating nodes according to
#' rules in the algorithm. The algorithm then selects some of the consolidated nodes as clusters, deselecting others. For example, if Node A has children
#' B and C, the algorithm might select A, and then all points under A, B, and C would be assigned to the same cluster. Or, it might deselect A,
#' and select B and C instead. In that case, the nodes under B would be assigned to one cluster, the nodes under C to a different cluster, and nodes
#' under A but not B or C would not be assigned to any cluster. This function returns a dendrogram of the middle stage, the hierarchy of consolidated
#' nodes. Whether a node was selected as as cluster is an attribute of each node.
#' @param object An \code{hdbscan} object.
#' @param includeNodes Whether individual nodes should be included in the dedrogram. Can cause a substantial increase in the size of the object.
#' @param ... For compatibility with \code{\link[stats]{as.dendrogram}}, and currently ignored.
#' @return A \code{dendrogram} object, where nodes have the following attributes:
#' \describe{
#' \item{'leaf'}{As in \code{\link[stats]{dendrogram}}.}
#' \item{'members'}{As in \code{\link[stats]{dendrogram}}.}
#' \item{'size'}{The number of nodes underneath the cluster.}
#' \item{'height'}{The core distance at which the cluster or node was merged.}
#' \item{'probability'}{The probability that the leaf is a true member of its assigned cluster.}
#' \item{'GLOSH'}{The leaf's GLOSH outlier score.}
#' \item{'stability'}{The node's determined stability, taking into account child-node stabilities. Missing for leaves.}
#' \item{'selected'}{Whether the node was selected as a cluster.  Missing for leaves.  Note that when a node is selected,
#' all points under child branches are assigned to the same cluster.}
#' \item{'cluster'}{The cluster number, for reference against the \code{hdbscan} object.}
#' }
#'
#' @note This function remains experimental in terms of finding the best way to represent an hdbscan object in a dendrogram.
#'
#' @importFrom stats as.dendrogram
#' @export
#' @examples
#' \dontrun{
#' data(iris)
#' vis <- largeVis(t(iris[,1:4]), K = 20, sgd_batches = 1)
#' hdbscanobj <- hdbscan(vis, minPts = 10, K = 5)
#' plot(as_dendrogram_hdbscan(hdbscanobj))
#' }
as.dendrogram.hdbscan <- function(object, includeNodes = FALSE, ...) {
	C <- length(object$hierarchy$parent)

	leafs <- NULL
	if (includeNodes) leafs <- lapply(1:length(object$hierarchy$nodemembership), FUN = function(i)
		structure(as.list(i),
							leaf = TRUE,
							members = 1L,
							label = i,
							GLOSH = object$glosh[i],
							cluster = object$clusters[i],
							probability = object$probabilities[i],
							midpoint = object$probabilities[i] - 0.5,
							height = 1 / object$hierarchy$lambda[i],
							class = "dendrogram"
		))

	clusters <- vector("list", C)
	h <- function(z) attr(z, "height")
	m <- function(z) attr(z, "members")
	for (i in C:1) {
		children <- NULL
		extraChildren <- 0
		whichcousins <- which(object$hierarchy$parent == i & (1:C != i))
		cousins <- clusters[whichcousins]
		cousins <- cousins[!sapply(cousins, is.null)]

		if (includeNodes) {
			children <- which(object$hierarchy$nodemembership == i)
			children <- leafs[children]
		} else {
			cnt <- sum(object$hierarchy$nodemembership == i)
			if (cnt > 0) children <- list(
				structure(
					.Data = as.list(cnt),
					leaf = TRUE,
					members = 1L,
					label = paste(cnt, "Fallen Points"),
					height = 1 / max(object$hierarchy$lambda[object$hierarchy$nodemembership == i], na.rm = TRUE),
					class = "dendrogram"
				)
			)
			extraChildren <- extraChildren + cnt
		}

		members <- sum(unlist(vapply(children, FUN.VALUE = 0L, FUN = m))) +
							 sum(unlist(vapply(cousins, FUN.VALUE = 0L, FUN = m)))
		height <- max(0, unlist(vapply(children, FUN.VALUE = 0, FUN = h)),
					           unlist(vapply(cousins, FUN.VALUE = 0, FUN = h)),
									   max(1 / object$hierarchy$lambda[object$hierarchy$nodemembership == i], na.rm = TRUE))
		midpoint <- 0
		if (!is.na(object$hierarchy$parent[i])) {
			midpoint <- ifelse(i == min(which(object$hierarchy$parent == object$hierarchy$parent[i])), 0, 1)
		}
		lb <- object$hierarchy$lamba_birth[i]

		childnodes <- c(children,cousins)
		heights <- unlist(vapply(childnodes, FUN.VALUE = 0, FUN = h))
		childnodes <- childnodes[order(heights)]

		newcluster <- structure(
			.Data = childnodes,
			members = members,
			leaf = ifelse(members == 0, TRUE, FALSE),
			size = members + extraChildren,
			height = height,
			midpoint = midpoint,
			selected = object$hierarchy$selected[i],
			cluster = i,
			lambda_birth = lb,
			lambda_death = object$hierarchy$lambda_death[i],
			stability = object$hierarchy$stability[i],
			label = paste("cluster", i, "stability",signif(object$hierarchy$stability[i], 4), ifelse(object$hierarchy$selected[i], "selected", "")),
			class = "dendrogram"
		)
		if (attr(newcluster, "members") > 0 || !includeNodes) {
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
