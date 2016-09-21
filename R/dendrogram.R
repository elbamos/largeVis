#' @title as.dendrogram.hdbscan
#' @description Convert hdbscan object to dendrogram.
#' @export
#' @examples
#' data(iris)
#' neighbors <- randomProjectionTreeSearch(t(iris[,1:4]), K = 20)
#' edges <- buildEdgeMatrix(t(iris[,1:4]), neighbors)
#' hdbscanobj <- hdbscan(edges, minPts = 10, K = 5)
#' plot(as.dendrogram(hdbscanobj))
as.dendrogram.hdbscan <- function(x) {
	counts <- tabulate(x$hierarchy$nodemembership)
	C <- length(counts)
	counts <- counts + tabulate(x$hierarchy$parent, nbins = C)

	leafs <- lapply(1:length(x$hierarchy$nodemembership), FUN = function(i) structure(as.list(i),
			leaf = TRUE,
			members = 1L,
			label = i,
			cluster = x$clusters[i],
			probability = x$probabilities[i],
			height = x$hierarchy$lambda[i],
			class = "dendrogram"
		))

	clusters <- vector("list", C)
	h <- function(z) attr(z, "height")
	m <- function(z) attr(z, "members")
	for (i in C:1) {
		clusters[[i]] <- structure(c(leafs[which(x$hierarchy$nodemembership == i)],
																 clusters[which(x$hierarchy$parent == i & (1:C != i))]),
			members = sum(unlist(vapply(leafs[which(x$hierarchy$nodemembership == i)], FUN.VALUE = 0L, FUN = m))) +
				        sum(unlist(vapply(clusters[which(x$hierarchy$parent == i & (1:C != i))], FUN.VALUE = 0L, FUN = m))),
			leaf = FALSE,
			height = 1.1 * max(unlist(vapply(leafs[which(x$hierarchy$nodemembership == i)], FUN.VALUE = 0, FUN = h)),
												 unlist(vapply(clusters[which(x$hierarchy$parent == i & (1:C != i))], FUN.VALUE = 0, FUN = h))),
			selected = x$hierarchy$selected[i],
			cluster = i,
			stability = x$hierarchy$stability[i],
			label = paste("cluster", i, "stability",x$hierarchy$stability[i], ifelse(x$hierarchy$selected[i], "selected", "")),
			class = "dendrogram"
		)
		clusters[(x$hierarchy$parent == i) & (1:C != i)] <- NULL
	}
	clusters <- clusters[! sapply(clusters, is.null)]
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
