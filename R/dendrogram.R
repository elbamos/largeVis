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
	C <- length(x$hierarchy$parent)
	counts <- tabulate(x$hierarchy$nodemembership, nbins = C) + tabulate(x$hierarchy$parent, nbins = C)

	leafs <- lapply(1:length(x$hierarchy$nodemembership), FUN = function(i)
		structure(as.list(i),
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
		children <- which(x$hierarchy$nodemembership == i)
		children <- leafs[children]
		cousins <- which(x$hierarchy$parent == i & (1:C != i))
		cousins <- clusters[cousins]
		cousins <- cousins[!sapply(cousins, is.null)]
		newcluster <- structure(
			.Data = c(children, cousins),
			members = sum(unlist(vapply(children, FUN.VALUE = 0L, FUN = m))) +
				        sum(unlist(vapply(cousins, FUN.VALUE = 0L, FUN = m))),
			leaf = FALSE,
			height = 1.1 * max(0, unlist(vapply(children, FUN.VALUE = 0, FUN = h)),
												    unlist(vapply(cousins, FUN.VALUE = 0, FUN = h))),
			selected = x$hierarchy$selected[i],
			cluster = i,
			stability = x$hierarchy$stability[i],
			label = paste("cluster", i, "stability",x$hierarchy$stability[i], ifelse(x$hierarchy$selected[i], "selected", "")),
			class = "dendrogram"
		)
		if (attr(newcluster, "members") > 0) {
			clusters[[i]] <- newcluster
		}
		clusters[(x$hierarchy$parent == i) & (1:C != i)] <- NA
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
