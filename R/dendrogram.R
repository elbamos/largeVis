#' @title as.dendrogram.hdbscan
#' @description Convert hdbscan object to dendrogram
#' @export
as.dendrogram.hdbscan <- function(x) {
	leaves <- list()
	for (i in 1:length(x$hierarchy$nodemembership)) {
		cluster <- as.character(x$hierarchy$nodemembership[i])
		if (!is.na(cluster)) {
			newleaf <- list(i)
			attr(newleaf, "members") <- 1
			attr(newleaf, "leaf") <- TRUE
			attr(newleaf, "height") <- x$hierarchy$lambda[i]
			class(newleaf) <- "dendrogram"
			if (cluster %in% names(leaves)) {
				leaves[[cluster]][[length(leaves[[cluster]]) + 1]] <- newleaf
			} else {
				newcluster <- list()
				newcluster[[1]] <- newleaf
				leaves[[cluster]] <- newcluster
			}
		}
	}
	roots <- list()
	for (i in 1:length(x$hierarchy$parent)) {
		cluster <- as.character(i)
		members <- 0
		for (j in leaves[[cluster]]) {
			members <- members + attr(j, "members")
		}
		if (! cluster %in% names(leaves)) leaves[[cluster]] <- list()
		attr(leaves[[cluster]], "members") <- members
		attr(leaves[[cluster]], "height") <- x$hierarchy$stability[i]
		attr(leaves[[cluster]], "selected") <- x$hierarchy$selected[i]
		attr(leaves[[cluster]], "leaf") <- FALSE
		class(leaves[[cluster]]) <- "dendrogram"
		print(paste(i, x$hierarchy$parent[i]))
		if (x$hierarchy$parent[i] + 1 == i) roots[[length(roots) + 1]] <- leaves[[cluster]]
		else {
			parent <- as.character(x$hierarchy$parent[i] + 1)
			leaves[[parent]][[length(leaves[[parent]]) + 1]] <- leaves[[cluster]]
			leaves[[cluster]] <- NULL
		}
	}
	if (length(roots) > 1) {
		class(roots) <- "dendrogram"
		members <- 0
		for (i in roots) members <- members + attr(i, "members")
		attr(roots, "members") <- members
		attr(roots, "height") <- 100
		return(roots)
	} else {
		return(roots[[names(roots)[1]]])}
}
