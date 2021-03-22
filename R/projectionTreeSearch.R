#' Find approximate k-Nearest Neighbors using random projection tree search.
#'
#' A fast and accurate algorithm for finding approximate k-nearest neighbors.
#'
#' Note that the algorithm does not guarantee that it will find K neighbors for each node. A
#' warning will be issued if it finds fewer neighbors than requested. If the input data contains
#' distinct partitionable clusters, try increasing the \code{tree_threshold} to increase the number
#' of returned neighbors.
#'
#' @param x A (potentially sparse) matrix, where examples are columnns and features are rows.
#' @param K How many nearest neighbors to seek for each node.
#' @param n_trees The number of trees to build.
#' @param max_iter Number of iterations in the neighborhood exploration phase.
#' @param distance_method One of "Euclidean" or "Cosine."
#' @param save_file If not null, the annoy index will be built on-disk at this path.
#' @param verbose Whether to print verbose logging using the \code{progress} package.
#'
#' @return A [K, N] matrix of the approximate K nearest neighbors for each vertex.
#' @export
randomProjectionTreeSearch <- function(x,
																			 K = 150,
																			 n_trees = 50,
																			 max_iter = 1,
																			 distance_method = "Euclidean",
																			 save_file = NULL,
																			 verbose = getOption("verbose", TRUE))
	UseMethod("randomProjectionTreeSearch")

#' @export
#' @rdname randomProjectionTreeSearch
randomProjectionTreeSearch.matrix <- function(x,
																							K = 150,
																							n_trees = 50,
																							max_iter = 1,
																							distance_method = "Euclidean",
																							save_file = NULL,
																							verbose = getOption("verbose", TRUE)) {
	if (verbose)
		cat("Searching for neighbors.\n")

	knns <- searchTrees(
		n_trees = as.integer(n_trees),
		K = as.integer(K),
		maxIter = as.integer(max_iter),
		data = x,
		distMethod = as.character(distance_method),
		saveFile = save_file,
		verbose = as.logical(verbose)
	)

	if (sum(colSums(knns$neighbors != -1) == 0) > 0)
		stop ("After neighbor search, no candidates for some nodes.")
	if (sum(is.na(knns$neighbors)) + sum(is.nan(knns$neighbors)) > 0)
		stop ("NAs or nans in neighbor graph.")
	if (verbose[1] && sum(knns$neighbors == -1) > 0)
		warning ("Wanted to find",
						 nrow(knns$neighbors) * ncol(knns$neighbors),
						 " neighbors, but only found",
						 ((nrow(knns$neighbors) * ncol(knns$neighbors)) - sum(knns$neighbors == -1)))

	attr(knns$edgematrix, "Metric") <- tolower(distance_method)

	return(knns)
}

#' @export
#' @rdname randomProjectionTreeSearch
randomProjectionTreeSearch.CsparseMatrix <- function(x,
																										 K = 150,
																										 n_trees = 50,
																										 max_iter = 1,
																										 distance_method = "Euclidean",
																										 save_file = NULL,
																										 verbose = getOption("verbose", TRUE)) {
	if (verbose)
		cat("Searching for neighbors.\n")

	knns <- searchTreesCSparse(
		n_trees = as.integer(n_trees),
		K = as.integer(K),
		maxIter = as.integer(max_iter),
		i = x@i,
		p = x@p,
		x = x@x,
		distMethod = as.character(distance_method),
		saveFile = save_file,
		verbose = as.logical(verbose)
	)

	if (sum(colSums(knns$neighbors != -1) == 0) > 0)
		stop ("After neighbor search, no candidates for some nodes.")
	if (sum(is.na(knns$neighbors)) + sum(is.nan(knns$neighbors)) > 0)
		stop ("NAs or nans in neighbor graph.")
	if (verbose[1] && sum(knns$neighbors == -1) > 0)
		warning ("Wanted to find",
						 nrow(knns$neighbors) * ncol(knns$neighbors),
						 " neighbors, but only found",
						 ((nrow(knns$neighbors) * ncol(knns$neighbors)) - sum(knns$neighbors == -1)))

	attr(knns$edgematrix, "Metric") <- tolower(distance_method)

	return(knns)
}

#' @export
#' @rdname randomProjectionTreeSearch
randomProjectionTreeSearch.TsparseMatrix <- function(x,
																										 K = 150,
																										 n_trees = 50,
																										 max_iter = 1,
																										 distance_method =
																										 	"Euclidean",
																										 save_file = NULL,
																										 verbose = getOption("verbose", TRUE)) {
	if (verbose)
		cat("Searching for neighbors.\n")

	knns <- searchTreesTSparse(
		n_trees = as.integer(n_trees),
		K = as.integer(K),
		maxIter = as.integer(max_iter),
		i = x@i,
		j = x@j,
		x = x@x,
		distMethod = as.character(distance_method),
		saveFile = save_file,
		verbose = as.logical(verbose)
	)

	if (sum(colSums(knns$neighbors != -1) == 0) > 0)
		stop ("After neighbor search, no candidates for some nodes.")
	if (sum(is.na(knns$neighbors)) + sum(is.nan(knns$neighbors)) > 0)
		stop ("NAs or nans in neighbor graph.")
	if (verbose[1] && sum(knns$neighbors == -1) > 0)
		warning ("Wanted to find",
						 nrow(knns$neighbors) * ncol(knns$neighbors),
						 " neighbors, but only found",
						 ((nrow(knns$neighbors) * ncol(knns$neighbors)) - sum(knns$neighbors == -1)))

	attr(knns$edgematrix, "Metric") <- tolower(distance_method)

	return(knns)
}


#' Run randomProjectionTreeSearch using a previously saved RccpAnnoy index
#' @param save_file The path to the index
#' @param D Dimensionality of the dataset
#' @param K The number of nearest neighbors to find
#' @param max_iter Number of iterations in the neighborhood exploration phase.
#' @param distance_method One of "Euclidean" or "Cosine" - must match the previously saved index.
#' @param verbose Whether to print verbose logging using the \code{progress} package.
#'
#' @return A [K, N] matrix of the approximate K nearest neighbors for each vertex.
#' @export
#' @rdname randomProjectionTreeSearch
randomProjectionTreeSearchFromIndex <- function(save_file,
																								 D,
																								 K = 150,
																								 max_iter = 1,
																								 distance_method = "Euclidean",
																								 verbose = getOption("verbose", TRUE)) {
	if (verbose)
		cat("Searching for neighbors.\n")


	knns <- searchTreesFromIndex(
		K = as.integer(K),
		D = as.integer(D),
		maxIter = as.integer(max_iter),
		distMethod = as.character(distance_method),
		saveFile = as.character(save_file),
		verbose = as.logical(verbose)
	)

	if (sum(colSums(knns$neighbors != -1) == 0) > 0)
		stop ("After neighbor search, no candidates for some nodes.")
	if (sum(is.na(knns$neighbors)) + sum(is.nan(knns$neighbors)) > 0)
		stop ("NAs or nans in neighbor graph.")
	if (verbose[1] && sum(knns$neighbors == -1) > 0)
		warning ("Wanted to find",
						 nrow(knns$neighbors) * ncol(knns$neighbors),
						 " neighbors, but only found",
						 ((nrow(knns$neighbors) * ncol(knns$neighbors)) - sum(knns$neighbors == -1)))

	attr(knns$edgematrix, "Metric") <- tolower(distance_method)

	return(knns)
}
