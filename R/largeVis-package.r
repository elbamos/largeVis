#' largeVis: high-quality visualizations for large, high-dimensionality datasets
#'
#' This is an implementation of the \code{largeVis} algorithm by Tang et al., and related functions and algorithms.
#'
#' \code{largeVis} estimates a low-dimensional embedding for high-dimensional data, where the distance between vertices
#' in the low-dimensional space is proportional to the distance between them in the high-dimensional space. The algorithm
#' works in 4 phases:
#'
#' \itemize{
#' \item  Estimate candidate nearest-neighbors for each vertex by building \code{n.trees} random projection trees.
#' \item  Estimate \code{K} nearest-neighbors for each vertex by visiting each vertex' 2d-degree neighbors (its neighbors' neighbors).
#' This is repeated \code{max.iter} times.  Note that the original paper suggested a \code{max.iter} of 1, however a larger number
#' may be appropriate for some datasets if the algorithm has trouble finding K neighbors for every vertex.
#' \item Estimate \eqn{p_{j|i}}, the conditional probability that each edge found in the previous step is actually to a
#' nearest neighbor of each of its nodes.
#' \item Using stochastic gradient descent, estimate an embedding for each vertex in the low-dimensional space.
#' }
#'
#' The nearest-neighbor search functionality is also available as a separate function, where it offers an extremely fast approximate
#' nearest-neighbor search.  (See the Benchmarks vignette for details.)
#'
#' The package also includes implementations of the HDBSCAN, DBSCAN, and OPTICS clustering algorithms, and LOF outlier detection, optimized to use
#' data generated by running \code{largeVis}.
#'
#' @references Jian Tang, Jingzhou Liu, Ming Zhang, Qiaozhu Mei. \href{https://arxiv.org/abs/1602.00370}{Visualizing Large-scale and High-dimensional Data.}
#' R. Campello, D. Moulavi, and J. Sander, Density-Based Clustering Based on Hierarchical Density Estimates In: Advances in Knowledge Discovery and Data Mining, Springer, pp 160-172. 2013
#' Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Jorg Sander (1999). OPTICS: Ordering Points To Identify the Clustering Structure. ACM SIGMOD international conference on Management of data. ACM Press. pp. 49-60.
#' Martin Ester, Hans-Peter Kriegel, Jorg Sander, Xiaowei Xu (1996). Evangelos Simoudis, Jiawei Han, Usama M. Fayyad, eds. A density-based algorithm for discovering clusters in large spatial databases with noise. Proceedings of the Second International Conference on Knowledge Discovery and Data Mining (KDD-96). AAAI Press. pp. 226-231. ISBN 1-57735-004-9.
#' @name largeVis-package
#' @docType package
#' @useDynLib largeVis
#' @aliases NULL
#' @import Rcpp
#' @importFrom RcppParallel RcppParallelLibs
"_PACKAGE"


#' @importFrom devtools r_env_vars
#' @importFrom RcppParallel setThreadOptions
checkCRAN <- function() {
	if (devtools::r_env_vars()[['NOT_CRAN']] != "true") RcppParallel::setThreadOptions(numThreads = 2)
}
