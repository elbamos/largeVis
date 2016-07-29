#' OPTICS
#'
#' An implementation of the OPTICS algorithm. 
#' 
#' @param data Input data, where examples are columns.
#' @param neighbors An adjacency matrix of the type produced by \code{\link{randomProjectionTreeSearch}}
#' @param edges A weighted graph of the type produced by \code{\link{buildEdgeMatrix}}.
#' @param eps See \code{\link[dbscan]{optics}}.
#' @param minPts See \code{\link[dbscan]{optics}}.
#' @param eps_cl See \code{\link[dbscan]{optics}}.
#' @param xi See \code{\link[dbscan]{optics}}.
#' 
#' @details This is a preliminary implementation of the OPTICS algorithm that attempts
#' to leverage the \code{largeVis} nearest-neighbor search. 
#' 
#' One of \code{neighbors} or \code{edges} must be specified. If \code{edges} is missing,
#' \code{data} must also be given. If \code{data} is given along with either \code{edges}
#' or \code{neighbors}, the algorithm will attempt a more thorough search.
#' 
#' @note Support for dbscan and optics are preliminary, and not fully tested for
#' correctness.
#'
#' @return An \code{\link[dbscan]{optics}} object.
#' 
#' @importFrom dbscan optics_cut opticsXi 
#' @export
optics <- function(data = NULL, 
                   neighbors = NULL, 
                   edges = NULL,
                   eps, 
                   minPts = nrow(data) + 1, 
                   eps_cl,
                   xi,
                   verbose = TRUE) {
  if (! is.null(edges) && ! is.null(data))
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
  
  if(!missing(eps_cl)) ret <-optics_cut(ret, eps_cl)
  if(!missing(xi)) ret <- opticsXi(ret, xi)
  
  ret
}

#' dbscan
#' 
#' An implementation of the dbscan algorithm.
#'
#' @param data Input data, where examples are columns.
#' @param neighbors An adjacency matrix of the type produced by \code{\link{randomProjectionTreeSearch}}
#' @param edges A weighted graph of the type produced by \code{\link{buildEdgeMatrix}}.
#' @param eps See \code{\link[dbscan]{dbscan}}.
#' @param minPts Minimum size of a cluster.
#' 
#' @details This is a preliminary implementation of the OPTICS algorithm that attempts
#' to leverage the \code{largeVis} nearest-neighbor search. 
#' 
#' One of \code{neighbors} or \code{edges} must be specified. If \code{edges} is missing,
#' \code{data} must also be given. If \code{data} is given along with either \code{edges}
#' or \code{neighbors}, the algorithm will attempt a more thorough search.
#' 
#' @note Support for dbscan and optics are preliminary, and not fully tested for
#' correctness.
#' 
#' @return An \code{\link[dbscan]{dbscan_fast}} object.
#' 
#' @export
dbscan <- function(data = NULL, 
                   neighbors = NULL,
                   edges = NULL,
                   eps, 
                   minPts = nrow(data) + 1,
                   verbose = TRUE) {
  
  if (! is.null(edges) && ! is.null(data))
    ret <- dbscan_e(edges = edges,  
                    eps = as.double(eps), minPts = as.integer(minPts),
                    verbose = verbose)
  else if (! is.null(edges)) 
    ret <- dbscan_ed(edges = edges, data = data, 
                     eps = as.double(eps), minPts = as.integer(minPts),
                     verbose = verbose)
  else 
    ret <- dbscan_nd(neighbors = neighbors, data = data, 
                     eps = as.double(eps), minPts = as.integer(minPts),
                     verbose = verbose)
  
  structure(list(cluster = ret, eps = eps, minPts = minPts),
            class = c("dbscan_fast", "dbscan"))
}
  