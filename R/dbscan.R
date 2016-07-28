#' OPTICS
#'
#' An implementation of the OPTICS algorithm. 
#' 
#' @param data Input data, where examples are columns.
#' @param neighbors An adjacency matrix of the type produced by \code{\link{randomProjectionTreeSearch}}
#' @param eps See \code{\link[dbscan]{optics}}.
#' @param minPts See \code{\link[dbscan]{optics}}.
#' @param eps_cl See \code{\link[dbscan]{optics}}.
#' @param xi See \code{\link[dbscan]{optics}}.
#'
#' @return An \code{\link[dbscan]{optics}} object.
#' 
#' @importFrom dbscan optics_cut opticsXi 
#' @export
optics <- function(data, 
                   neighbors, 
                   eps, 
                   minPts = nrow(data) + 1, 
                   eps_cl,
                   xi) {
  
  ret <- optics_int(as.matrix(data), neighbors, 
                    as.double(eps), as.integer(minPts))
  
  ret$minPts <- minPts
  ret$eps <- eps
  ret$eps_cl <- NA
  class(ret) <- "optics"
  
  ### find clusters
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
#' @param eps See \code{\link[dbscan]{optics}}.
#' @param minPts 
#'
#' @return An \code{\link[dbscan]{dbscan_fast}} object.
#' 
#' @export
dbscan <- function(data, 
                   neighbors, 
                   eps, 
                   minPts = nrow(data) + 1) {
  
  ret <- dbscan_int(as.matrix(data), 
                    as.matrix(neighbors), 
                    as.double(eps), 
                    as.integer(minPts))
  
  structure(list(cluster = ret, eps = eps, minPts = minPts),
            class = c("dbscan_fast", "dbscan"))
}
  