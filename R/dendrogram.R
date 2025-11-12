#' @title as.dendrogram.hdbscan
#' @description Convert an hdbscan object into a dendrogram compatible with the \code{stats} package.
#'
#' @details The hdbscan algorithm works by first building a hierarchy based on a minimal spanning tree,
#' then consolidating nodes according to rules in the algorithm. The algorithm then selects some of the
#' consolidated nodes as clusters, deselecting others. For example, if Node A has children B and C, the
#' algorithm might select A, and then all points under A, B, and C would be assigned to the same cluster.
#' Or, it might deselect A, and select B and C instead. In that case, the nodes under B would be assigned
#' to one cluster, the nodes under C to a different cluster, and nodes under A but not B or C would not be
#' assigned to any cluster.
#'
#' This function returns a dendrogram of the consolidated cluster hierarchy. Whether a node was selected
#' as a cluster is stored as an attribute of each dendrogram node.
#'
#' @param object An \code{hdbscan} object.
#' @param includeNodes Whether individual nodes should be included in the dendrogram as leaves.
#'   If \code{FALSE} (default), only the cluster hierarchy is shown, with "fallen points" (points that
#'   belong directly to a cluster without being in a sub-cluster) aggregated. Can cause a substantial
#'   increase in the size of the object if \code{TRUE}.
#' @param ... For compatibility with \code{as.dendrogram}, currently ignored.
#'
#' @return A \code{dendrogram} object, where nodes have the following attributes:
#' \describe{
#'   \item{'leaf'}{As in \code{\link[stats]{dendrogram}}.}
#'   \item{'members'}{As in \code{\link[stats]{dendrogram}}. The total number of observations (points) under this node.}
#'   \item{'height'}{The distance at which this cluster was formed (1/lambda_birth for clusters, 1/lambda for points).}
#'   \item{'label'}{A descriptive label for the node.}
#'   \item{'midpoint'}{As in \code{\link[stats]{dendrogram}}, used for plotting.}
#'   \item{'cluster'}{The cluster ID number. Missing for individual points when includeNodes=TRUE.}
#'   \item{'stability'}{The cluster's stability score. Missing for leaves.}
#'   \item{'selected'}{Boolean indicating whether the cluster was selected in the final clustering. Missing for leaves.}
#'   \item{'lambda_birth'}{The lambda value at which the cluster was born. Missing for leaves.}
#'   \item{'lambda_death'}{The lambda value at which the cluster died/split. Missing for leaves.}
#'   \item{'GLOSH'}{The GLOSH outlier score (only for individual point leaves when includeNodes=TRUE).}
#'   \item{'probability'}{The membership probability (only for individual point leaves when includeNodes=TRUE).}
#' }
#'
#' @importFrom stats as.dendrogram
#' @export
#' @examples
#' \dontrun{
#' data(iris)
#' dat <- as.matrix(iris[, 1:4])
#' vis <- largeVis(t(dat), K = 20, sgd_batches = 1, verbose = FALSE)
#' hdbscanobj <- hdbscan(vis, minPts = 10, K = 5, verbose = FALSE)
#'
#' # Create dendrogram showing only cluster hierarchy
#' dend <- as.dendrogram(hdbscanobj)
#' plot(dend)
#'
#' # Create dendrogram including individual points as leaves
#' dend_full <- as.dendrogram(hdbscanobj, includeNodes = TRUE)
#' plot(dend_full)
#' }
as.dendrogram.hdbscan <- function(object, includeNodes = FALSE, ...) {

  hierarchy <- object$hierarchy
  n_clusters <- length(hierarchy$parent)
  n_points <- length(hierarchy$nodemembership)

  # Build dendrogram recursively for each cluster
  build_cluster_dend <- function(cluster_id) {

    # Find child clusters (clusters whose parent is this cluster)
    child_cluster_ids <- which(hierarchy$parent == cluster_id & seq_along(hierarchy$parent) != cluster_id)

    # Find points that belong directly to this cluster
    point_ids <- which(hierarchy$nodemembership == cluster_id)

    # Build child dendrograms
    children <- list()

    # Add child clusters
    if (length(child_cluster_ids) > 0) {
      for (child_id in child_cluster_ids) {
        children <- c(children, list(build_cluster_dend(child_id)))
      }
    }

    # Add individual points if requested
    if (includeNodes && length(point_ids) > 0) {
      for (pt_idx in point_ids) {
        point_dend <- structure(
          pt_idx,
          leaf = TRUE,
          members = 1L,
          label = pt_idx,
          height = 1 / hierarchy$lambda[pt_idx],
          GLOSH = object$GLOSH[pt_idx],
          probability = object$probabilities[pt_idx],
          cluster_membership = object$clusters[pt_idx],
          class = "dendrogram"
        )
        children <- c(children, list(point_dend))
      }
    }

    # If not including nodes, create a single aggregate leaf for fallen points
    if (!includeNodes && length(point_ids) > 0) {
      n_fallen <- length(point_ids)
      max_lambda <- max(hierarchy$lambda[point_ids], na.rm = TRUE)
      fallen_dend <- structure(
        n_fallen,
        leaf = TRUE,
        members = n_fallen,
        label = paste(n_fallen, "point(s)"),
        height = 1 / max_lambda,
        class = "dendrogram"
      )
      children <- c(children, list(fallen_dend))
    }

    # Calculate total members
    total_members <- if (length(children) > 0) {
      sum(vapply(children, function(x) attr(x, "members"), 0L))
    } else {
      0L
    }

    # Sort children by height for better dendrogram structure
    if (length(children) > 1) {
      heights <- vapply(children, function(x) attr(x, "height"), 0)
      children <- children[order(heights)]
    }

    # Calculate midpoint for plotting
    midpoint <- if (length(children) > 1) {
      # Standard dendrogram midpoint calculation
      (attr(children[[1]], "members") +
         sum(vapply(children[-length(children)], function(x) attr(x, "members"), 0L))) / 2
    } else if (length(children) == 1) {
      0
    } else {
      0.5
    }

    # Determine height for this cluster
    # Use lambda_birth if available, otherwise use max child height
    cluster_height <- if (!is.na(hierarchy$lambda_birth[cluster_id]) &&
                          hierarchy$lambda_birth[cluster_id] > 0) {
      1 / hierarchy$lambda_birth[cluster_id]
    } else if (length(children) > 0) {
      max(vapply(children, function(x) attr(x, "height"), 0)) * 1.1
    } else {
      0
    }

    # Build the cluster label
    label <- paste("Cluster", cluster_id)
    if (!is.na(hierarchy$stability[cluster_id])) {
      label <- paste0(label, " (stab:", signif(hierarchy$stability[cluster_id], 3), ")")
    }
    if (!is.na(hierarchy$selected[cluster_id]) && hierarchy$selected[cluster_id]) {
      label <- paste0(label, " *")
    }

    # Create the dendrogram node
    cluster_dend <- structure(
      children,
      members = total_members,
      leaf = (total_members == 0),
      height = cluster_height,
      midpoint = midpoint,
      label = label,
      cluster = cluster_id,
      stability = hierarchy$stability[cluster_id],
      selected = hierarchy$selected[cluster_id],
      lambda_birth = hierarchy$lambda_birth[cluster_id],
      lambda_death = hierarchy$lambda_death[cluster_id],
      class = "dendrogram"
    )

    return(cluster_dend)
  }

  # Find root clusters (clusters with no parent or whose parent is themselves)
  root_ids <- which(hierarchy$parent == seq_along(hierarchy$parent) |
                    is.na(hierarchy$parent))

  # If there are no roots identified, cluster 1 is usually the root
  if (length(root_ids) == 0) {
    root_ids <- 1
  }

  # Build dendrograms for all root clusters
  root_dends <- lapply(root_ids, build_cluster_dend)

  # If there's only one root, return it directly
  if (length(root_dends) == 1) {
    return(root_dends[[1]])
  }

  # If there are multiple roots, create a super-root
  # Sort by height
  heights <- vapply(root_dends, function(x) attr(x, "height"), 0)
  root_dends <- root_dends[order(heights)]

  total_members <- sum(vapply(root_dends, function(x) attr(x, "members"), 0L))
  max_height <- max(heights)

  midpoint <- (attr(root_dends[[1]], "members") +
                 sum(vapply(root_dends[-length(root_dends)], function(x) attr(x, "members"), 0L))) / 2

  super_root <- structure(
    root_dends,
    members = total_members,
    leaf = FALSE,
    height = max_height * 1.2,
    midpoint = midpoint,
    label = "Root",
    class = "dendrogram"
  )

  return(super_root)
}
