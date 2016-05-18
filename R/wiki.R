#' Voting data on wikipedia from inception until January, 2008.
#'
#' @format A symmetric sparse matrix in C-compressed format. Weights for present edges are either 1,
#' indicating that each node case a vote for the other, or 0.5.  Nodes with fewer than 5 votes were
#' removed from the dataset.
#'
#' @source \url{https://snap.stanford.edu/data/wiki-Vote.html}
"wiki"
