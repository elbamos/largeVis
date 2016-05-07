
# returns a vector positioning each vertex in tomatch above or below a hyperplane equidistant between x1 and x2
#hyperplane <- function(x1, x2, tomatch) rowSums(tomatch * 2 * (x1 - x2)) + sum(x2^2) - sum(x1^2)
hyperplane <- function(x1, x2, tomatch) {
  v <- (x2 - x1)
  m <- (x1 + x2) / 2
  return((tomatch %*% v) - (m %*% v)[1])
  # return(cbind(
  #   (tomatch %*% v) - (m %*% v)[1],
  #   t(v %*% t(tomatch - m)),
  #   rowSums(tomatch * 2 * (x1 - x2)) + sum(x2^2) - sum(x1^2)
  # ))
}

# function:
#   randomly sample two points from node
#   use those two points to define hyperplane equidistant between them
#   divide points into the two spaces
#   for each space -- if number of points > tree.threshold, recurse
partition <- function(indices, .threshold, .data) {
  if (length(indices) == 0) return(NULL)
  if (length(indices) <= .threshold) return(list(indices))
  xs <- sample(indices,2)
  x1 <- .data[xs[1],]
  x2 <- .data[xs[2],]
  selections <- hyperplane(x1, x2, .data[indices,]) >= 0
  # recurse
  c(partition(indices[selections], .threshold, .data), partition(indices[! selections], .threshold, .data))
}
