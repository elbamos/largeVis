hyperplane <- function(x1, x2, tomatch) {
  v <- (x2 - x1)
  m <- (x1 + x2) / 2
  return((tomatch %*% v) - (m %*% v)[1])
}

partition <- function(indices, .threshold, .data) {
  if (length(indices) == 0) return(NULL)
  if (length(indices) <= .threshold) return(list(indices))
  xs <- sample(indices,2)
  x1 <- .data[xs[1],]
  x2 <- .data[xs[2],]
  selections <- hyperplane(x1, x2, .data[indices,]) >= 0
  c(partition(indices[selections], .threshold, .data), partition(indices[! selections], .threshold, .data))
}
