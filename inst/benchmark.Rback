benchmark <- function(path,
                           K = 40,
                           tree_range = c(10,20,50,100),
                           thresholds = c(10, 20, 50, 100),
                           iters = c(1,2,5),
                           n = 10000) {
  data <- readr::read_delim(path, delim = " ", col_names = F)
  data <- as.matrix(data)
  data <- scale(data)
  cat("Getting actual neighbors...\n")

  samples <- sample(nrow(data), n, replace = F)

  actualneighbors <- FNN::get.knnx(
    data, data[samples,], K, algorithm = 'brute'
  )$nn.index

  print(str(actualneighbors))
  data <- t(data)

  results <- data.frame(time = numeric(0), precision = numeric(0), n_trees = numeric(0),
                        max_iterations = numeric(0), tree_threshold = numeric(0))

  lapply(tree_range, FUN = function(n_trees) {
    lapply(iters, FUN = function(max_iters) {
      lapply(thresholds, FUN = function(threshold) {
        print(paste(n_trees, max_iters, threshold))
        time <- system.time(
          knns <- randomProjectionTreeSearch(data,
                                             K, n_trees, tree_threshold, max_iter,
                                             verbose = TRUE)
        )
        precision <- lapply(1:n, FUN = function(x)  sum(knns[,samples[x]] %in% actualneighbors[x,]))
        print(time)
        print(precision)
        results <- cbind(results,
                         data.frame(
                           time = time$user,
                           precision = precision,
                           n_trees = n_trees,
                           max_iterations = max_iters,
                           tree_threshold = threshold))
      })
    })
  })
  return(results)
}

require(largeVis)
path <- "/mnt/hfsshare/DATASETS/sift/siftknns.txt"

results <- benchmark(path)
print(results)
save(results, "benchmark.Rda")
