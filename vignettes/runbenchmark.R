require(largeVis)
tree_range <- c(10, 20, 50, 100)
iters <- c(1, 2, 5, 10)
thresholds <- c(10, 20, 50, 100)
path <- "/mnt/hfsshare/DATASETS/sift/siftknns.txt"
K <- 40
n <- 10000

source("../inst/benchmark.R")

results <- data.frame(time = numeric(0), precision = numeric(0), n_trees = numeric(0), 
                      max_iterations = numeric(0), tree_threshold = numeric(0))

lapply(tree_range, FUN = function(n_trees) {
  lapply(iters, FUN = function(max_iters) {
    lapply(thresholds, FUN = function(threshold) {
      print(paste(n_trees, max_iters, threshold))
      one_result <- benchmark(path, 
                              K = K, 
                              n = n, 
                              n_trees = n_trees, 
                              max_iter = max_iter,
                              tree_threshold = tree_threshold)
      results <- cbind(results, 
                       data.frame(
                         time = one_result$time$user, 
                         precision = one_result$precision, 
                         n_trees = n_trees, 
                         max_iterations = max_iters, 
                         tree_threshold - threshold))
    })
  })
})
