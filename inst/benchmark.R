benchmark <- function(path,
                      samplepath,
                           K = 40,
                           tree_range = c(10, 20),
                           thresholds = c(10, 20, 30, 40, 50),
                           iters = c(3),
                           n = 10) {
  data <- readr::read_delim(path, delim = " ", col_names = F)
  data <- as.matrix(data)
  data <- scale(data)
  cat("Getting actual neighbors...\n")

  if (file.exists(samplepath)) {
    load(samplepath)
  } else {
    samples <- sample(nrow(data), n, replace = F)

    actualneighbors <- RANN::nn2(
      data, data[samples, ], k = K, treetype='kd',
    )$nn.idx - 1

    savedsamples <- list(samples = samples, neighbors = actualneighbors)
    save(savedsamples, file = samplepath)
    rm(samples)
    rm(actualneighbors)
  }
  data <- t(data)

  results <- data.frame(time = numeric(0),
                        precision = numeric(0),
                        n_trees = numeric(0),
                        max_iterations = numeric(0),
                        tree_threshold = numeric(0))

  for (n_trees in tree_range) {
	  for (max_iters in iters) {
		  for (threshold in thresholds) {
    		print(paste(n_trees, max_iters, threshold))
    		time <- system.time(
      		knns <- randomProjectionTreeSearch(data,
                                         K, n_trees, threshold, max_iters,
                                         verbose = TRUE)
    		)
    		precision <- lapply(1:n,
                        FUN = function(x)  sum(knns[, savedsamples$samples[x]] %in% savedsamples$neighbors[x, ]))

   			one_result <- data.frame(
                       time = time[1] + time[5],
                       precision = sum(as.numeric(precision)) / n,
                       n_trees = n_trees,
                       max_iterations = max_iters,
                       tree_threshold = threshold)
    		print(one_result)
    		readr::write_csv(one_result, path = "results.csv", append = TRUE)
  		}
  	}
  }
  return(results)
}

require( largeVis )
path <- "/mnt/hfsshare/DATASETS/sift/siftknns.txt"
samplepath <- "./samples.Rda"

results <- benchmark(path, samplepath, n = 10000, K = 100, tree_range = 10, thresholds = 20, iters = 2)
print(results)
