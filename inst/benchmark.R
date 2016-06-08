benchmark <- function(path,
                      samplepath,
                      K = 40,
                      tree_range = c(10, 20, 50),
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
      data, data[samples, ], k = K, treetype = "kd",
    )$nn.idx - 1

    savedsamples <- list(samples = samples, neighbors = actualneighbors)
    save(savedsamples, file = samplepath)
    rm(samples)
    rm(actualneighbors)
  }
  data <- t(data)

  for (n_trees in tree_range) {
    for (max_iters in iters) {
      for (threshold in thresholds) {
        print(paste(n_trees, max_iters, threshold))
        gc()
        start <- Sys.time()
          knns <- randomProjectionTreeSearch(data,
                                             K, n_trees, threshold, max_iters,
                                             verbose = TRUE)
        time <- Sys.time() - start
        units(time) <- "mins"
        precision <- lapply(1:n,
                            FUN = function(x)
          sum(knns[, savedsamples$samples[x]] %in% savedsamples$neighbors[x, ]))
        precision <- sum(as.numeric(precision)) / n

        one_result <- data.frame(
          time = time,
          precision = precision,
          n_trees = n_trees,
          max_iterations = max_iters,
          tree_threshold = threshold,
          method = "largeVis-emplace",
          tree_type = "",
          searchtype = "",
          eps = 0)
        print(one_result)
        readr::write_csv(one_result, path = "results.csv", append = TRUE)
        if (precision == K) break
      }
    }
  }
}

benchmarkRANN <- function(path,
                          samplepath,
                          K = 40,
                          tree_types = "kd", # c("kd", "bd"),
                          searchtypes = c("priority", "standard"),
                          epss = c(0.1, 0.2, 0.5),
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
      data, data[samples, ], k = K, treetype = "kd",
    )$nn.idx - 1

    savedsamples <- list(samples = samples, neighbors = actualneighbors)
    save(savedsamples, file = samplepath)
    rm(samples)
    rm(actualneighbors)
  }
  library(RANN)
  for (tree_type in tree_types) {
    for (searchtype in searchtypes) {
      for (eps in epss) {
        print(paste(tree_type, searchtype, eps))
        start <- Sys.time()
          knns <- nn2(data, k = K, query = data[savedsamples$samples, ],
                      treetype = tree_type,
                      searchtype = searchtype,
                      eps = eps)$nn.idx
        time <- Sys.time() - start
        units(time) <- "mins"
        precision <- lapply(1:n,
                            FUN = function(x)
          sum(knns[x, ] %in% (savedsamples$neighbors[x, ] + 1)))

        one_result <- data.frame(
          time = time,
          precision = sum(as.numeric(precision)) / n,
          n_trees = 0,
          max_iterations = 0,
          tree_threshold = 0,
          method = "RANN",
          tree_type = tree_type,
          searchtype = searchtype,
          eps = eps)
        print(one_result)
        readr::write_csv(one_result, path = "results.csv", append = TRUE)
      }
    }
  }
}

benchmarkAnnoy <- function(path,
                           samplepath,
                           K = 40,
                           tree_range = c(10, 20, 50, 100),
                           n = 10,
                           full = FALSE) {
  data <- readr::read_delim(path, delim = " ", col_names = F)
  data <- as.matrix(data)
  data <- scale(data)
  cat("Getting actual neighbors...\n")

  if (file.exists(samplepath)) {
    load(samplepath)
  } else {
    samples <- sample(nrow(data), n, replace = F)

    actualneighbors <- RANN::nn2(
      data, data[samples, ], k = K, treetype = "kd",
    )$nn.idx - 1

    savedsamples <- list(samples = samples, neighbors = actualneighbors)
    save(savedsamples, file = samplepath)
    rm(samples)
    rm(actualneighbors)
  }
  library(RcppAnnoy)
  for (n_trees in tree_range) {
    knns <- list()
    print(n_trees)
    gc()
    start <- Sys.time()
    a <- new(AnnoyEuclidean, ncol(data))
    for (i in 1:nrow(data)) a$addItem(i - 1, data[i, ])
    a$build(n_trees)
    if (! full) for (i in 1:n)
      knns[[i]] <- a$getNNsByItem(item = savedsamples$samples[i] - 1,
                                               size = K)
    else for (i in 1:nrow(data))
      knns[[i]] <- a$getNNsByItem(item = i - 1, size = K)
    time <- Sys.time() - start
    units(time) <- "mins"

    if (full) knns <- knns[savedsamples$samples]

    precision <- lapply(1:n,
                        FUN = function(x)
                          sum(knns[[x]] %in% savedsamples$neighbors[x, ]))

    if (full) method <- "RcppAnnoy-Full"
    else method <- "RcppAnnoy"

    one_result <- data.frame(
      time = time,
      precision = sum(as.numeric(precision)) / n,
      n_trees = n_trees,
      max_iterations = 0,
      tree_threshold = 0,
      method = method,
      tree_type = "",
      searchtype = "",
      eps = 0)
    print(one_result)
    readr::write_csv(one_result, path = "results.csv", append = TRUE)
  }
}

require( largeVis )
path <- "/mnt/hfsshare/DATASETS/sift/siftknns.txt"
samplepath <- "./samples.Rda"

# Annoyresults <- benchmarkAnnoy(path,
#                                samplepath,
#                                tree_range = c(10, 20, 50, 100, 200, 400),
#                                n = 10000,
#                                K = 100,
#                                full = FALSE)
# results <- benchmark(path,
#                      samplepath,
#                      n = 10000,
#                      tree_range = c(10, 20, 50, 100, 200),
#                      thresholds = c(128),
#                      iters = c(1, 0, 2, 3),
#                      K = 100)
results2 <- benchmark(path,
                     samplepath,
                     n = 10000,
                     tree_range = c(10, 20, 50),
                     thresholds = c(10, 20, 50, 80, 256, 512),
                     iters = c(1, 0, 2, 3),
                     K = 100)
Annoyresults2 <- benchmarkAnnoy(path,
                               samplepath,
                               tree_range = c(10, 20, 50, 100, 200, 400),
                               n = 10000,
                               K = 100,
                               full = TRUE)
# RANNresults <- benchmarkRANN(path,
#                              samplepath,
#                              epss = c(.1, .5, 1,2,5),
#                              n = 10000,
#                              K = 100)

print(results)
