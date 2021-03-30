library(largeVis)
library(plyr)
library(dplyr)
library(magrittr)



get_true_nns <- function(data, k) {
#	bests <- RANN::nn2(data, k = k, eps = 0.0)
#	bests <- t(bests$nn.idx) - 1
 bests <- randomProjectionTreeSearch(data, K = k, n_trees = 200, max_iter = 5, verbose = T)
 bests$neighbors
}

calc_precision <- function(bests, tst) {
	precision <- lapply(1:ncol(bests),
											FUN = function(x)
												sum(bests[, x] %in% tst[, x]))

	sum(as.numeric(precision)) / (ncol(bests) * nrow(bests))
}


benchmark_one <- function(bests, data, k, trees, iters) {
	gc()
	start <- Sys.time()
	neighbors <- randomProjectionTreeSearch(data, K = k, n_trees = trees, max_iter = iters)
	end <- Sys.time()
	local_bests <- bests[1:k, ]
	score <- calc_precision(local_bests, neighbors$neighbors)
	time <- end - start
	units(time) <- "secs"
	data.frame(time = time, precision = score)
}

sift_path <- "/Volumes/Datasets2/DATASETS/sift/siftknns.txt"
sift <- readr::read_delim(sift_path, delim = " ", col_names = F)
sift <- t(as.matrix(sift))

bests <- get_true_nns(sift, 200)
gc()



test_range <- expand.grid(nns = c(50, 100, 150), n_trees = c(7, 14, 21, 42, 84), n_iters = c(0, 1, 2))

results <- adply(test_range, .margins = 1, .fun = function(x) {
	res <- benchmark_one(bests, sift, k = x$nns, trees = x$n_trees, iters = x$n_iters)
	cbind(res, x)
})


library(ggplot2)

results %>%
	mutate(speed = 1000000 / as.numeric(time)) %>%
	filter(n_iters < 5) %>%
	group_by(nns, n_trees) %>%
	arrange(n_iters) %>%
	ggplot(aes(x = speed, y = precision, color = factor(n_iters), group = factor(n_trees))) +
	geom_point(size = 0.6) +
  geom_line(size = 0.3, alpha = 0.3, color = "black") +
	facet_wrap(. ~ nns, scales = "free") +
	scale_y_sqrt("Accuracy") +
	scale_x_log10("Speed (points / second)") +
	theme_minimal()

save(results, file = "./inst/extdata/benchmark.Rda")
