## Visualization of effect of the alpha and gamma parameters

mnist_path <- "/Volumes/Datasets2/DATASETS/mnist/raw"
library(deepnet)
mnist <- deepnet::load.mnist(mnist_path)
data <- t(mnist$train$x)
labels <- mnist$train$y
visObj <- largeVis(data, sgd_batches = 1)

library(tidyr)
library(dplyr)
library(magrittr)
library(purrr)
inputs <- expand_grid(gamma = c(5,7,14), alpha = c(0.1, 1, 10))

agcoords <- map2(inputs$alpha, inputs$gamma, function(a, g) {
	obj <- projectKNNs(visObj$wij, dim = 2, alpha = a, gamma = g)
	data.frame(t(obj)) %>%
		set_colnames(c("x", "y")) %>%
		mutate(alpha = a, gamma = g)
}) %>% bind_rows()

agcoords$label = rep(labels, nrow(inputs))

library(ggplot2)

agcoords %>%
	ggplot(aes(x = x, y = y, color = factor(label))) +
	geom_point(size = 0.1, alpha = 0.1) +
	facet_grid(alpha ~ gamma, scales = "free") +
	theme_minimal()

Ks <- c(75, 150, 300)
Ms <- c(1, 5, 10)

coord_list <- c()

for (k in Ks) {
	visObj <- largeVis(data, K = k, sgd_batches = 1)
	for (m in Ms) {
		coords <- projectKNNs(visObj$wij, M = m)
		new_data <- data.frame(t(coords)) %>%
			set_colnames(c("x", "y")) %>%
			mutate(label = labels, M = m, K = k)
		coord_list <- c(coord_list, list(new_data))
	}
}

kmcoords <- bind_rows(coord_list)

kmcoords %>%
	ggplot(aes(x = x, y = y, color = factor(label))) +
	geom_point(size = 0.1, alpha = 0.1) +
	facet_grid(K ~ M) +
	theme_minimal() +
	guides(color = F)

pathToGraphFile <- "/Volumes/Datasets2/DATASETS/YouTubeCommunities/com-youtube.ungraph.txt"
pathToCommunities <- "/Volumes/Datasets2/DATASETS/YouTubeCommunities/com-youtube.top5000.cmty.txt"
youtube <- readr::read_tsv(pathToGraphFile, skip=4, col_names=FALSE)
youtube <- as.matrix(youtube)
youtube <- Matrix::sparseMatrix(i = youtube[, 1],
																j = youtube[, 2],
																x = rep(1, nrow(youtube)),
																dims = c(max(youtube), max(youtube)))
youtube <- youtube + t(youtube)
communities <- readr::read_lines(pathToCommunities)
communities <- lapply(communities,
											FUN = function(x) as.numeric(unlist(strsplit(x, "\t"))))
community_assignments <- rep(0,
														 nrow(youtube))
for (i in 1:length(communities)) community_assignments[communities[[i]]] <- i

wij <- buildWijMatrix(youtube)
youTube_coordinates <- projectKNNs(youtube, verbose=T)
youTube_coordinates <- data.frame(scale(t(youTube_coordinates)))
colnames(youTube_coordinates) <- c("x", "y")
youTube_coordinates$community <- factor(community_assignments)
youTube_coordinates$alpha <- factor(ifelse(youTube_coordinates$community == 0, 0.05, 0.2))
ggplot(youTube_coordinates, aes( x = x,
																 y = y,
																 color = community,
																 alpha = alpha,
																 size = alpha)) +
	geom_point() +
	scale_alpha_manual(values = c(0.005, 0.2), guide = FALSE) +
	scale_size_manual(values = c(0.03, 0.15), guide = FALSE) +
	scale_x_continuous("",
										 breaks = NULL, limits = c(-2.5,2.5)) +
	scale_y_continuous("",
										 breaks = NULL, limits = c(-2.5,2.5)) +
	ggtitle("YouTube Communities") +
	guides(color = FALSE) +
	theme_minimal()

