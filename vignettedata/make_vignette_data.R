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

save(agcoords, file = "./inst/vignettedata/agcoords.Rda")

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

save(kmcoords, file = "./inst/vignettedata/kmcoords.Rda")
