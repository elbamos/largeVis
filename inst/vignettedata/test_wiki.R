library(word2vec)
library(largeVis)
dataPath <- "/Volumes/Datasets2/DATASETS/Wiki_embedding/"

vectors <- t(read.wordvectors(paste0(dataPath, "wiki_word_vec_100d.txt"), type="bin"))
neighbors <- randomProjectionTreeSearch(vectors, distance_method = "Cosine", verbose = T, save_file = tempfile(pattern = "wikidocs"))
save(neighbors, file=paste0(dataPath, "wiki_word_neighbors.Rda"))

wij <- buildWijMatrix(neighbors$edgematrix)
coords <- projectKNNs(wij)
save(coords, file=paste0(dataPath, "wiki_word_coords.Rda"))

library(magrittr)

labels <- readr::read_tsv(paste0(dataPath, "wiki_word_label.txt"), col_names = FALSE) %>%
	set_colnames(c("word", "label"))

library(ggplot2)
library(dplyr)

data.frame(t(coords)) %>%
	set_colnames(c("x", "y")) %>%
	mutate(word = colnames(vectors)) %>%
	left_join(labels) %>%
	ggplot(aes(x = x, y = y, color = factor(label))) +
	geom_point(size = 0.1, alpha = 0.1) +
	theme_minimal() +
	guides(color = FALSE)
