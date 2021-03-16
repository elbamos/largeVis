library(word2vec)
library(largeVis)
dataPath <- "/Volumes/Datasets2/DATASETS/Wiki_embedding/"

vectors <- t(read.wordvectors(paste0(dataPath, "wiki_doc_vec_100d.txt"), type="bin"))
neighbors <- randomProjectionTreeSearch(vectors, verbose = T)
save(neighbors, paste0(dataPath, "wiki_doc_neighbors.Rda"))