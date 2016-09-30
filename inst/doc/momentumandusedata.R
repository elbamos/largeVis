## ----setupvignette,eval=T,echo=F,warning=F,error=F,message=F-------------
require(ggplot2, 
        quietly = TRUE)
require(RColorBrewer, 
        quietly = TRUE)
require(wesanderson, 
        quietly = TRUE)
library(dplyr)
knitr::opts_chunk$set(collapse = TRUE, 
                      comment = "#>",
                      cache = FALSE, 
											echo = FALSE)
colors_discrete <- function(x) 
  rep(wes_palette("Darjeeling", n = min(x, 5)), 2)[1:x]

nacol <- colors_discrete(4)[4]
theme_set(
  theme_bw() %+replace%
  theme(
    legend.title = element_text(size = rel(0.8),
                              face = "bold"),
    legend.margin = unit(0, "cm"),
    legend.position = "right",
    legend.key.size = unit(2, "lines"),
    legend.text = element_text(size = unit(8, "points")), 
    axis.title.y = element_text(angle = 90),
    axis.text = element_text(size = rel(0.7)),
    plot.margin = unit(c(0, 0.5, 1, 0), "lines"), 
    axis.title = element_text(size = rel(0.8),
                              face = "bold"),
    title = element_text(size = rel(0.9))
  ) 
)
rebuild <- FALSE
par(mar = c(0.3, 0.3, 0.3, 0.3))
require(largeVis,quietly = TRUE)

## ----prepmnist-----------------------------------------------------------
if (rebuild) {
	load("../../largeVisData/mnist/test.RData")
	load("../../largeVisData/mnist/train.RData")
	mnist <- t(rbind(trainData, testData)) - 0.5
	preVis <- largeVis(mnist, K = 100, tree_threshold = 200, n_trees = 100, sgd_batches = 100, verbose = TRUE)
	set.seed(1974)
	labels <- rbind(trainLabels, testLabels)
	labels <- apply(labels, 1, function(x) which(x == 1))
	labels <- labels - 1
}

## ----usedegree-----------------------------------------------------------
if (rebuild) {
	noDegree <- projectKNNs(preVis$wij, seed = 1974, verbose = TRUE)
	degree <- projectKNNs(preVis$wij, useDegree = TRUE, verbose = TRUE)
	degreeCoords <- data.frame(rbind(t(noDegree), t(degree)), 
														 label = factor(c(labels, labels)), 
														 degree = factor(rep(c("weights", "degree"), 
														 										each = ncol(degree))))
	save(degreeCoords, file = paste(sep = "/", 
																	system.file(package = "largeVis", "extdata"), 
																	"degreeCoords.Rda"))
}

## ----drawdegree,fig.width=5,fig.height=4,fig.align='center',fig.cap='Effect of useDegree'----
load(system.file(package = "largeVis", "extdata/degreeCoords.Rda"))
degreeCoords %>%
	ggplot(aes(x = X1, y = X2, color = label)) + 
	geom_point(size = 0.05, alpha = 0.3) + 
	facet_grid(. ~ degree) +
	scale_x_continuous("", labels = NULL, breaks = NULL) + 
	scale_y_continuous("", labels = NULL, breaks = NULL) +
	guides(color = FALSE) +
	ggtitle("Effect of useDegree")

## ----momentum------------------------------------------------------------
if (rebuild) {
	starterCoords <- matrix(runif(n = 2 * ncol(mnist)) - 0.5, ncol = 2)
	timex <- system.time(firstCoords <- data.frame(t(
            projectKNNs(preVis$wij, coords = t(starterCoords), sgd_batches = 0.1, verbose = TRUE)), lambda = 0, batches = 0.1, label = labels))
	timex <- data.frame(elapsed = timex["elapsed"], batches = 1)
	for (batches in c(0.3, 0.8)) {
		newtime <- system.time(newCoords <- data.frame(t(projectKNNs(preVis$wij, verbose = TRUE, sgd_batches = batches, coords = t(starterCoords)))))
		print(str(newCoords))
		timex <- rbind(timex, data.frame(elapsed = newtime["elapsed"], batches = batches))
		newCoords$lambda <- 0
		newCoords$batches <- batches
		newCoords$label <- labels
		firstCoords <<- rbind(firstCoords, newCoords)
	}
	for (lambda in c(0.4, 0.9)) {
			for (batches in c(0.1, 0.3, 0.5)) {
				newtime <- system.time(newCoords <- data.frame(t(projectKNNs(preVis$wij, verbose = TRUE, sgd_batches = batches, momentum = lambda, coords = t(starterCoords)))))
				timex <<- rbind(timex, data.frame(elapsed = newtime["elapsed"], batches = batches))
				newCoords$lambda <- lambda
				newCoords$batches <- batches
				newCoords$label <- labels
				firstCoords <<- rbind(firstCoords, newCoords)
			}
	}
	momentumCoords <- firstCoords
	momentumCoords$label <- factor(momentumCoords$label)
	save(momentumCoords, file = paste(sep = "/", system.file(package = "largeVis", "extdata"), "momentumCoords.Rda"))
	save(timex, file = paste(sep = "/", system.file(package = "largeVis", "extdata"), "momentumTimes.Rda"))
}

## ----drawmomentum,warning=FALSE,fig.width=6,fig.height=5,fig.align='center',fig.cap='Effect of Momentum'----
load(system.file(package = "largeVis", "extdata/momentumCoords.Rda"))
momentumCoords %>%
	ggplot(aes(x = X1, y = X2, color = label)) + 
	geom_point(size = 0.05, alpha = 0.3) + 
	facet_grid(batches ~ lambda, scales = "free", labeller = label_bquote(cols = lambda == .(lambda), 
																																				rows = b == .(batches))) +
	scale_x_continuous("", limits = c(-40, 40), labels = NULL, breaks = NULL) + 
	scale_y_continuous("", limits = c(-40, 40), labels = NULL, breaks = NULL) +
	guides(color = FALSE) +
	ggtitle("Effect of Momentum and Reduced Training Batches")

## ----dbscan--------------------------------------------------------------
library(clusteringdatasets)
data(spiral)
dat <- spiral[, 1:2]
neighbors <- randomProjectionTreeSearch(t(dat), K = 20)
edges <- buildEdgeMatrix(t(dat), neighbors = neighbors)
par(mfrow = c(3, 3), mar = c(0.1, 0.1, 0.1, 0.1))
set <- rbind(Map(f = function(y) {
	rbind(Map(f = function(x) {
		clust = lv_dbscan(edges = edges, neighbors = neighbors, eps = x, minPts = y)$cluster
		data.frame(cluster = clust, eps = x, minPts = y)
	}, c(1, 3, 5)))
}, c(5, 10, 20)))
set <- lapply(set, FUN = bind_rows)
set <- bind_rows(set)
set$x <- rep(dat$x, 9)
set$y <- rep(dat$y, 9)
set$cluster <- factor(set$cluster)
set$eps <- factor(paste("epsilon == ", set$eps, sep = ""))
set$minPts <- factor(paste("minPts == ", set$minPts, sep = ""))
ggplot(data = set, aes(x = x, y = y, color = cluster)) +
	geom_point(size = 0.5, alpha = 0.7) +
	facet_grid(minPts ~ eps, labeller = label_parsed) + 
	scale_x_continuous("", breaks = NULL) +
	scale_y_continuous("", breaks = NULL) +
	guides(color = FALSE) +
	ggtitle("Effect of eps and minPts on DBSCAN results")

## ----optics,fig.width=5--------------------------------------------------
library(dbscan, quietly = TRUE)
optClust <- lv_optics(edges = edges, neighbors = neighbors, eps = 5, minPts = 5)
par(mfrow = c(1, 1))
plot(optClust)

## ----opticsvsdbscan------------------------------------------------------
suppressWarnings(opticsPoints <-	bind_rows(Map(f = function(x) {
		clust = thiscut <- optics_cut(optClust, x)$cluster
		data.frame(cluster = clust, eps = x)
	}, c(1, 3, 5))))
opticsPoints$cluster <- factor(opticsPoints$cluster)
opticsPoints$x <- rep(dat$x, 3)
opticsPoints$y <- rep(dat$y, 3)
opticsPoints$eps <- factor(paste("epsilon ==" , opticsPoints$eps, sep=""))

ggplot(data = opticsPoints, aes(x = x, y = y, color = cluster)) +
	geom_point(size = 0.5, alpha = 0.7) +
	facet_grid(. ~ eps, labeller = label_parsed) + 
	scale_x_continuous("", breaks = NULL) +
	scale_y_continuous("", breaks = NULL) +
	guides(color = FALSE) +
	ggtitle("OPTICS Clusters With EPS Cuts")

## ----hdbscan-------------------------------------------------------------
suppressWarnings(set <- rbind(Map(f = function(y) {
	rbind(Map(f = function(x) {
		hdclust <- hdbscan(edges = edges, neighbors = neighbors, K = y, minPts = x, threads = 1)$cluster
		data.frame(cluster = hdclust, K = x, minPts = y)
	}, c(6, 10, 20)))
}, c(4, 6, 8))))
set <- lapply(set, FUN = bind_rows)
set <- bind_rows(set)
set$x <- rep(dat$x, 9)
set$y <- rep(dat$y, 9)
set$cluster <- factor(set$cluster)
set$K <- factor(paste("K=", set$K))
set$minPts <- factor(paste("minPts=", set$minPts))

ggplot(data = set, aes(x = x, y = y, color = cluster)) +
	geom_point(size = 0.5, alpha = 0.7) +
	facet_grid(minPts ~ K) + 
	scale_x_continuous("", breaks = NULL) +
	scale_y_continuous("", breaks = NULL) +
	guides(color = FALSE) +
	ggtitle("HDBSCAN Is Robust To Hyperparameter Changes")

