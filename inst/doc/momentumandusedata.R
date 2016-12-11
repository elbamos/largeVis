## ----setupvignette,eval=T,echo=F,warning=F,error=F,message=F-------------
require(ggplot2, 
        quietly = TRUE)
library(png)
library(grid)
knitr::opts_chunk$set(collapse = TRUE, 
                      comment = "#>",
                      cache = FALSE, 
											echo = FALSE)
theme_set(
  theme_bw() %+replace%
  theme(
    legend.title = element_text(size = rel(0.8), face = "bold"),
    legend.margin = margin(0, "cm"),
    legend.position = "right",
    legend.key.size = unit(2, "lines"),
    legend.text = element_text(size = unit(8, "points")), 
    axis.title.y = element_text(angle = 90),
    axis.text = element_text(size = rel(0.7)),
    plot.margin = margin(0, 0.5, 1, 0, "lines"), 
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
}

## ----drawdegree----------------------------------------------------------
if (rebuild) {
	degreeplot <- ggplot(degreeCoords, aes(x = X1, y = X2, color = label)) + 
		geom_point(size = 0.02, alpha = 0.1) + 
		facet_grid(. ~ degree) +
		scale_x_continuous("", labels = NULL, breaks = NULL) + 
		scale_y_continuous("", labels = NULL, breaks = NULL) +
		guides(color = FALSE) +
		ggtitle("Effect of useDegree")
	ggsave(degreeplot, 
			 filename = system.file(package = "largeVis", "vignettedata/degreeplot.png"),
			 width = 5, height = 4)
}

## ----drawdegreeimage,fig.width=5,fig.height=4----------------------------
img <- readPNG(system.file(package = "largeVis", "vignettedata/degreeplot.png"))
grid.raster(img)

## ----momentum------------------------------------------------------------
if (rebuild) {
	starterCoords <- matrix(runif(n = 2 * ncol(mnist)) - 0.5, ncol = 2)
	firstCoords <- data.frame(
		t(projectKNNs(preVis$wij, coords = t(starterCoords), sgd_batches = 0.1, verbose = TRUE)), 
    lambda = 0, batches = 0.1, label = labels)
	for (batches in c(0.3, 0.8)) {
		newCoords <- data.frame(t(projectKNNs(preVis$wij, 
																					verbose = TRUE, 
																					sgd_batches = batches, 
																					coords = t(starterCoords))))
		newCoords$lambda <- 0
		newCoords$batches <- batches
		newCoords$label <- labels
		firstCoords <<- rbind(firstCoords, newCoords)
	}
	for (lambda in c(0.4, 0.9)) {
			for (batches in c(0.1, 0.3, 0.5)) {
				newtime <- system.time(newCoords <- data.frame(t(projectKNNs(preVis$wij, 
																																		 verbose = TRUE,
																																		 sgd_batches = batches, 
																																		 momentum = lambda, 
																																		 coords = t(starterCoords)))))
				newCoords$lambda <- lambda
				newCoords$batches <- batches
				newCoords$label <- labels
				firstCoords <<- rbind(firstCoords, newCoords)
			}
	}
	momentumCoords <- firstCoords
	momentumCoords$label <- factor(momentumCoords$label)
}

## ----drawmomentum,warning=FALSE------------------------------------------
if (rebuild) {
	momentumPlot <- ggplot(momentumCoords, aes(x = X1, y = X2, color = label)) + 
		geom_point(size = 0.01, alpha = 0.1) + 
		facet_grid(batches ~ lambda, scales = "free", labeller = label_bquote(cols = lambda == .(lambda), 
																																					rows = b == .(batches))) +
		scale_x_continuous("", limits = c(-40, 40), labels = NULL, breaks = NULL) + 
		scale_y_continuous("", limits = c(-40, 40), labels = NULL, breaks = NULL) +
		guides(color = FALSE) +
		ggtitle("Effect of Momentum and Reduced Training Batches")
	ggsave(momentumPlot, 
				 filename = system.file(package = "largeVis", "vignettedata/momentumplot.png"),
				 width = 6, height = 4)
}

## ----drawmomentumimage,fig.width=6,fig.height=4--------------------------
img <- readPNG(system.file(package = "largeVis", "vignettedata/momentumplot.png"))
grid.raster(img)

## ----dbscan,fig.width=6,fig.height=6-------------------------------------
load(system.file(package = "largeVis", "vignettedata/spiral.Rda"))
dat <- spiral
neighbors <- randomProjectionTreeSearch(t(dat), K = 20)
edges <- buildEdgeMatrix(t(dat), neighbors = neighbors)
set <- rbind(Map(f = function(y) {
	rbind(Map(f = function(x) {
		clust = lv_dbscan(edges = edges, neighbors = neighbors, eps = x, minPts = y)$cluster
		data.frame(cluster = clust, eps = x, minPts = y)
	}, c(1, 3, 5)))
}, c(5, 10, 20)))
lbind <- function(x) do.call(rbind, x)
set <- lapply(set, FUN = lbind)
set <- lbind(set)
set$x <- rep(dat[, 1], 9)
set$y <- rep(dat[, 2], 9)
set$cluster <- factor(set$cluster)
set$eps <- ordered(set$eps, labels = paste("epsilon == ", c(1, 3, 5), sep = ""))
set$minPts <- ordered(set$minPts, labels = paste("minPts == ", c(5, 10, 20), sep = ""))
ggplot(data = set, aes(x = x, y = y, color = cluster)) +
	geom_point(size = 0.5, alpha = 0.7) +
	facet_grid(minPts ~ eps, labeller = label_parsed) + 
	scale_x_continuous("", breaks = NULL) +
	scale_y_continuous("", breaks = NULL) +
	guides(color = FALSE) +
	ggtitle("Effect of eps and minPts on DBSCAN results")

## ----optics,fig.width=5,message=FALSE,warning=FALSE----------------------
library(dbscan, quietly = TRUE)
optClust <- lv_optics(edges = edges, neighbors = neighbors, eps = 5, useQueue = FALSE, minPts = 5)
optClust2 <- lv_optics(edges = edges, neighbors = neighbors, eps = 5, useQueue = TRUE, minPts = 5)
ggplot(data.frame(
	o = c(optClust$order, optClust2$order), 
	d = c(optClust$reachdist, optClust2$reachdist), 
	useQueue = rep(c("No Queue", "Use Queue"), each = length(optClust$order))
), aes(x = o, y = d)) + 
		geom_point(stat = 'identity', size = 0.1) + 
		geom_bar(stat = 'identity', alpha = 0.3) +
		facet_grid(useQueue ~ .) +
		scale_x_continuous("Order") + 
		scale_y_continuous("Reachability Distance")

## ----opticsvsdbscan,fig.width=2,fig.width=6------------------------------
suppressWarnings(opticsPoints <- do.call(rbind, Map(f = function(x) {
		clust = thiscut <- optics_cut(optClust, x)$cluster
		data.frame(cluster = clust, eps = x)
	}, c(1, 3, 5))))
opticsPoints$cluster <- factor(opticsPoints$cluster)
opticsPoints$x <- rep(dat[, 1], 3)
opticsPoints$y <- rep(dat[, 2], 3)
opticsPoints$eps <- factor(paste("epsilon ==" , opticsPoints$eps, sep = ""))

ggplot(data = opticsPoints, aes(x = x, y = y, color = cluster)) +
	geom_point(size = 0.5, alpha = 0.7) +
	facet_grid(. ~ eps, labeller = label_parsed) + 
	scale_x_continuous("", breaks = NULL) +
	scale_y_continuous("", breaks = NULL) +
	guides(color = FALSE) +
	ggtitle("OPTICS Clusters With EPS Cuts")

## ----hdbscan,fig.width=6,fig.height=6------------------------------------
suppressWarnings(set <- do.call(rbind, Map(f = function(y) {
	rbind(Map(f = function(x) {
		hdclust <- hdbscan(edges = edges, neighbors = neighbors, K = y, minPts = x)$cluster
		data.frame(cluster = as.numeric(hdclust), K = x, minPts = y)
	}, c(6, 10, 20)))
}, c(2, 6, 12))))
lbind <- function(x) do.call(rbind, x)
set <- lbind(set)
set$x <- rep(dat[, 1], 9)
set$y <- rep(dat[, 2], 9)
set$cluster <- factor(set$cluster)
set$K <- factor(paste("K=", set$K))
set$minPts <- factor(paste("minPts=", set$minPts))

ggplot(data = set, aes(x = x, y = y, color = cluster)) +
	geom_point(size = 0.5, alpha = 0.7) +
	facet_grid(minPts ~ K) + 
	scale_x_continuous("", breaks = NULL) +
	scale_y_continuous("", breaks = NULL) +
	guides(color = FALSE) +
	ggtitle("HDBSCAN Is Robust\nTo Hyperparameter Changes")

