## ----setup,eval=T,echo=F,warning=F,error=F,message=F---------------------
# Note to reader:  Please don't steal the semi-distinctive visual style I spent several minutes creating for myself.
library(RColorBrewer,quietly=T)
library(wesanderson,quietly=T)
colors_discrete <- function(x) rep(wes_palette("Darjeeling", n = min(x,5)), 
                                   2)[1:x]
colors_divergent_discrete <- function(x) grDevices::colorRampPalette(RColorBrewer::brewer.pal(x, "Spectral"))
colors_continuous <-  function(x) wes_palette(name= "Zissou",n = x, type= "continuous")

nacol <- colors_discrete(4)[4]
require(ggplot2,quietly = T)
theme_set(
  theme_bw() %+replace%
  theme(
    legend.key.size=unit(4,"mm"), 
    legend.title=element_text(size=rel(0.8), face = "bold"),
    legend.margin=unit(0,"cm"),
    legend.position="bottom",
    legend.key.size=unit(0.5,"lines"),
    legend.text=element_text(size = unit(8, "points")), 
    axis.title.y = element_text(angle=90),
    axis.text = element_text(size=rel(0.7)),
    plot.margin = unit(c(0, 0.5, 1, 0), "lines"), 
    axis.title=element_text(size=rel(0.8),face="bold"),
  title = element_text(size=rel(0.9))
                          ) 
)
require(largeVis)

## ----MNIST,echo=F,message=F,warning=F,results='hide',eval=F--------------
#  darch::provideMNIST(download=T)
#  load("data/train.RData")
#  
#  mnistCoords <- vis(t(trainData) - 0.5, K = 40, tree_threshold = 700,
#                     n_trees = 40, max_iter = 2, verbose=F)
#  mnistCoords <- mnistCoords$coords
#  mnistCoords <- scale(t(mnistCoords))
#  mnistCoords <- data.frame(mnistCoords)
#  colnames(mnistCoords) <- c("x", "y")
#  labs <- apply(trainLabels, MARGIN=1, FUN=function(x) which(x == 1))
#  mnistCoords$labels <- factor(labs - 1)

## ----drawmnist,echo=F,warning=F,fig.width=3.5,fig.height=4,fig.align='center',fig.show='hold'----
load(system.file("extdata", "mnistcoords.Rda", package="largeVis"))
ggplot(mnistCoords, aes(x = x, y = y, color = labels)) +
  geom_point(size = 0.1, alpha = 0.3) +
  scale_x_continuous(name = "", limits = c(-2.5, 2), breaks = NULL) +
  scale_y_continuous(name = "", limits = c(-2, 2.5), breaks = NULL) +
  scale_color_manual(values = colors_divergent_discrete(10)(10)) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  ggtitle("MNIST")

## ----ldafromldavis,echo=F,eval=F-----------------------------------------
#  library(LDAvis)
#  data("TwentyNewsgroups")
#  theta <- scale(t(TwentyNewsgroups$theta))
#  visObj <- vis(theta, K = 100, n_trees = 20, tree_threshold = 100,
#                max_iter = 2)
#  
#  ngcoords <- scale(t(visObj$coords))
#  ngcoords <- data.frame(ngcoords)
#  colnames(ngcoords) <- c("x", "y")
#  library(lda)
#  data("newsgroup.train.labels")
#  ngcoords$label <- factor(newsgroup.train.labels)[-1]

## ----draw20ng,fig.align='center',echo=F,fig.width=3.5,fig.height=4,eval=T,warning=FALSE,error=FALSE,message=FALSE,fig.show='hold'----
load(system.file("extdata", "ngcoords.Rda", package="largeVis"))
ggplot(ngcoords, 
       aes(x = x, y = y, color = label)) +
  geom_point(size = 0.4, alpha = 0.5) + 
  scale_color_manual(values = colors_divergent_discrete(20)(20),
                     guide=FALSE) +
  scale_x_continuous(name = "", limits = c(-2, 2.5), breaks = NULL) +
  scale_y_continuous(name = "", limits = c(-2, 2.5), breaks = NULL) +
  ggtitle("20 Newsgroups")

## ----3draw,webgl=TRUE,echo=F,eval=F,results='asis'-----------------------
#  # d3coords <- projectKNNs(visObj$wij, dim = 3)
#  # d3coords <- data.frame(scale(t(d3coords)))
#  # colnames(d3coords) <- c("x", "y", "z")
#  # d3coords$label <- factor(newsgroup.train.labels)[-1]
#  # library(threejs)
#  # rgl::plot3d(x = d3coords[,1],
#  #             y = d3coords[,2],
#  #             z = d3coords[,3],
#  #            main = "20 Newsgroups",
#  #            type = "p",
#  #            col = c(newsgroup.train.labels,
#  #                      newsgroup.test.labels))

## ----performance,echo=F,eval=F-------------------------------------------
#  benchmark <- readr::read_csv(system.file("extdata", "results.csv", package="largeVis"))
#  colnames(benchmark) <- c("time",
#                         "precision",
#                         "n_trees",
#                         "max_iters",
#                         "threshold")
#  benchmark$series <- factor(paste(benchmark$n_trees, "trees,",
#                                 benchmark$max_iters, "iterations."))

## ----plotpeformance,echo=F,fig.width=3.5,fig.height=4,fig.align='center'----
load(system.file("extdata", "benchmark.Rda", package = "largeVis"))
ggplot(benchmark, aes(x = time, y = precision / 100, 
                    group = series, color = series, 
                    shape = series,
                    label =threshold)) +
  geom_point(size = 1) + geom_line(size = 0.5) + 
  geom_text(vjust = 1, hjust = -0.1, size = 2.5) +
  scale_x_continuous("Time (relative)") + 
  scale_y_log10("Precision", limits = c(0.1,1), 
                breaks = c(.1, .25, .5, .8, .9, .99)) +
  scale_color_manual(values =      colors_divergent_discrete(nlevels(benchmark$series))(nlevels(benchmark$series))) +
  guides(color = guide_legend(nrow=3)) +
  ggtitle(expression(
    atop("Time vs. Precision (K = 1000)",
         atop(italic("Labelled by Tree Threshold"))
         )
    ))

## ----wikihyperparameters,echo=F,eval=F-----------------------------------
#  data(wiki)
#  
#  inputs <- data.frame(
#    g = rep(c(.5,1,7,14), 4),
#    a = rep(c(.1,1,5,10), each = 4)
#  )
#  
#  agcoords <- do.call(rbind, lapply(1:nrow(inputs), FUN = function(x) {
#    a <- inputs[x, 'a']
#    g <- inputs[x, 'g']
#    localcoords <- projectKNNs(wiki, alpha =  a, gamma = g,verbose=FALSE)
#    localcoords <- data.frame(scale(t(localcoords)))
#    colnames(localcoords) <- c("x", "y")
#    localcoords$a <- a
#    localcoords$g <- g
#    localcoords$activity <- log(Matrix::colSums(wiki))
#    localcoords
#  }))

## ----drawhyperparameters,echo=F,fig.width=3.5,fig.height=4,fig.align='center'----
load(system.file("extdata", "agcoords.Rda", package="largeVis"))
ggplot(agcoords,
       aes(x = x, y = y, color = activity)) +
  geom_point(alpha = 0.2, size = 0.05) +
  facet_grid(a ~ g,
             labeller = label_bquote(alpha == .(a), gamma == .(g)),
             scales = 'free') +
  scale_x_continuous(breaks=NULL,name="") +
  scale_y_continuous(breaks=NULL,name = "") +
  scale_color_gradientn(colors = colors_continuous(10), guide=FALSE) +
  ggtitle(expression(paste("Effect of", alpha, "vs.", gamma, sep = "  ")))

## ----iris,echo=F,fig.width=5,fig.height=5,eval=F-------------------------
#  data(iris)
#  Ks <- c(5, 10, 20, 40)
#  Ms <- c(1, 5, 10, 20)
#  data(iris)
#  dat <- iris[,1:4]
#  dupes <- duplicated(dat)
#  dat <- dat[-dupes,]
#  labels <- iris$Species[-dupes]
#  dat <- scale(dat)
#  dat <- as.matrix(dat)
#  dat <- t(dat)
#  
#  inputs <- data.frame(
#    K = rep(Ks, length(Ms)),
#    M = rep(Ms, each = length(Ks))
#  )
#  iriscoords <- do.call(rbind, lapply(1:nrow(inputs), FUN = function(x) {
#    K <- inputs[x, 'K']
#    M <- inputs[x, 'M']
#    visO <- vis(dat, K = K, M = M, verbose=FALSE)
#    localcoords <- data.frame(scale(t(visO$coords)))
#    colnames(localcoords) <- c("x", "y")
#    localcoords$K <- K
#    localcoords$M <- M
#    localcoords$Species <- as.integer(labels)
#    localcoords
#    }))
#  iriscoords$Species <- factor(iriscoords$Species)
#  levels(iriscoords$Species) <- levels(iris$Species)

## ----drawiriscoords,echo=F,fig.width=4,fig.height=4.5,fig.align='center'----
load(system.file("extdata", "iriscoords.Rda", package="largeVis"))
ggplot(iriscoords,
       aes(x = x,
           y = y,
           color =Species)) +
         geom_point(size = 0.5) +
  scale_x_continuous("", breaks = NULL) +
  scale_y_continuous("", breaks = NULL) +
  facet_grid(K ~ M, scales = 'free', labeller = label_bquote(K == .(K), M == .(M))) +
  scale_color_manual(values = colors_discrete(3)) +
  ggtitle("Effect of M and K on Iris Dataset")

## ----loadmnistimages,eval=F,echo=F---------------------------------------
#  load("data/train.RData")

## ----drawmanifoldmap,echo=T,fig.width=8,fig.height=8,message=F,warning=F,fig.align='center'----
if (exists("trainData")) {
  dim(trainData) <- c(60000, 28, 28)
  manifoldMap(mnistCoords[,1:2],
      n = 5000,
      scale = 0.003,
      transparency = F,
      images = trainData,
      xlab="", ylab="",
      xlim = c(-2, 2),
      ylim = c(-2, 2))
} 

## ----tdm,echo=F,eval=F---------------------------------------------------
#  library(stm)
#  data("poliblog5k")
#  p <- c(0, cumsum(as.numeric(lapply(poliblog5k.docs, function(x) ncol(x)))))
#  i <- do.call("c", lapply(poliblog5k.docs, function(x) x[1,]))
#  p[length(p)] <- length(i)
#  j <- rep(0:(length(diff(p)) - 1), diff(p))
#  v <- do.call("c", lapply(poliblog5k.docs, function(x) x[2,]))
#  poli <- Matrix::sparseMatrix(i = i + 1, j = j + 1, x = v)
#  dupes <- duplicated(slam::as.simple_triplet_matrix(Matrix::t(poli)))
#  poli <- poli[, ! dupes]
#  poli <- poli / log(Matrix::rowSums(poli > 0)) # tf-idf weight
#  policoords <- vis(poli, K = 100, n_trees = 20,
#                tree_threshold = 100, max_iter = 10,
#                M=10,gamma=15,
#              distance_method = 'Cosine',verbose=F)
#  polidata <- data.frame(scale(t(policoords$coords)))
#  colnames(polidata) <- c('x', 'y')
#  polidata$rating <- poliblog5k.meta$rating[!dupes]
#  polidata$blog <- poliblog5k.meta$blog[!dupes]

## ----drawtdm,echo=F,fig.height=4,fig.width=7-----------------------------
load(system.file("extdata", "polidata.Rda", package="largeVis"))
ggplot(polidata, aes(x = x, y = y, color = blog)) +
  geom_point(size = 0.3, alpha = 0.8) +
  scale_color_manual(values = colors_divergent_discrete(6)(6)) +
  facet_grid(. ~ rating, scale = 'free') +
  scale_x_continuous("", breaks = NULL) +
  scale_y_continuous("", breaks = NULL) +
  ggtitle("Visualization of a tf-idf Matrix")

## ----eval=F,echo=T-------------------------------------------------------
#  neighbors <- randomProjectionTreeSearch(largeDataset)
#  neighborIndices <- neighborsToVectors(neighbors)
#  rm(neighbors)
#  distances <- distance(neighborIndices$i,
#                        neighborIndices$j,
#                        largeDataset)
#  rm(largeDataset)
#  wij <- buildEdgeMatrix(i = neighborIndices$i,
#                         j = neighborIndices$j,
#                         d = distances)
#  rm(distances, neighborIndices)
#  coords <- projectKNNs(wij$wij)

