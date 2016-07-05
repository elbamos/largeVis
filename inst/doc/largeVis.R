## ----setupvignette,eval=T,echo=F,warning=F,error=F,message=F-------------
# Note to reader:  Please don't steal the semi-distinctive visual style I spent several minutes creating for myself.
require(ggplot2, 
        quietly = TRUE)
require(RColorBrewer, 
        quietly = TRUE)
require(wesanderson, 
        quietly = TRUE)
library(jpeg,
        quietly = TRUE)
knitr::opts_chunk$set(collapse = TRUE, 
                      comment = "#>")
colors_discrete <- function(x) rep(wes_palette("Darjeeling", 
                                               n = min(x, 5)), 
                                   2)[1:x]
colors_divergent_discrete <- function(x) 
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(x, "Spectral"))
colors_continuous <-  function(x) wes_palette(name = "Zissou",
                                              n = x, 
                                              type = "continuous")

nacol <- colors_discrete(4)[4]
theme_set(
  theme_bw() %+replace%
  theme(
    legend.key.size = unit(4, "mm"), 
    legend.title = element_text(size = rel(0.8),
                              face = "bold"),
    legend.margin = unit(0, "cm"),
    legend.position = "bottom",
    legend.key.size = unit(0.5, "lines"),
    legend.text=element_text(size = unit(8, "points")), 
    axis.title.y = element_text(angle = 90),
    axis.text = element_text(size = rel(0.7)),
    plot.margin = unit(c(0, 0.5, 1, 0), "lines"), 
    axis.title = element_text(size = rel(0.8),
                              face = "bold"),
    title = element_text(size = rel(0.9))
  ) 
)
rebuild <- FALSE

require(largeVis,quietly = TRUE)

## ----buildhyperparameters,echo=F,eval=rebuild----------------------------
#  data(wiki)
#  
#  inputs <- data.frame(
#    g = rep(c(.5,1,7,14), 5),
#    a = rep(c(0,.1,1,5,10), each = 4)
#  )
#  set.seed(1974)
#  initialcoords <- matrix(rnorm(ncol(wiki) * 2), nrow = 2)
#  
#  agcoords <- do.call(rbind,
#                      lapply(1:nrow(inputs),
#                             FUN = function(x) {
#    a <- inputs[x, 'a']
#    g <- inputs[x, 'g']
#    newcoords <- initialcoords
#    localcoords <- projectKNNs(wiki,
#                               alpha =  a,
#                               gamma = g,
#                               verbose = FALSE,
#                               coords = newcoords)
#    localcoords <- data.frame(scale(t(localcoords)))
#    colnames(localcoords) <- c("x", "y")
#    localcoords$a <- a
#    localcoords$g <- g
#    localcoords$activity <- log(Matrix::colSums(wiki))
#    localcoords
#  }))

## ----drawhyperparameters,echo=F,fig.width=3.5,fig.height=4,fig.align='center',results='asis'----
load(system.file("extdata/agcoords.Rda", package = "largeVis"))
ggplot(agcoords,
       aes(x = x, 
           y = y, 
           color = activity)) +
  geom_point(alpha = 0.2, 
             size = 0.05) +
  facet_grid(a ~ g,
             labeller = label_bquote(alpha == .(a), 
                                     gamma == .(g)),
             scales = 'free') +
  scale_x_continuous(breaks = NULL, 
                     name = "") +
  scale_y_continuous(breaks = NULL, 
                     name = "") +
  scale_color_gradientn(colors = colors_continuous(10), 
                        guide=FALSE) +
  ggtitle(expression(paste("Effect of ", alpha, " vs. ", gamma, sep = "  ")))

## ----iris_mkhyperparams,echo=F,eval=rebuild------------------------------
#  data(iris)
#  Ks <- c(5, 10,20,30)
#  Ms <- c(5, 10, 20)
#  dat <- iris[,1:4]
#  dupes <- duplicated(dat)
#  dat <- dat[-dupes,]
#  labels <- iris$Species[-dupes]
#  dat <- as.matrix(dat)
#  dat <- t(dat)
#  
#  set.seed(1974)
#  coordsinput <- matrix(rnorm(ncol(dat) * 2), nrow = 2)
#  neighbors <- randomProjectionTreeSearch(dat,
#                                        K = max(Ks),
#                                        verbose = FALSE)
#  
#  iriscoords <- do.call(rbind, lapply(Ks, FUN = function(K) {
#    neighborIndices <- neighborsToVectors(neighbors[1:K,])
#    distances <- largeVis::distance(x = dat,
#                                    neighborIndices$i,
#                                    neighborIndices$j,
#                                    verbose = FALSE)
#    wij <- buildEdgeMatrix(i = neighborIndices$i,
#                         j = neighborIndices$j,
#                         d = distances, verbose = FALSE)
#    do.call(rbind, lapply(Ms, FUN = function(M) {
#      coords <- projectKNNs(wij = wij$wij, M = M,
#                            coords = coordsinput,
#                            verbose = FALSE)
#      coords <- scale(t(coords))
#      coords <- data.frame(coords)
#      colnames(coords) <- c("x", "y")
#      coords$K <- K
#      coords$M <- M
#      coords$Species <- as.integer(labels)
#      coords
#    }))
#  }))
#  iriscoords$Species <- factor(iriscoords$Species)
#  levels(iriscoords$Species) <- levels(iris$Species)

## ----drawiris,echo=F,fig.width=4,fig.height=4.5,fig.align='center',results='asis'----
load(system.file("extdata/iriscoords.Rda", package = "largeVis"))

ggplot(iriscoords,
       aes(x = x,
           y = y,
           color = Species)) +
         geom_point(size = 0.5) +
  scale_x_continuous("", 
                     breaks = NULL) +
  scale_y_continuous("", 
                     breaks = NULL) +
  facet_grid(K ~ M, 
             scales = 'free', 
             labeller = label_bquote(K == .(K), M == .(M))) +
  scale_color_manual(values = colors_discrete(3)) +
  ggtitle("Effect of M and K on Iris Dataset")

## ----echomanifold,echo=T,eval=F------------------------------------------
#  dim(trainData) <- c(60000, 28, 28)
#  aperm(trainData, perm = c(1,3,2), resize = FALSE)
#  set.seed(1974)
#  manifoldMap(mnistCoords[,1:2],
#      n = 5000,
#      scale = 0.1,
#      images = trainData,
#      xlab = "",
#      ylab = "")

## ----youtube,eval=F,echo=T-----------------------------------------------
#  pathToGraphFile <-
#    "./YouTubeCommunities/com-youtube.ungraph.txt"
#  pathToCommunities <-
#    "./YouTubeCommunities/com-youtube.top5000.cmty.txt"
#  
#  youtube <- readr::read_tsv(pathToGraphFile, skip=4, col_names=FALSE)
#  youtube <- as.matrix(youtube)
#  youtube <- Matrix::sparseMatrix(i = youtube[, 1],
#                                  j = youtube[, 2],
#                                  x = rep(1, nrow(youtube)),
#                                  dims = c(max(youtube), max(youtube)))
#  youtube <- youtube + t(youtube)
#  communities <- readr::read_lines(pathToCommunities)
#  communities <- lapply(communities,
#                        FUN = function(x) as.numeric(unlist(strsplit(x, "\t"))))
#  community_assignments <- rep(0,
#                               nrow(youtube))
#  for (i in 1:length(communities)) community_assignments[communities[[i]]] <- i
#  
#  youTube_coordinates <- projectKNNs(youtube)
#  youTube_coordinates <- data.frame(scale(t(youTube_coordinates)))
#  colnames(youTube_coordinates) <- c("x", "y")
#  youTube_coordinates$community <- factor(community_assignments)
#  youTube_coordinates$alpha <- factor(ifelse(youTube_coordinates$community == 0, 0.05, 0.2))
#  ggplot(youTube_coordinates, aes( x = x,
#                        y = y,
#                        color = community,
#                        alpha = alpha,
#                        size = alpha)) +
#    geom_point() +
#    scale_color_manual(values =
#                         c("black", colors_continuous(5000)),
#                       guide = FALSE) +
#    scale_alpha_manual(values = c(0.005, 0.2), guide = FALSE) +
#    scale_size_manual(values = c(0.03, 0.15), guide = FALSE) +
#    scale_x_continuous("",
#                       breaks = NULL, limits = c(-2.5,2.5)) +
#    scale_y_continuous("",
#                       breaks = NULL, limits = c(-2.5,2.5)) +
#    ggtitle("YouTube Communities")

## ----lowmemexample,eval=F,echo=T-----------------------------------------
#  neighbors <- randomProjectionTreeSearch(largeDataset)
#  neighborIndices <- neighborsToVectors(neighbors)
#  rm(neighbors)
#  gc()
#  distances <- distance(x = largeDataset,
#                        i = neighborIndices$i,
#                        j =neighborIndices$j)
#  rm(largeDataset)
#  gc()
#  wij <- buildEdgeMatrix(i = neighborIndices$i,
#                         j = neighborIndices$j,
#                         d = distances)
#  rm(distances, neighborIndices)
#  gc()
#  coords <- projectKNNs(wij$wij)

