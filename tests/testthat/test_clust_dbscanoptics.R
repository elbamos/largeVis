cutoptics <- function(x) {
	clust <- 1
	eps <- x$eps
	minPts <- x$minPts
	clusters <- rep(0, length(x$order))
	for (i in 1:length(x$order)) {
		thisNode <- x$order[i]
		if (is.na(x$predecessor[thisNode]) | x$reachdist[thisNode] > eps) clust <- clust + 1
		clusters[thisNode] <- clust
	}
	badClusters <- which(tabulate(clusters) < minPts)
	clusters[clusters %in% badClusters] <- NA
	remainders <- unique(clusters)
	match(clusters, remainders)
}

getDefects <- function(x, y, threshold) {
	if (length(x) != length(y)) return("bad lengths")
	dif <- (x - y)^2
	under <- dif > threshold
	sum(under, na.rm = TRUE)
}

context("dbscan-iris")

set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
K <- 147
neighbors <- randomProjectionTreeSearch(dat, K = K, tree_threshold = 80, n_trees = 10,  max_iter = 4, threads = 2, verbose = FALSE)
edges <- buildEdgeMatrix(data = dat,
												 neighbors = neighbors,
												 verbose = FALSE)

test_that("dbscan doesn't crash on iris", {
	expect_silent(lv_dbscan(edges = edges, neighbors = neighbors, eps = 1, minPts = 10, verbose = FALSE))
})

test_that("dbscan matches iris", {
	load(system.file(package = "largeVis", "extdata/irisdbscan.Rda"))
	dbclusters <- lv_dbscan(edges = edges, neighbors = neighbors, eps = 1, minPts = 10, verbose = FALSE)
	expect_lte(sum(dbclusters$cluster != irisclustering$cluster), 1)
})

context("dbscan-jain")

test_that("dbscan matches dbscan on jain when the neighborhoods are complete", {
	skip_on_cran()
	skip_on_travis()
	load(system.file(package = "largeVis", "extdata/jaindata.Rda"))
	jainclusters <- lv_dbscan(edges = jaindata$edges, neighbors = jaindata$neighbors, eps = 2.5, minPts = 10, verbose = FALSE)
	expect_equal(jainclusters$cluster, jaindata$dbclusters25$cluster)
})

context("optics-iris")

test_that("optics doesn't crash on iris", {
  expect_silent(lv_optics(edges = edges, neighbors = neighbors, eps = 10, minPts = 10, verbose = FALSE))
})

load(system.file(package = "largeVis", "extdata/irisoptics.Rda"))
opclusters <- lv_optics(edges = edges, neighbors = neighbors, eps = 1, minPts = 10, verbose = FALSE)

test_that("optics matches optics core infinities", {
	expect_equal(which(is.infinite(opclusters$coredist)), which(is.infinite(irisoptics$coredist)))
})

test_that("optics matches optics core dist not infinities", {
	expect_equal(opclusters$coredist[! is.infinite(opclusters$coredist)], irisoptics$coredist[! is.infinite(irisoptics$coredist)])
})

test_that("opticis iris cut to dbscan matches dbscan", {
	cl <- cutoptics(opclusters)
	dbclusters <- lv_dbscan(edges = edges, neighbors = neighbors, eps = 1, minPts = 10, verbose = FALSE)
	expect_equal(cl, dbclusters$cluster)
})

test_that("optics less than or equal to optics reachdist", {
	expect_equal(getDefects(opclusters$reachdist, irisoptics$reachdist, 1e-3), 0)
})

context("optics-jain")

load(system.file(package = "largeVis", "extdata/jaindata.Rda"))
jainclusters <- lv_optics(edges = jaindata$edges,
													neighbors = jaindata$neighbors,
													eps = 2.5, minPts = 5,
													verbose = FALSE)

test_that("optics matches optics core on jain when the neighborhoods are complete", {
	expect_equal(is.infinite(jainclusters$coredist), is.infinite(jaindata$optics$coredist))
	selections <- ! is.infinite(jainclusters$coredist) & ! is.infinite(jaindata$optics$coredist)
	expect_equal(jainclusters$coredist[selections], jaindata$optics$coredist[selections])
})

test_that("optics matches dbscan on jain when the neighborhoods are complete", {
	cl <- cutoptics(jainclusters)
	print(table(cl))
	print(table(jaindata$dbclusters5$cluster))
	expect_equal(cl, jaindata$dbclusters5$cluster + 1)
})

test_that("optics matches optics reachdist on jain when the neighborhoods are complete", {
	expect_equal(is.infinite(jainclusters$reachdist), is.infinite(jaindata$optics$reachdist))
	selections <- complete.cases(jainclusters$reachdist) & complete.cases(jaindata$optics$reachdist)
	expect_equal(getDefects(jainclusters$reachdist[selections], jaindata$optics$reachdist[selections], 1e-3), 0)
})

context("optics-elki")

load(system.file("extdata/opttest.Rda", package = "largeVis"))

x <- opttest$test_data
neighbors <- randomProjectionTreeSearch(t(opttest$test_data), K = 300, tree_threshold = 100, max_iter = 5, seed = 1974)
edges <- buildEdgeMatrix(t(opttest$test_data), neighbors = neighbors, threads = 1)

eps <- .1
#eps <- .06
eps_cl <- .1
minPts <- 10
res <- lv_optics(edges, neighbors, eps = eps,  minPts = minPts)

test_that("optics output format is correct", {
	expect_identical(length(res$order), nrow(x))
	expect_identical(length(res$reachdist), nrow(x))
	expect_identical(length(res$coredist), nrow(x))
	expect_identical(res$eps, eps)
	expect_identical(res$minPts, minPts)
})

test_that("optics result matches elki after cut", {
	optcut <- cutoptics(res)
	refcut <- cutoptics(opttest$elki)
	expect_equal(optcut, refcut)
})