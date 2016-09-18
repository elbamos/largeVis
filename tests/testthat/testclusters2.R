context("dbscan")

set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
K <- 149
neighbors <- randomProjectionTreeSearch(dat, K = K,  threads = 2, verbose = FALSE)
edges <- buildEdgeMatrix(data = dat,
												 neighbors = neighbors,
												 verbose = FALSE)

test_that("dbscan doesn't crash on iris", {
	expect_silent(dbscan(edges = edges, neighbors = neighbors, eps = 1, minPts = 10, verbose = FALSE))
})

test_that("dbscan matches iris", {
	load(system.file(package = "largeVis", "extdata/irisdbscan.Rda"))
	dbclusters <- dbscan(edges = edges, neighbors = neighbors, eps = 1, minPts = 10, verbose = FALSE)
	expect_lt(sum(dbclusters$cluster != irisclustering$cluster), 3)
})

test_that("dbscan matches dbscan when the neighborhoods are complete", {
	load(system.file(package = "largeVis", "extdata/jaindata.Rda"))
	jainclusters <- dbscan(edges = jaindata$edges, neighbors = jaindata$neighbors, eps = 2, minPts = 5, verbose = FALSE)
	expect_equal(jainclusters$cluster, jaindata$dbclusters)
})

context("optics")

test_that("optics doesn't crash on iris", {
  expect_silent(optics(edges = edges, neighbors = neighbors, eps = 10, minPts = 10, verbose = FALSE))
})

test_that("optics matches optics", {
	load(system.file(package = "largeVis", "extdata/irisoptics.Rda"))
	opclusters <- optics(edges = edges, neighbors = neighbors, eps = 1, minPts = 10, verbose = FALSE)
	expect_equal(opclusters$order, irisoptics$order)
	expect_equal(opclusters$reachdist, irisoptics$reachdist)
	expect_equal(opclusters$coredist, irisoptics$coredist)
})

test_that("optics matches optics on jain when neighbors are complete", {
	load(system.file(package = "largeVis", "extdata/irisoptics.Rda"))
	opclusters <- optics(edges = edges, neighbors = neighbors, eps = 1, minPts = 10, verbose = FALSE)
	expect_equal(opclusters$order, irisoptics$order)
	expect_equal(opclusters$reachdist, irisoptics$reachdist)
	expect_equal(opclusters$coredist, irisoptics$coredist)
})

test_that("optics matches optics when the neighborhoods are complete", {
	load(system.file(package = "largeVis", "extdata/jaindata.Rda"))
	jainclusters <- optics(edges = jaindata$edges, neighbors = jaindata$neighbors, eps = 4, minPts = 5, verbose = FALSE)
	expect_equal(jainclusters$coredist, jaindata$optics$coredist)
	expect_equal(jainclusters$order, jaindata$optics$order)
	expect_equal(jainclusters$reachdist, jaindata$optics$reachdist)
})

context("LOF")

set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
K <- 80
neighbors <- randomProjectionTreeSearch(dat, K = K,  threads = 2, verbose = FALSE)
edges <- buildEdgeMatrix(data = dat,
												 neighbors = neighbors,
												 verbose = FALSE)

test_that(paste("LOF is consistent", 20), {
	load(system.file("extdata/truelof20.Rda", package = "largeVis"))
	ourlof <- largeVis:::lof(edges)
	expect_lt(sum(truelof20 - ourlof)^2 / ncol(dat), 0.4)
})

test_that("LOF is consistent 10", {
	edges <- buildEdgeMatrix(data = dat,
													 neighbors = neighbors[1:10,],
													 verbose = FALSE)
	load(system.file("extdata/truelof10.Rda", package = "largeVis"))
	ourlof <- lof(edges)
	expect_lt(sum(truelof10 - ourlof)^2 / ncol(dat), 0.4)
})

context("hdbscan")

test_that("hdbscan finds 3 clusters and outliers in spiral", {
	load(system.file("extdata/spiral.Rda", package = "largeVis"))
	clustering <- hdbscan(spiral, K = 3, minPts = 20, threads = 1)
	expect_equal(length(unique(clustering$clusters)), 4)
})

set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
K <- 20
neighbors <- randomProjectionTreeSearch(dat, K = K,  threads = 2, verbose = FALSE)
edges <- buildEdgeMatrix(data = dat, neighbors = neighbors, verbose = FALSE)

test_that("hdbscan doesn't crash with neighbors", {
	tst <- hdbscan(edges, minPts = 10, neighbors = neighbors, K = 3, threads = 2, verbose = FALSE)
})

test_that("hdbscan is correct", {
  clustering <- hdbscan(edges, minPts = 10, K = 3,  threads = 2, verbose = FALSE)
  expect_equal(length(unique(clustering$clusters)), 3)
})

test_that("hdbscan is less correct with neighbors", {
  clustering <- hdbscan(edges, neighbors = neighbors, minPts = 10, K = 3,  threads = 2, verbose = FALSE)
  expect_equal(length(unique(clustering$clusters)), 2)
})

test_that("hdbscan doesn't crash on glass edges", {
	load(system.file("extdata/glassEdges.Rda", package = "largeVis"))
	clustering <- hdbscan(edges, threads = 2, verbose = FALSE)
	expect_equal(length(unique(clustering$clusters)), 3)
})

test_that("hdbscan doesn't crash on big bad edges", {
	skip("skipping big bad edges test because the data is too big for cran")
	load(system.file("extdata/kddneighbors.Rda", package = "largeVis"))
	load(system.file("extdata/kddedges.Rda", package = "largeVis"))
	expect_silent(clusters <- hdbscan(edges, neighbors = neighbors, threads = 2, verbose = FALSE))
})

test_that("hdbscan doesn't crash on big bad edges", {
	skip_on_cran()
	skip("skipping long test")
	load(system.file("extdata/badedges.Rda", package = "largeVis"))
	expect_silent(clustering <- hdbscan(badedges, threads = 2, verbose = FALSE))
})
