context("LOF")

set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
K <- 80
neighbors <- randomProjectionTreeSearch(dat, K = K, verbose = FALSE)


test_that(paste("LOF is consistent", 20), {
	load(system.file("testdata/truelof20.Rda", package = "largeVis"))
	edges <- neighbors$edgematrix
	ourlof <- lof(edges)
	expect_lt(sum(truelof20 - ourlof)^2 / ncol(dat), 0.4)
})

context("hdbscan")

test_that("hdbscan finds 3 clusters and outliers in spiral with a large Vis object", {
	load(system.file("testdata/spiral.Rda", package = "largeVis"))
	expect_silent(clustering <- hdbscan(spiral, K = 3, minPts = 20))
	expect_equal(length(unique(clustering$clusters)), 3)
	expect_equal(sum(as.numeric(clustering$clusters) < 1, na.rm = TRUE), 0)
})

test_that("hdbscan finds outliers", {
	load(system.file("testdata/spiral.Rda", package = "largeVis"))
	expect_silent(clustering <- hdbscan(spiral, K = 10, minPts = 20))
	expect_true(any(is.na(clustering$clusters)))
})

test_that("hdbscan is fine with minpts < 6", {
	load(system.file("testdata/spiral.Rda", package = "largeVis"))
	expect_silent(clustering <- hdbscan(spiral, K = 10, minPts = 3))
	expect_true(any(is.na(clustering$clusters)), 0)
})

test_that("hdbscan finds 3 clusters and outliers in spiral", {
	load(system.file("testdata/spiral.Rda", package = "largeVis"))
	expect_silent(clustering <- hdbscan(spiral, K = 3, minPts = 20))
	expect_equal(length(unique(clustering$clusters)), 3)
})


set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
K <- 20
neighbors <- randomProjectionTreeSearch(dat, K = K, verbose = FALSE)

test_that("hdbscan doesn't crash without 3 neighbors and is correct", {
	expect_silent(clustering <- hdbscan(neighbors = neighbors, minPts = 20, K = 3, verbose = FALSE))
	expect_equal(length(unique(clustering$clusters)), 3)
})

test_that("hdbscan doesn't crash on glass edges", {
	skip_on_travis()
	load(system.file("testdata/glassEdges.Rda", package = "largeVis"))
	expect_silent(clustering <- hdbscan(glassEdges, verbose = FALSE))
})

test_that("failing example doesn't fail", {
	data(iris)
	expect_silent(vis <- largeVis(t(iris[,1:4]), K = 20, sgd_batches = 1))
	expect_silent(hdbscanobj <- hdbscan(vis, minPts = 10, K = 5))
})

test_that("glosh is in range", {
	data(iris)
	expect_silent(vis <- largeVis(t(iris[,1:4]), K = 20, sgd_batches = 1))
	expect_silent(hdbscanobj <- hdbscan(vis, minPts = 10, K = 5))
	expect_equal(sum(hdbscanobj$glosh < 0), 0)
	expect_equal(sum(hdbscanobj$glosh > 1), 0)
})

context("as.dendrogram")

set.seed(1974)
data(iris)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
dat <- dat[-dupes, ]
dat <- t(dat)
K <- 20
neighbors <- randomProjectionTreeSearch(dat, K = K, verbose = FALSE)
hdobj <- hdbscan(neighbors = neighbors, minPts = 10, K = 4, verbose = FALSE)

test_that("as.dendrogram returns a dendrogram object (includeNodes = FALSE)", {
	expect_silent(dend <- as.dendrogram(hdobj, includeNodes = FALSE))
	expect_true(inherits(dend, "dendrogram"))
	expect_true(!is.null(attr(dend, "members")))
	expect_true(!is.null(attr(dend, "height")))
	expect_true(!is.null(attr(dend, "leaf")))
})

test_that("as.dendrogram returns a dendrogram object (includeNodes = TRUE)", {
	expect_silent(dend <- as.dendrogram(hdobj, includeNodes = TRUE))
	expect_true(inherits(dend, "dendrogram"))
	expect_true(!is.null(attr(dend, "members")))
	expect_true(!is.null(attr(dend, "height")))
})

test_that("as.dendrogram has correct member count (includeNodes = TRUE)", {
	dend <- as.dendrogram(hdobj, includeNodes = TRUE)
	expect_equal(attr(dend, "members"), ncol(dat))
})

test_that("as.dendrogram attributes are correct", {
	dend <- as.dendrogram(hdobj, includeNodes = FALSE)
	# Check basic dendrogram attributes exist
	expect_true(!is.null(attr(dend, "members")))
	expect_true(!is.null(attr(dend, "height")))
	expect_true(!is.null(attr(dend, "label")))
	expect_true(is.numeric(attr(dend, "height")))
	expect_true(is.integer(attr(dend, "members")) || is.numeric(attr(dend, "members")))
})

test_that("as.dendrogram can be plotted", {
	dend <- as.dendrogram(hdobj, includeNodes = FALSE)
	expect_silent({
		pdf(file = NULL)
		plot(dend)
		dev.off()
	})
})

test_that("as.dendrogram handles spiral dataset", {
	load(system.file("testdata/spiral.Rda", package = "largeVis"))
	clustering <- hdbscan(spiral, K = 3, minPts = 20, verbose = FALSE)
	expect_silent(dend <- as.dendrogram(clustering, includeNodes = FALSE))
	expect_true(inherits(dend, "dendrogram"))
})