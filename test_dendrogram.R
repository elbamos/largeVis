#!/usr/bin/env Rscript
# Test script for as.dendrogram.hdbscan function
# Run this after installing the package to verify the function works correctly

library(largeVis)

cat("Testing as.dendrogram.hdbscan function\n")
cat("======================================\n\n")

# Test with iris dataset (similar to old tests)
cat("Loading iris dataset...\n")
data(iris)
dat <- as.matrix(iris[, 1:4])
dat <- scale(dat)
dupes <- which(duplicated(dat))
if (length(dupes) > 0) {
  dat <- dat[-dupes, ]
}
dat <- t(dat)

cat("Building neighbor graph...\n")
K <- 20
neighbors <- randomProjectionTreeSearch(dat, K = K, verbose = FALSE)
edges <- buildEdgeMatrix(data = dat, neighbors = neighbors, verbose = FALSE)

cat("Running HDBSCAN...\n")
hdobj <- hdbscan(edges, neighbors = neighbors, minPts = 10, K = 4, verbose = FALSE)

cat("\nHDBSCAN results:\n")
cat("  Number of clusters found:", length(unique(hdobj$clusters[!is.na(hdobj$clusters)])), "\n")
cat("  Number of outliers:", sum(is.na(hdobj$clusters)), "\n")
cat("  Hierarchy - number of cluster nodes:", length(hdobj$hierarchy$parent), "\n")

# Test 1: Create dendrogram without individual nodes
cat("\n--- Test 1: as.dendrogram(hdobj, includeNodes = FALSE) ---\n")
tryCatch({
  dend <- as.dendrogram(hdobj, includeNodes = FALSE)
  cat("SUCCESS: Dendrogram created\n")
  cat("  Class:", class(dend), "\n")
  cat("  Number of members:", attr(dend, "members"), "\n")
  cat("  Height:", attr(dend, "height"), "\n")
  cat("  Is leaf:", attr(dend, "leaf"), "\n")

  # Try to plot it
  cat("  Attempting to plot...\n")
  pdf(file = NULL)  # Suppress plot output
  plot(dend, main = "HDBSCAN Dendrogram (iris, without nodes)")
  dev.off()
  cat("  Plot SUCCESS\n")

}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
})

# Test 2: Create dendrogram with individual nodes
cat("\n--- Test 2: as.dendrogram(hdobj, includeNodes = TRUE) ---\n")
tryCatch({
  dend_full <- as.dendrogram(hdobj, includeNodes = TRUE)
  cat("SUCCESS: Full dendrogram created\n")
  cat("  Class:", class(dend_full), "\n")
  cat("  Number of members:", attr(dend_full, "members"), "\n")
  cat("  Expected members:", ncol(dat), "\n")

  if (attr(dend_full, "members") == ncol(dat)) {
    cat("  PASS: Member count matches number of data points\n")
  } else {
    cat("  WARNING: Member count doesn't match\n")
  }

  # Try to plot it
  cat("  Attempting to plot...\n")
  pdf(file = NULL)  # Suppress plot output
  plot(dend_full, main = "HDBSCAN Dendrogram (iris, with nodes)")
  dev.off()
  cat("  Plot SUCCESS\n")

}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
})

# Test 3: Check dendrogram attributes
cat("\n--- Test 3: Checking dendrogram attributes ---\n")
tryCatch({
  dend <- as.dendrogram(hdobj, includeNodes = FALSE)

  # Check if key attributes exist
  attrs_to_check <- c("members", "height", "leaf", "label", "class")
  for (attr_name in attrs_to_check) {
    val <- attr(dend, attr_name)
    if (!is.null(val)) {
      cat("  ", attr_name, ": ", val, "\n")
    } else {
      cat("  WARNING:", attr_name, "is NULL\n")
    }
  }

  # Check cluster-specific attributes
  cluster_attrs <- c("cluster", "stability", "selected")
  cat("  Cluster-specific attributes:\n")
  for (attr_name in cluster_attrs) {
    val <- attr(dend, attr_name)
    if (!is.null(val)) {
      cat("    ", attr_name, ": ", val, "\n")
    }
  }

}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
})

# Test 4: Test with different K parameter
cat("\n--- Test 4: Testing with K = 3 ---\n")
tryCatch({
  hdobj3 <- hdbscan(edges, neighbors = neighbors, minPts = 10, K = 3, verbose = FALSE)
  dend3 <- as.dendrogram(hdobj3, includeNodes = FALSE)
  cat("SUCCESS: Dendrogram created for K=3\n")
  cat("  Height:", attr(dend3, "height"), "\n")

  dend3_full <- as.dendrogram(hdobj3, includeNodes = TRUE)
  cat("SUCCESS: Full dendrogram created for K=3\n")
  cat("  Members:", attr(dend3_full, "members"), "\n")

}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
})

cat("\n======================================\n")
cat("Testing complete!\n")
cat("\nTo use the function in your own code:\n")
cat("  library(largeVis)\n")
cat("  # ... create hdbscan object ...\n")
cat("  dend <- as.dendrogram(hdbscanobj)\n")
cat("  plot(dend)\n")
