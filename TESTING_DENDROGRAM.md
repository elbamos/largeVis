# Testing as.dendrogram.hdbscan

This document explains how to test the `as.dendrogram.hdbscan` function.

## Running Unit Tests

The comprehensive test suite has been added to `tests/testthat/test_h_hdbscan.R` under the context "as.dendrogram".

### Run all tests
```r
# Install devtools if needed
install.packages("devtools")

# Run all tests
devtools::test()
```

### Run only dendrogram tests
```r
library(testthat)
library(largeVis)

# Run tests from a specific context
test_file("tests/testthat/test_h_hdbscan.R", filter = "as.dendrogram")
```

## Standalone Test Script

A standalone test script `test_dendrogram.R` is also available for manual testing:

```bash
Rscript test_dendrogram.R
```

This script tests:
- Creating dendrograms with `includeNodes = FALSE`
- Creating dendrograms with `includeNodes = TRUE`
- Verifying dendrogram attributes
- Testing with different K parameters
- Plotting dendrograms

## Test Coverage

The test suite covers:

1. **S3 Method Registration**: Verifies `as.dendrogram.hdbscan` is properly registered
2. **Return Type**: Ensures the function returns a valid dendrogram object
3. **Structure Validation**:
   - Correct member counts
   - Valid height values
   - Proper attributes (leaf, members, height, label)
4. **Both Modes**:
   - `includeNodes = FALSE`: Cluster hierarchy only
   - `includeNodes = TRUE`: Full tree with individual points
5. **HDBSCAN-specific Attributes**:
   - Cluster IDs
   - Stability scores
   - Selected status
   - Lambda values
6. **Plotting**: Verifies dendrograms can be plotted without errors
7. **Multiple Datasets**: Tests on both iris and spiral datasets
8. **Parameter Variations**: Tests with different K values

## Expected Test Results

All tests should pass. The key assertions are:

- `as.dendrogram` is recognized as an S3 method for hdbscan objects
- Returned objects inherit from the "dendrogram" class
- When `includeNodes = TRUE`, member count equals the number of data points
- Dendrograms can be plotted using base R `plot()` function
- HDBSCAN-specific attributes are preserved in the dendrogram nodes

## Manual Verification Example

```r
library(largeVis)
data(iris)

# Prepare data
dat <- as.matrix(iris[, 1:4])
vis <- largeVis(t(dat), K = 20, sgd_batches = 1, verbose = FALSE)

# Run HDBSCAN
hdbscanobj <- hdbscan(vis, minPts = 10, K = 5, verbose = FALSE)

# Create dendrogram
dend <- as.dendrogram(hdbscanobj)

# Verify attributes
str(attributes(dend))

# Plot
plot(dend, main = "HDBSCAN Dendrogram")

# With individual nodes
dend_full <- as.dendrogram(hdbscanobj, includeNodes = TRUE)
plot(dend_full, leaflab = "none", main = "Full HDBSCAN Dendrogram")
```

## Debugging Failed Tests

If tests fail, check:

1. **NAMESPACE**: Ensure `S3method(as.dendrogram,hdbscan)` is present
2. **Function Location**: Verify `R/dendrogram.R` is in the package
3. **Dependencies**: Make sure the package is properly installed with `devtools::install()`
4. **Data**: Ensure test data files in `inst/testdata/` are accessible

## Continuous Integration

These tests will run automatically as part of the package's CI/CD pipeline on:
- Travis CI
- GitHub Actions (if configured)
- R CMD check

To run R CMD check locally:
```bash
R CMD build .
R CMD check largeVis_*.tar.gz
```
