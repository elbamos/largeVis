### largeVis 0.1.6

* Revisions for CRAN release, including verifying correctness by reproducing paper examples, and timing tests/benchmarks
  + Tested against the paper authors' wiki-doc and wiki-word datasets
  + Tested with up to 2.5m rows, 100m edges (processed in 12 hours). 
* Neighbor search:
  + Dense search is much, much faster and more efficient
  + Tree search for cosine distances uses normalized vectors
* projectKNNs 
  + Should be 10x faster for small datasets
  + Replaced binary search ( O(n log n) ) with the alias algorithm for weighted sampling ( O(1) )
  + Clips and smooths gradients, per discussion with paper authors
  + Optimized implementation for alpha == 1
  + Removed option for mixing weights into loss function - doesn't make sense if gradients are being clipped. 
  + Fixed OpenMP-related bug which caused visualizations to be "fuzzy" if threads were scheduled in an
    unexpected order. 
* Vignettes:
  + Reuse initialization matrices and neighbors, to make it easier to see the effect of hyperparameters
  + Benchmarks now a separate vignette, more detailed
  + Examples removed from vignettes and moved to readme
  + Added examples of manifold map with color faces using OpenFace vectors
* Sigms, P_ij matrix, w_ij matrix
  + Replaced C++ code entirely with new code based on reference implementation 
  + Refactored R code into `buildEdgeMatrix()` and `buildWijMatrix()`, which are simpler. 
* Visualization
  + Color manifold maps work
  + Ported Karpathy's function for non-overlapping embeddings (experimental)
  + Removed transparency parameter
  + Added ggManifoldMap function for adding a manifold map to a ggplot2 plot
* vis
  + Whether to return neighbors and sigmas now adjustable parameters, for memory reasons
  + Runs gc() periodically
* Data
  + Removed most data and extdata that had been included before; this is to reduce size for CRAN submission
* Dependencies & Build
  + Many misc changes to simplify dependencies for CRAN
  + Re-added ARMA_64BIT_WORD; otherwise, could exceed the limitation on size of an arma sparse matrix with moderately sized datasets (~ 1 M rows, K = 100)
  + Now depends on R >= 3.0.2, so RcppProgress and RcppArmadillo could be moved from the Depends section of the DESCRIPTION file
  + Will now compile on systems that lack OpenMP (e.g., OS X systems with old versions of xcode). 
* Correctness and Testing
  + Tests are separated by subject
  + Additional, more extensive tests with greater code coverage
  + Added travis testing against OSX

### largeVis 0.1.5

* Handles substantially larger datasets
* Support for sparse matrices (for *much* larger datasets)
* Added better error reporting for tree search
* Handle situation in tree search where nodes are equidistant from the hyperplane
* Broke-out several components as separate functions, which makes a more-memory-efficient mode of operation possible
* Removed some unnecessary checking when processing neighbor graph
* RcppArmadillo 0.7.100.3.0 is now required (this was necessary for support for larger datasets)
* Added appveyor to check Windows compatibility

### largeVis 0.1.4

* Added option of Euclidean or Cosine distance. 
* Now using the median in random projection trees, to make splits more even. This should eliminate the need for the
max_depth parameter. 
* Benchmarks
* Vignette
* Vastly improved multi-threading performance

### largeVis 0.1.3

* Rewrote neighbor-search code to improve memory performance. 

* Performance of wji calculation improved by an order of magnitude.

### largeVis 0.1.2

* Added visualization function.

### largeVis 0.1.1

* Many changes because of OpenMP compatibility issues. 

* Moved from Rcpp to RcppArmadillo objects, etc. 

### largeVis 0.1.0

* Initial development releases.  Focused on correctness, performance, testing against larger datasets.

* Added a `NEWS.md` file to track changes to the package.
