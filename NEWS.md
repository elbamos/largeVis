### largeVis 0.1.6

* Neighbor search:
  + Dense search is much, much faster and more efficient
  + Tree search for cosine distances uses normalized vectors
* projectKNNs 
  + faster because of more efficient scanning for negative samples
  + Skips negative samples with already-large distances
* Vignettes:
  + Reuse initialization matrices and neighbors, to make it easier to see the effect of hyperparameters
  + Benchmarks now a separate vignette, more detailed
  + Added examples of manifold map with color faces using OpenFace vectors
* Sigmas
  + Fixed longstanding bug
  + Changed optimization method; pre-scan
  + Sigma estimation could, in some cases involving very large datasets, could become unstable and settle into 
  an edge of the interval. This should be resolved.  
* Visualization
  + Color manifold maps work
  + Ported Karpathy's function for non-overlapping embeddings (experimental)
  + Removed transparency parameter
  + Added ggManifoldMap function for adding a manifold map to a ggplot2 plot
* vis
  + Whether to return neighbors and sigmas now adjustable parameters, for memory reasons
  + Runs gc() periodically
* Dependencies & Build
  + Many misc changes to simplify dependencies for CRAN
  + Re-added ARMA_64BIT_WORD; otherwise, could exceed the limitation on size of an arma sparse matrix with moderately sized datasets (~ 1 M rows, K = 100)
  + Consolidated C++ code into a single file to reduce library size
  + Now depends on R >= 3.0.2, so RcppProgress and RcppArmadillo could be moved from the Depends section of the DESCRIPTION file
* Correctness and Testing
  + Tested against the paper authors' wiki-doc and wiki-word datasets
  + Tested with up to 2.5mm rows. 
  + Tests are separated by subject
  + Additional, more extensive tests with greater code coverage
* OpenMP
  + Will now compile on systems that lack OpenMP (e.g., OS X systems with old versions of xcode). 

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
