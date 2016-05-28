
### largeVis 0.1.5

* Handles substantially larger datasets
* Support for sparse matrices (for *much* larger datasets)
* Added better error reporting for tree search
* Handle situation in tree search where nodes are equidistant from the hyperplane
* Broke-out several components as separate functions, which makes a more-memory-efficient mode of operation possible
* Removed some unnecessary checking when processing neighbor graph
* C++11 and RcppArmadillo 0.7.100.3.0 are now required (this was necessary for support for larger datasets)
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
