## Submission

This is the first submission of a minor version update to the package.  Changes in this update include:
* Fixes for bugs that emerged after 0.1.9.1 was accepted.
* Performance improvements.
* The addition of momentum to the stochastic gradient descent routine.
* Reimplementation of the dbscan and optics algorithms.
* Addition of a number of helper functions intended for compatibility with other R packages. E.g., there is now an `as.dendrogram` function for `hdbscan` objects; edge matrices can be converted to `dist` objects; etc. 
* An additional vignette covering momentum. 

## Test environments
* local OS X install, R 3.3.1
* ubuntu 12.04 (on travis-ci), R 3.3.1
* ubuntu 14.04 (on travis-ci), R 3.3.1 and R-devel
* OS X (on travis-ci), R-devel
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

## Reverse dependencies

There are no reverse dependencies reported on CRAN. 