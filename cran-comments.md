## Resubmission
 
This is a resubmission of a new release. 

The last submission failed to compile on Windows.  This was a reversion of a fix I'd put in to deal with an upstream change that happened earlier in the week.  Very sorry about that.  

In this version I have:
* Bumped the version to 0.1.9
* Changed certain variable types to better support 32-bit systems.
* Added C++ unit tests.
* Added package startup messages if the package was compiled for 32-bit or without OpenMP
* Add regularization to certain functions, to handle an edge-case where a dataset contains a large number of duplicates.
* Implemented Pairing Heap for hdbscan implementation of Prim's algorithm.
* Worked around apparent bug in RcppArmadillo that caused intermittent seg faults when returning a sparse matrix from C++ to R.

* NOTE:  "largeVis" is the correct name of the package and function implementing the LargeVis algorithm.
 
## Test environments
* local OS X install, R 3.3.1
* ubuntu 14.04 (on travis-ci), R 3.3.1
* ubuntu 16.04 (locally) R 3.3.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* 1 NOTE relates to the inclusion of covr in "suggests".

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---
  
