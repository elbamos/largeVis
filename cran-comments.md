## Resubmission
 
This is a resubmission of a new release. 

The original submission produced this error:  "ERROR Package required but not available: ‘dbscan’ "

In this version I have:
* Removed references to dbscan in the tests, replacing them with a new data file of test data. 
* Added implementations of the OPTICS and HDBSCAN clustering algorithms, bumping the version to 0.1.8.
 
## Test environments
* local OS X install, R 3.3.1
* ubuntu 14.04 (on travis-ci), R 3.3.1
* ubuntu 16.04 (locally) R 3.3.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 2 notes

* 1 NOTE relates to the install size - I would like to include a data file. 
* 1 NOTE relates to the inclusion of covr in "suggests".

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---
  
