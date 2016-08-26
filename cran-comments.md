## Resubmission
 
This is a resubmission of a new release. 

The last submission produced this response:  
It failed to compile on Windows.  This was a reversion of a fix I'd put in to deal with an upstream change that happened earlier in the week.  Very sorry about that.  

In this version I have:
* Changed variable type for parameter that caused casting error on some Win32 systems. 
* Added tests and test data to resolve several rare bugs that appeared when testing on large datasets.

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
  
