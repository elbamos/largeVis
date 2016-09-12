## Resubmission
 
This is a resubmission of a new release. 

On the last submission: 1) A concern was raised that tests were disabled on 32-bit Windows; and 2) It was requested to add an authorship entry for the `lof` function.

In this submission, I have:
* Re-written the `lof` function to avoid any need for attribution, and removed the license header from that file.  
* No tests are disabled on Windows 32-bit.  (Tests that require datasets that would be too large for CRAN are, however, disabled on every OS.)

* NOTE:  "largeVis" is the correct name of the package and function implementing the LargeVis algorithm.

## Test environments
* local OS X install, R 3.3.1
* ubuntu 14.04 (on travis-ci)
	- R 3.3.1 & devel
	- 32 and 64 bit Armadillo
* OS X (on travis-ci)
	- R 3.3.1 & devel
* ubuntu 16.04 (locally) R 3.3.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---
  
