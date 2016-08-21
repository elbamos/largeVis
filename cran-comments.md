## Resubmission
 
This is a resubmission of a new release. 

The last submission produced this response:  
    "FOSS licence with BuildVignettes: false
    Package has a VignetteBuilder field but no prebuilt vignette index."


In this version I have:
* Changed BuildVignettes: to TRUE
* Rebuilt vignettes
* Fixed a bug in the routine for plotting hdbscan results
 
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
  
