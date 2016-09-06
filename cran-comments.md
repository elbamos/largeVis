## Resubmission
 
This is a resubmission of a new release. 

The last submission had a space between DOI and the DOI number.  
The last submission also showed an exception in a test run on winbuilder on i386 with R-devel.

I could not replicate the winbuilder error.  (See http://win-builder.r-project.org/6YruaSH7GHJV http://win-builder.r-project.org/Zy8AVTN130Kv http://win-builder.r-project.org/usCkIgX5m4F1 http://win-builder.r-project.org/AugD9S5NGEgM http://win-builder.r-project.org/NZWg67o6dKd6 http://win-builder.r-project.org/Xe0ibnk9oAz9).  I have disabled the offfending test on i386. 

In this version I have:
* Removed the space before the DOI number. 
* Disabled a test that threw an error on Windows i386.

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
  
