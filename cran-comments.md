## Resubmission
 
This is a resubmission of a new release. 

The last submission produced this response:  
"You incorrectly have  OS_type: unix, windows  in your DESCRIPTION, which needs to be removed. Once done, we're down 
to  Possibly mis-spelled words in DESCRIPTION:  largeVis (7:29) 
Size of tarball: 8398002 bytes 
Examples with CPU or elapsed time > 5s 
user system elapsed 
largeVis 41.712 0.012 11.009 
which all need fixing. The last would seem to indicate that you make 
use of more than 2 cores, in violation of the CRAN Policy."

In this version I have:
* Removed the OS_type: line from Description (Thank you Kurt Hornik!)
* Added a parameter to all functions that use OpenMP, to limit the number of threads, and set the number of threads used in tests to 2.
* Removed the facevectors dataset to reduce size. 
* Limited the number of batches run for one example that used excessive CPU.

* NOTE:  "largeVis" is the correct name of the package and function implementing the LargeVis algorithm.
 
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
  
