## Changes from Prior Version

* Re-implemented DBSCAN and OPTICS. These had been taken out pending determination of whether I needed to include licensing information for the package on which I'd bsaed the code. On inspection, it turned out that the code at issue was virtually a cut-and-paste from the pseudocode for the relevant algorithms found on Wikipedia. The code has been completely re-written for this version and is substantially improved.  

* Added momentum to the largeVis stochastic gradient descent function.  In testing, this sped-up the very slow sgd phase of the algorithm by up to 10x without loss of fidelity. 
 
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

* This is a minor version update.

## Reverse dependencies

As of yet, none. 

---
  
