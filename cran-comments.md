## Resubmission
 
This is a resubmission of a new release. 

On the last submission, a SourceTree temp file was inadvertantly left in the directory.

In this submission I have:
* Removed it.

*** I am also including the submission notes from the last submission. The reason is, it seems like UL and KH are switching-off by the day, that would make today UL's turn, and I don't want him to think that his last comments were not addressed.

On the last submission, Uwe Ligges found a crash running on Windows 32. 

In this submission, I have:
* Run the full test suite through valgrind
* Fixed a bug in the ReferenceEdges object in which, under certain circumstances, a read could take place from beyond the end of a vector. This is likely what caused Ligges' seg fault. Thank you for pointing it out! It had not shown up for any users or in any of the test environments I'd tried. 
* Reenabled two tests on 32-bit. Please note that at no point were all tests disabled on 32-bit. Kurt Hornik had seen a segfault in a single test on 32-bit through winbuilder with R-devel. Repeatedly resubmitting the same code to winbuilder did not reproduce the error. I then disabled only the individual test context on 32-bit; that test context related to an edge case where HDBSCAN would fail on certain artificial datasets with pathological features. No other tests were disabled. I noted this in cran-comments, along with links to the winbuilder logs from resubmitting the same code that KH had seen.

* NOTE:  "largeVis" is the correct name of the package and function implementing the LargeVis algorithm.

## Test environments
* local OS X install, R 3.3.1
* ubuntu 14.04 (on travis-ci)
	- R 3.3.1 & devel
	- clang & gcc
	- ARMA_64BIT_WORD & ARMA_32BIT_WORD
* OS X (on travis-ci)
	- R 3.3.1 & devel
	- Xcode 7.3, 7.4 & 8
	- ARMA_64BIT_WORD & ARMA_32BIT_WORD
* ubuntu 16.04 (locally) R 3.3.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---
  
