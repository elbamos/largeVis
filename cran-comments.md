## Submission 
 
This is a fix to an issue that appeared in 0.1.10 after submission. In particular, if largeVis was compiled to use 32-bit Armadillo words, with certain compilers including the most recent Apple compiler, the parameter specifying the number of training batches was not being passed from R to C++ properly. This caused the tests to time-out and checks to fail on those systems. (At the same time, this release fixes a reported bug where the neighbor search could fail with sparse matrices under certain conditions.)
 
## Test environments 
* local OS X install, R 3.3.1 
*  OS X (with Valgrind), R 3.3.1 
* ubuntu 12.04 (on travis-ci), R 3.3.1 
* ubuntu 14.04 (on travis-ci), R 3.3.1 and R-devel 
* Solaris 11 x86 (via Virtual Box), R 3.3.1 
* OS X (on travis-ci), R-devel 
* win-builder (devel and release) 
 
## R CMD check results 
 
0 errors | 0 warnings | 1 note 
 
* The Note concerns the installation size, which I've been able to reduce since the prior version. 
 
## Note regarding test errors 
 
The submission dialog asks me to confirm that CI errors have been fixed. This submission is intended to fix those errors, except that one of the errors concerns Windows with the old release of R. The log shows that the compiler in this test is forcing the use of C++ standard 0X even though `Makevars` requires C++11. I believe requiring C++11 is permissible under the CRAN standards and this is specified in the DESCRIPTION file. 
 
## Reverse dependencies 
 
There are no reverse dependencies reported on CRAN. 