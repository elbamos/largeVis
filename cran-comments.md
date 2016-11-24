## Submission 
 
This is a resubmisison. The last submission was slightly larger than the 5MB limit on tar file sizes.  This submission eliminates a large data file that was used to test the clustering algorithms. 
 
## Test environments 
* local OS X install, R 3.3.2 
* OS X (with Valgrind), R 3.3.2 
* ubuntu 12.04 (on travis-ci), R 3.3.2 
* ubuntu 14.04 (on travis-ci), R 3.3.2 and R-devel 
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