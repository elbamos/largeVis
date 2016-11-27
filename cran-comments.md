## Submission 
 
This is a resubmisison of a minor version update intended as a hotfix. The response to the last submission was to inquire about test errors with version 0.1.10 and builds of R devel 3.4.

I’ve tested the current submission against R-devel 3.3.2 on linux https://travis-ci.org/elbamos/largeVis/jobs/178220156, R-devel using win builder, and 3.4 on OS X on my own machine.  I’ve tested with valgrind on my own machine using both 32-bit and 64-bit builds. 

I believe the error that was showing in the 0.1.10-r-3.4 log is related to the bug fix that is the reason for this minor version update. The error was occuring on an example that is very similar to an existing test.  I have moved the example to a test intended to replicate it exactly and encased the example version in \dontrun{}. 
 
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
 
The submission dialog asks me to confirm that CI errors have been fixed. This submission is intended to fix those errors, except that one of the errors concerns Windows with the old release of R. The log for the old-Windows build shows that the compiler in this test is forcing the use of C++ standard 0X even though `Makevars` requires C++11. I believe requiring C++11 is permissible under the CRAN standards and this is specified in the DESCRIPTION file. 
 
## Reverse dependencies 
 
There are no reverse dependencies reported on CRAN. 