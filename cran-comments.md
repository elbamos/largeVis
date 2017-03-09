## Submission 
 
This is the a resubmission of a minor version update, to 0.2. 
This resubmission addresses issues related to the new checks for `R_registerRoutines`

I am seeing warnings with R_devel regarding "Found ‘abort’, possibly from ‘abort’ (C)." I believe the warning is a false alarm that arises out of linking the `RcppProgress` package, which has a function called `abort()`. I am working with the author of that package on an update to resolve the issue.  

The notes for the original submission were:
There are several feature and performance changes, which are detailed in NEWS.md
This update also fixes an installation error created by an update to the `dbscan` package that changed the name of a function. The dependency on `dbscan` has been removed entirely. 
 
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
 
* The Note concerns the installation size. 
 
## Note regarding test errors 
 
The submission dialog asks me to confirm that CI errors have been fixed. This submission fixes all errors related to the update to `dbscan`. This submission does not address errors related to old-Windows or Solaris. The log for the old-Windows build shows that the compiler in this test is forcing the use of C++ standard 0X even though `Makevars` requires C++11. I believe requiring C++11 is permissible under the CRAN standards and this is specified in the DESCRIPTION file. Regarding Solaris, I am unable to reproduce the error using gcc on Solaris. I believe the error is tied to Solaris Studio and to a buggy implementation of OpenMP in particular. 
 
## Reverse dependencies 
 
There are no reverse dependencies reported on CRAN. 