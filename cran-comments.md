## Submission 
 
This is the first submission of a bug-fix update, to 0.2.1.
 
## Test environments 
* local OS X install, R 3.4
* OS X (with Valgrind), R 3.4
* ubuntu 14.04 (on travis-ci), R 3.4 and R-devel 
* OS X (on travis-ci), R-devel 
* win-builder (devel and release) 
 
## R CMD check results 
 
0 errors | 0 warnings | 2 notes 
 
* The Notes concern installation size, and a false-alarm regarding misspelled word.
 
## Note regarding test errors 
 
The submission dialog asks me to confirm that CI errors have been fixed. This submission does not address errors related to Solaris.  The Solaris errors appear to be related to Solaris Studio, and I have not been able to reproduce it using gcc on Solaris.  
 
## Reverse dependencies 
 
There are no reverse dependencies reported on CRAN. 