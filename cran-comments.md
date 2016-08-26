## Resubmission
 
This is a resubmission of a new release. 

The last submission received this response:

Can you pls provide a reference for this algorithm?

The paper is:  Tang, et al. (2016) <DOI: 10.1145/2872427.2883041> available at https://arxiv.org/abs/1602.00370

In this version I have:
* Added the DOI to the DESCRIPTION file. Note that I am not one of the paper authors. I did communicate extensively with the paper authors in preparing this implementation. They provided the datasets used in the paper, which I used to verify the correctness of the implementation. They have also reviewed and commented on the code, but not in any "official" capacity. I believe they also include this R implementation in their presentations.

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
  
