## Test environments
* local OS X install, R 3.6.3
* ubuntu 14.04 (on travis-ci), R 3.6.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* Old version was removed from CRAN without warning
* I assume it was because of problems with the vignettes
* I have updated them so that they operate properly with missing packages
* Also added BH to Suggests
