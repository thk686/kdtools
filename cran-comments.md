## Test environments
* local R installation, R 4.1.0
* ubuntu 16.04 (on travis-ci), R 4.1.0
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Release notes

* I have added a configure/configure.win so that the package now should pass checks when CXX17 is not defined
* I have tested this on rhub and winbuilder -- hopefully I got it right this time!

