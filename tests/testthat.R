library(testthat)
library(kdtools)

if (has_cxx17()) test_check("kdtools")
