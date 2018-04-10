
<!-- README.md is generated from README.Rmd. Please edit that file -->
kdtools
=======

The kdtools package exports a C++ header implementing sorting and searching on ranges of tuple-like objects without using trees. It is based on a kd-tree-like recursive sorting algorithm. Once sorted, one can perform a range- or nearest-neighbor- query. More details are [here](https://thk686.github.io/kdtools/). Methods and benchmarks are [here](https://thk686.github.io/kdtools/articles/methods.html).

``` r
library(kdtools)
x = kd_sort(matrix(runif(400), 200))
plot(x, type  = 'l', asp = 1, axes = FALSE, xlab = NA, ylab = NA)
points(x, pch = 19, col = rainbow(200, alpha = 0.25), cex = 2)
y = kd_range_query(x, c(1/4, 1/4), c(3/4, 3/4))
points(y, pch = 19, cex = 0.5, col = "red")
```

<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />

[![CRAN status](https://www.r-pkg.org/badges/version/kdtools)](https://cran.r-project.org/package=kdtools) [![Travis build status](https://travis-ci.org/thk686/kdtools.svg?branch=master)](https://travis-ci.org/thk686/kdtools) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/thk686/kdtools?branch=master&svg=true)](https://ci.appveyor.com/project/thk686/kdtools) [![Coverage status](https://codecov.io/gh/thk686/kdtools/branch/master/graph/badge.svg)](https://codecov.io/github/thk686/kdtools?branch=master) [![DOI](https://zenodo.org/badge/125262786.svg)](https://zenodo.org/badge/latestdoi/125262786)
