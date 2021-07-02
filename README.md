
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/thk686/kdtools/workflows/R-CMD-check/badge.svg)](https://github.com/thk686/kdtools/actions)
[![CircleCI build
status](https://circleci.com/gh/thk686/kdtools.svg?style=svg)](https://circleci.com/ghr&svg=true)
[![Codecov test
coverage](https://codecov.io/gh/thk686/kdtools/branch/master/graph/badge.svg)](https://codecov.io/gh/thk686/kdtools?branch=master)
[![DOI](https://zenodo.org/badge/125262786.svg)](https://zenodo.org/badge/latestdoi/125262786)
[![CRAN
status](https://www.r-pkg.org/badges/version/kdtools)](https://CRAN.R-project.org/package=kdtools)

<!-- badges: end -->

# kdtools

The kdtools package exports a C++ header implementing sorting and
searching on ranges of tuple-like objects without using trees. **Note
that searching and sorting are supported on mixed-types.** It is based
on a kd-tree-like recursive sorting algorithm. Once sorted, one can
perform a range- or nearest-neighbor- query. More details are
[here](https://thk686.github.io/kdtools/). Methods and benchmarks are
[here](https://thk686.github.io/kdtools/articles/methods.html).

``` r
library(kdtools)
x = kd_sort(matrix(runif(400), 200))
plot(x, type  = 'l', asp = 1, axes = FALSE, xlab = NA, ylab = NA)
points(x, pch = 19, col = rainbow(200, alpha = 0.25), cex = 2)
y = kd_range_query(x, c(1/4, 1/4), c(3/4, 3/4))
points(y, pch = 19, cex = 0.5, col = "red")
```

<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />

## Native Data Frame Support

The core C++ header implements sorting and searching on vectors of
tuples with the number of dimensions determined at compile time. I have
generalized the package code to work on an arbitrary data frame (or any
list of equal-length vectors). This sorting and search works on any
times that are equality-comparable and less-than-comparable in the C++
STL sense.

``` r
df <- kd_sort(data.frame(a = runif(12),
                         b = as.integer(rpois(12, 1)),
                         c = sample(month.name),
                         stringsAsFactors = FALSE))
print(df)
#>             a b         c
#> 3  0.20677902 0  February
#> 6  0.15005525 1   October
#> 5  0.03177597 0 September
#> 10 0.12168593 2      July
#> 9  0.02825657 2      June
#> 2  0.22717825 2  November
#> 1  0.23413967 0  December
#> 8  0.98255966 0     April
#> 11 0.72837338 0    August
#> 12 0.61216673 0     March
#> 4  0.94480631 2   January
#> 7  0.45118667 2       May
lower <- list(0.1, 1L, "August")
upper <- list(0.9, 4L, "September")
i <- kd_rq_indices(df, lower, upper)
print(i)
#> [1]  2  4  6 12
df[i, ]
#>            a b        c
#> 6  0.1500553 1  October
#> 10 0.1216859 2     July
#> 2  0.2271783 2 November
#> 7  0.4511867 2      May
```
