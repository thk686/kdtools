
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/thk686/kdtools/workflows/R-CMD-check/badge.svg)](https://github.com/thk686/kdtools/actions)
[![CircleCI build
status](https://circleci.com/gh/thk686/kdtools.svg?style=svg)](https://circleci.com/ghr&svg=true)\]
[![Coverage
status](https://codecov.io/gh/thk686/kdtools/branch/master/graph/badge.svg)](https://codecov.io/github/thk686/kdtools?branch=master)
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
#>            a b         c
#> 10 0.3791399 0     April
#> 2  0.5831830 0     March
#> 5  0.2780069 1  November
#> 1  0.5442125 1  February
#> 9  0.1364754 2    August
#> 7  0.1016411 1      July
#> 12 0.6025303 0   October
#> 6  0.8171731 0   January
#> 4  0.6560180 0 September
#> 11 0.8235863 0       May
#> 3  0.9286726 0  December
#> 8  0.6423910 1      June
lower <- list(0.1, 1L, "August")
upper <- list(0.9, 4L, "September")
i <- kd_rq_indices(df, lower, upper)
print(i)
#> [1]  3  4  5  6 12
df[i, ]
#>           a b        c
#> 5 0.2780069 1 November
#> 1 0.5442125 1 February
#> 9 0.1364754 2   August
#> 7 0.1016411 1     July
#> 8 0.6423910 1     June
```
