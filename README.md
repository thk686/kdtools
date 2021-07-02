
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/kdtools)](https://cran.r-project.org/package=kdtools)
[![Travis build
status](https://travis-ci.com/thk686/kdtools.svg?branch=master)](https://travis-ci.com/thk686/kdtools)
[![CircleCI build
status](https://circleci.com/gh/thk686/kdtools.svg?style=svg)](https://circleci.com/gh/thk686/kdtools)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/thk686/kdtools?branch=master&svg=true)](https://ci.appveyor.com/project/thk686/kdtools)
[![Coverage
status](https://codecov.io/gh/thk686/kdtools/branch/master/graph/badge.svg)](https://codecov.io/github/thk686/kdtools?branch=master)
[![DOI](https://zenodo.org/badge/125262786.svg)](https://zenodo.org/badge/latestdoi/125262786)

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
#> 8  0.11918000 0  December
#> 1  0.41809739 0      June
#> 12 0.37167230 0       May
#> 10 0.44655391 1    August
#> 9  0.09920249 2  February
#> 2  0.39759275 2   October
#> 5  0.50210280 2      July
#> 11 0.57355531 2   January
#> 7  0.57271534 2  November
#> 4  0.76871993 2     April
#> 6  0.95399264 2     March
#> 3  0.66634571 3 September
lower <- list(0.1, 1L, "August")
upper <- list(0.9, 4L, "September")
i <- kd_rq_indices(df, lower, upper)
print(i)
#> [1] 4 6 7 8 9
df[i, ]
#>            a b        c
#> 10 0.4465539 1   August
#> 2  0.3975927 2  October
#> 5  0.5021028 2     July
#> 11 0.5735553 2  January
#> 7  0.5727153 2 November
```
