
<!-- README.md is generated from README.Rmd. Please edit that file -->

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
#> 5  0.07985814 1    August
#> 12 0.23057540 1  February
#> 11 0.02447241 0  November
#> 7  0.15839982 1      July
#> 2  0.33232334 1       May
#> 1  0.06414051 1   October
#> 3  0.36617135 2     March
#> 6  0.70503899 0     April
#> 8  0.55539120 0  December
#> 9  0.73120968 1   January
#> 4  0.99288002 3      June
#> 10 0.88029940 1 September
lower <- list(0.1, 1L, "August")
upper <- list(0.9, 4L, "September")
i <- kd_rq_indices(df, lower, upper)
print(i)
#> [1]  2  4  5  7 10
df[i, ]
#>            a b        c
#> 12 0.2305754 1 February
#> 7  0.1583998 1     July
#> 2  0.3323233 1      May
#> 3  0.3661713 2    March
#> 9  0.7312097 1  January
```

**Note: Switching to C++17 broke CI (I’d be grateful for a pull-request
fixing this as I cannot make heads-nor-tails of the required
configuration). The code builds and checks fine on my systems.**

[![CRAN
status](https://www.r-pkg.org/badges/version/kdtools)](https://cran.r-project.org/package=kdtools)
[![CircleCI build
status](https://circleci.com/gh/thk686/kdtools.svg?style=svg)](https://circleci.com/gh/thk686/kdtools)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/thk686/kdtools?branch=master&svg=true)](https://ci.appveyor.com/project/thk686/kdtools)
[![Coverage
status](https://codecov.io/gh/thk686/kdtools/branch/master/graph/badge.svg)](https://codecov.io/github/thk686/kdtools?branch=master)
[![DOI](https://zenodo.org/badge/125262786.svg)](https://zenodo.org/badge/latestdoi/125262786)
