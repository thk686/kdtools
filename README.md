
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/thk686/kdtools/workflows/R-CMD-check/badge.svg)](https://github.com/thk686/kdtools/actions)
[![Codecov test
coverage](https://codecov.io/gh/thk686/kdtools/branch/master/graph/badge.svg)](https://codecov.io/gh/thk686/kdtools?branch=master)
[![DOI](https://zenodo.org/badge/125262786.svg)](https://zenodo.org/badge/latestdoi/125262786)
[![CRAN
status](https://www.r-pkg.org/badges/version/kdtools)](https://CRAN.R-project.org/package=kdtools)
[![CircleCI build
status](https://circleci.com/gh/thk686/kdtools.svg?style=svg)](https://circleci.com/gh/thk686/kdtools)

<!-- badges: end -->

# kdtools

This package emerged from a fascination with spatial partitioning
algorithms. I was curious if there was a simple and fast way to move
multidimensional data elements near to each other in memory so that they
could be retrieved quickly. I also admire the structure of the C++
standard library and was curious if the search functions for sorted
sequences could be extended to higher dimensions. At the same time, I
was interested to learn some of the newer feature of C++17.

The result is this package, which collects several ideas, including:

1.  A C++17 header file (`kdtools.h`) that can be used separately from
    this package. It has multidimensional analogs to the C++ standard
    library functions for searching sorted lists (e.g., upper\_bound,
    binary\_search) but that operate on sequences of tuples rather than
    scalar values. It uses template metaprogramming to generate code for
    each dimension at compile time. Sorting multidimensional arrays is
    fully threaded and very fast. It produces an implicit kd-tree split
    on medians that supports fast range queries and nearest-neighbor
    searches.
2.  Along the way, I needed to manipulate tuple-like objects to, for
    example, compute the euclidean distance between fixed-length arrays.
    I was curious about variadic templates and wondered if I could
    figure out nested variadic templates. The result is the included
    `tuplemapr.h` header, which allows one to apply an arbitrary
    callable over the dimensions of a collection of tuple-like objects.
    It is entirely `constexpr` (except when calling certain standard
    library functions) and employs fold-expressions where convenient.
3.  I generalized the code to work natively on data frames. This is
    slower, but still pretty fast (approx. 0.18 seconds to sort the
    nycflights13 flights data with 337k rows on an older macbook air).
    This code works on any R data type that supports comparison
    operators, including strings. Range and nearest-neighbors queries
    also support mixed types with different weights in each dimension.
4.  To support nearest-neighbor searches on strings, I coded a fast
    implementation of Levenshtein- or edit- distance. Allocating scratch
    space for the tableau is amortized constant. Most implementations
    allocate on every call, which is dramatically slower.

More details are [here](https://thk686.github.io/kdtools/). Methods and
benchmarks are
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

The following demonstrates a mixed-type range query on a data frame.

``` r
df <- kd_sort(data.frame(a = runif(12),
                         b = as.integer(rpois(12, 1)),
                         c = sample(month.name),
                         stringsAsFactors = FALSE))
print(df)
#>             a b         c
#> 11 0.07035400 0  December
#> 8  0.17601206 1  February
#> 12 0.32910504 1       May
#> 5  0.40281805 1 September
#> 10 0.05568524 2      June
#> 6  0.40731102 2     March
#> 7  0.54364982 2      July
#> 9  0.88907691 0     April
#> 4  0.57183950 0    August
#> 3  0.76644501 0   January
#> 2  0.92866102 0  November
#> 1  0.91877685 0   October
lower <- list(0.1, 1L, "August")
upper <- list(0.9, 4L, "September")
i <- kd_rq_indices(df, lower, upper)
print(i)
#> [1] 2 3 6 7
df[i, ]
#>            a b        c
#> 8  0.1760121 1 February
#> 12 0.3291050 1      May
#> 6  0.4073110 2    March
#> 7  0.5436498 2     July
```
