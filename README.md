
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
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
    library functions for searching sorted lists (e.g., upper_bound,
    binary_search) but that operate on sequences of tuples rather than
    scalar values. It uses template metaprogramming to generate code for
    each dimension at compile time at the cost of some binary bloat.
    Sorting multidimensional arrays is fully threaded and very fast. It
    produces an implicit kd-tree split on medians that supports fast
    range queries and nearest-neighbor searches.
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
#>            a b         c
#> 9  0.3958621 0      July
#> 4  0.2823962 0       May
#> 3  0.5700275 0 September
#> 10 0.3420624 1     April
#> 7  0.2514446 1  February
#> 5  0.5407698 3      June
#> 8  0.5898558 0   October
#> 6  0.8562382 0   January
#> 11 0.8751834 0     March
#> 12 0.9086955 1    August
#> 1  0.7749514 1  December
#> 2  0.9731244 1  November
lower <- list(0.1, 1L, "August")
upper <- list(0.9, 4L, "September")
i <- kd_rq_indices(df, lower, upper)
print(i)
#> [1]  5  6 11
df[i, ]
#>           a b        c
#> 7 0.2514446 1 February
#> 5 0.5407698 3     June
#> 1 0.7749514 1 December
```
