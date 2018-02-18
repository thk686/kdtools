library(kdtools)
context("Nearest neighbor")

r_nn <- function(x, y) {
  which.min(sapply(1:nrow(x), function(i) dist(rbind(x[i, ], y))))
}

test_that("nearest neighbors works", {
  for (ignore in 1:10)
  {
    for (n in 1:9)
    {
      x <- matrix(runif(n * 100), nc = n)
      x <- kd_sort(x)
      y <- runif(n)
      i <- kd_nearest_neighbor(x, y)
      j <- r_nn(x, y)
      expect_equal(i, j)
    }
  }
})
