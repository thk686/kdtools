library(kdtools)
context("Indexvec")

test_that("Indexvec works", {
  nr <- 1e4
  for (nc in 1:9)
  {
    x <- matrix(1:(nc * nr), nr)
    y <- matrix_to_indexed(x)
    expect_equal(dim(y), dim(x) + 0:1)
    expect_equal(ncol(y), nc + 1)
    expect_equal(nrow(y), nr)
    expect_equal(x[1, ], y[1, 1:nc])
    expect_equal(x[, 1], y[, 1])
    expect_equal(x[nr, ], y[nr, 1:nc])
    expect_equal(x[, nc], y[, nc])
    i <- sample(1:nr, 3, replace = TRUE)
    j <- sample(1:nc, 3, replace = TRUE)
    expect_equal(x[i, j], y[i, j])
    z <- indexed_to_matrix(y)
    expect_equal(cbind(x, 1:nrow(x)), z)
  }
})
