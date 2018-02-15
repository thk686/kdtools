library(kdtools)
context("Tuple conversion")

test_that("can convert matrix to tuples", {
  nr = 1e4
  for (nc in 1:9)
  {
    x <- matrix(1:(nc * nr), nr)
    y <- matrix_to_tuples(x)
    expect_equal(ncol(y), nc)
    expect_equal(nrow(y), nr)
    z <- tuples_to_matrix(y)
    expect_equal(x, z)
  }
})
