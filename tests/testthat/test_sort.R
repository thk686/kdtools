library(kdtools)
context("Sorting")

mk_ties <- function(nc) {
  x <- double()
  for (i in 1:nc)
    x <- cbind(x, sample(1:5))
  i <- sample(1:5, 100, replace = TRUE)
  return(as.matrix(x[i,]))
}

check_median <- function(x, j = 1) {
  if (nrow(x) == 1) return(TRUE)
  i <- nrow(x) %/% 2 + 1
  if (x[i, j] != sort(x[, j])[i]) return(FALSE)
  left_ans <- ifelse(i > 1, check_median(x[1:(i - 1), , drop = FALSE], j %% ncol(x) + 1), TRUE)
  right_ans <- ifelse(i < nrow(x), check_median(x[(i + 1):nrow(x), , drop = FALSE], j %% ncol(x) + 1), TRUE)
  return(left_ans & right_ans)
}

test_that("sort works on single point", {
  for (i in c(0, 10))
    expect_error(kd_sort(matrix(i, nc = i)))
  for (i in 1:9)
    expect_equal(kd_sort(matrix(i, nc = i)), matrix(i, nc = i))
})


if (using_circular_lexicographical_compare()) {
  test_that("handles circular tie breaking", {
    x <- rnorm(10)
    expect_equal(kd_sort(cbind(0, x)), cbind(0, sort(x)))
    expect_equal(kd_sort(cbind(0, 1, x)), cbind(0, 1, sort(x)))
    expect_equal(kd_sort(cbind(0, 1, 2, x)), cbind(0, 1, 2, sort(x)))
    expect_equal(kd_sort(cbind(0, x, 1)), cbind(0, sort(x), 1))
  })
}

test_that("correct sort order", {
  nr <- 1e2
  for (nc in 1:9)
  {
    x <- matrix(runif(nr * nc), nr)
    y <- kd_sort(x)
    expect_false(kd_is_sorted(x))
    expect_true(kd_is_sorted(y))
    expect_false(check_median(x))
    expect_true(check_median(y))
  }
  for (nc in 1:9)
  {
    x <- mk_ties(nc)
    y <- kd_sort(x)
    expect_false(kd_is_sorted(x))
    expect_true(kd_is_sorted(y))
    expect_false(check_median(x))
    expect_true(check_median(y))
  }
})

test_that("correct sort order parallel", {
  nr <- 1e2
  for (nc in 1:9)
  {
    x <- matrix(runif(nr * nc), nr)
    expect_equal(kd_sort(x, parallel = TRUE), kd_sort(x, parallel = FALSE))
  }
})

kd_order_sort <- function(x) x[kd_order(x),, drop = FALSE]

test_that("correct kd_order works", {
  nr <- 1e2
  for (nc in 1:9)
  {
    x <- matrix(runif(nr * nc), nr)
    y <- kd_order_sort(x)
    expect_false(kd_is_sorted(x))
    expect_true(kd_is_sorted(y))
    expect_false(check_median(x))
    expect_true(check_median(y))
  }
  for (nc in 1:9)
  {
    x <- mk_ties(nc)
    y <- kd_order_sort(x)
    expect_false(kd_is_sorted(x))
    expect_true(kd_is_sorted(y))
    expect_false(check_median(x))
    expect_true(check_median(y))
  }
})

test_that("kd_order inplace is correct", {
  nr <- 1e2
  for (nc in 1:9) {
    x <- matrix(runif(nr * nc), nr)
    x.av <- matrix_to_tuples(x)
    i <- kd_order(x.av, inplace = TRUE)
    y <- tuples_to_matrix(x.av)
    expect_true(kd_is_sorted(x.av))
    expect_true(kd_is_sorted(y))
    expect_true(check_median(y))
    expect_equal(x[i,, drop = FALSE], y)
  }
  for (nc in 1:9) {
    x <- mk_ties(nc)
    x.av <- matrix_to_tuples(x)
    i <- kd_order(x.av, inplace = TRUE)
    y <- tuples_to_matrix(x.av)
    expect_true(kd_is_sorted(x.av))
    expect_true(kd_is_sorted(y))
    expect_true(check_median(y))
    expect_equal(x[i,, drop = FALSE], y)
  }
})




