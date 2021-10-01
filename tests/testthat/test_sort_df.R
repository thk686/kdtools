library(kdtools)
context("Sorting")

nci <- seq(1, 9, 2)

mk_ties <- function(nc) {
  x <- double()
  for (i in 1:nc)
    x <- cbind(x, sample(1:5))
  i <- sample(1:5, 100, replace = TRUE)
  return(as.data.frame(x[i,]))
}

check_median <- function(x, j = 1) {
  if (nrow(x) == 1) return(TRUE)
  i <- nrow(x) %/% 2 + 1
  if (x[i, j] != sort(x[, j])[i]) return(FALSE)
  left_ans <- ifelse(i > 1, check_median(x[1:(i - 1), , drop = FALSE], j %% ncol(x) + 1), TRUE)
  right_ans <- ifelse(i < nrow(x), check_median(x[(i + 1):nrow(x), , drop = FALSE], j %% ncol(x) + 1), TRUE)
  return(left_ans & right_ans)
}

test_that("handles circular tie breaking", {
  x <- rnorm(10)
  cbind_df <- function(...) as.data.frame(cbind(...))
  f <- function(x) {
    row.names(x) <- NULL
    names(x) <- NULL
    return(x)
  }
  expect_equal(f(kd_sort(cbind_df(0, x))),       f(cbind_df(0, sort(x))))
  expect_equal(f(kd_sort(cbind_df(0, 1, x))),    f(cbind_df(0, 1, sort(x))))
  expect_equal(f(kd_sort(cbind_df(0, 1, 2, x))), f(cbind_df(0, 1, 2, sort(x))))
  expect_equal(f(kd_sort(cbind_df(0, x, 1))),    f(cbind_df(0, sort(x), 1)))
})

test_that("correct sort order", {
  nr <- 1e2
  for (nc in nci)
  {
    x <- as.data.frame(matrix(runif(nr * nc), nr))
    y <- kd_sort(x)
    expect_false(kd_is_sorted(x))
    expect_true(kd_is_sorted(y))
    expect_false(check_median(x))
    expect_true(check_median(y))
  }
  for (nc in nci)
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
  for (nc in nci)
  {
    x <- matrix(runif(nr * nc), nr)
    expect_equal(kd_sort(x, parallel = TRUE), kd_sort(x, parallel = FALSE))
  }
})

kd_order_sort <- function(x) x[kd_order(x),, drop = FALSE]

test_that("correct kd_order works", {
  nr <- 1e2
  for (nc in nci)
  {
    x <- matrix(runif(nr * nc), nr)
    y <- kd_order_sort(x)
    expect_false(kd_is_sorted(x))
    expect_true(kd_is_sorted(y))
    expect_false(check_median(x))
    expect_true(check_median(y))
  }
  for (nc in nci)
  {
    x <- mk_ties(nc)
    y <- kd_order_sort(x)
    expect_false(kd_is_sorted(x))
    expect_true(kd_is_sorted(y))
    expect_false(check_median(x))
    expect_true(check_median(y))
  }
})





