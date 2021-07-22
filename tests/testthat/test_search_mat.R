library(kdtools)
context("Search matrix")

reps <- 5
nci <- seq(1, 9, 2)

mk_ties <- function(nc) {
  x <- double()
  for (i in 1:nc)
    x <- cbind(x, sample(1:5))
  i <- sample(1:5, 100, replace = TRUE)
  return(as.matrix(x[i,]))
}

r_lower_bound <- function(x, y) {
  for (i in seq_len(nrow(x)))
    if (all(x[i, ] >= y)) return(i)
  return(as.integer(NA))
}

test_that("correct lower bound", {
  for (ignore in 1:reps)
  {
    for (n in nci)
    {
      x <- kd_sort(matrix(runif(n * 100), ncol = n))
      y <- rep(0.5, n)
      i <- kd_lower_bound(x, y)
      j <- r_lower_bound(x, y)
      expect_equal(i, j)
    }
    for (n in nci)
    {
      x <- kd_sort(mk_ties(n))
      y <- apply(x, 2, mean)
      i <- kd_lower_bound(x, y)
      j <- r_lower_bound(x, y)
      expect_equal(i, j)
    }
  }
})

r_upper_bound <- function(x, y) {
  for (i in seq_len(nrow(x)))
    if (all(x[i, ] > y)) return(i)
  return(as.integer(NA))
}

test_that("correct upper bound", {
  for (ignore in 1:reps)
  {
    for (n in nci)
    {
      x <- kd_sort(matrix(runif(n * 100), ncol = n))
      y <- rep(0.5, n)
      i <- kd_upper_bound(x, y)
      j <- r_upper_bound(x, y)
      expect_equal(i, j)
    }
  }
  for (n in nci)
  {
    x <- kd_sort(mk_ties(n))
    y <- apply(x, 2, mean)
    i <- kd_upper_bound(x, y)
    j <- r_upper_bound(x, y)
    expect_equal(i, j)
  }
})

r_contains <- function(x, a, b) {
  res <- matrix(nrow = 0, ncol = ncol(x))
  for (i in seq_len(nrow(x)))
    if (all(x[i, ] >= a) && all(x[i, ] < b)) {
      res <- rbind(res, x[i, ])
    }
  return(res)
}

test_that("range query works", {
  for (ignore in 1:reps)
  {
    for (n in nci)
    {
      x <- matrix(runif(n * 100), ncol = n)
      y <- kd_sort(x)
      l <- rep(0.25, n)
      u <- rep(0.75, n)
      z1 <- kd_range_query(y, l, u)
      z2 <- r_contains(x, l, u)
      z1 <- kd_sort(z1)
      z2 <- kd_sort(z2)
      expect_equal(nrow(z1), nrow(z2))
      if (nrow(z1) > 0) expect_equal(z1, z2)
    }
    for (n in nci)
    {
      x <- mk_ties(n)
      y <- kd_sort(x)
      l <- rep(0.25, n)
      u <- rep(0.75, n)
      z1 <- kd_range_query(y, l, u)
      z2 <- r_contains(x, l, u)
      z1 <- kd_sort(z1)
      z2 <- kd_sort(z2)
      expect_equal(nrow(z1), nrow(z2))
      if (nrow(z1) > 0) expect_equal(z1, z2)
    }
  }
})

r_contains_indices <- function(x, a, b) {
  res <- integer()
  for (i in seq_len(nrow(x)))
    if (all(x[i, ] >= a) && all(x[i, ] < b)) {
      res <- c(res, i)
    }
  return(res)
}

test_that("range query works", {
  for (ignore in 1:reps)
  {
    for (n in nci)
    {
      x <- matrix(runif(n * 100), ncol = n)
      y <- kd_sort(x)
      l <- rep(0.25, n)
      u <- rep(0.75, n)
      z1 <- kd_rq_indices(y, l, u)
      z2 <- r_contains_indices(y, l, u)
      z1 <- sort(z1)
      z2 <- sort(z2)
      expect_equal(z1, z2)
    }
    for (n in nci)
    {
      x <- mk_ties(n)
      y <- kd_sort(x)
      l <- rep(0.25, n)
      u <- rep(0.75, n)
      z1 <- kd_rq_indices(y, l, u)
      z2 <- r_contains_indices(y, l, u)
      z1 <- sort(z1)
      z2 <- sort(z2)
      expect_equal(z1, z2)
    }
  }
})

r_search <- function(x, y) {
  for (i in seq_len(nrow(x)))
    if (all(x[i, ] == y)) {
      return(TRUE)
    }
  return(FALSE)
}

test_that("binary search works", {
  for (ignore in 1:reps)
  {
    for (n in nci)
    {
      x <- kd_sort(matrix(runif(n * 100), ncol = n))
      y <- x[sample(seq_len(nrow(x)), 1), , drop = FALSE]
      expect_equal(r_search(x, y), kd_binary_search(x, y))
      expect_equal(r_search(x, rep(-1, n)), kd_binary_search(x, rep(-1, n)))
    }
  }
})
