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

test_that("approx. nearest neighbors works", {
  for (ignore in 1:10)
  {
    for (n in 1:9)
    {
      x <- matrix(runif(n * 100), nc = n)
      x <- kd_sort(x)
      y <- runif(n)
      eps <- rlnorm(1, -2)
      i <- kd_approx_nn(x, y, eps)
      j <- r_nn(x, y)
      d <- rep(0, 2)
      for (k in 1:n)
      {
        d[1] <- d[1] + (x[j, k] - y[k]) ^ 2
        d[2] <- d[2] + (x[i, k] - y[k]) ^ 2
      }
      d <- sqrt(d)
      expect_true(diff(d) < eps)
    }
  }
})

r_nns <- function(x, y, n) {
  x[which(rank(sapply(1:nrow(x), function(i) dist(rbind(x[i, ], y)))) <= n),, drop = FALSE]
}

test_that("nearest neighbors works", {
  for (ignore in 1:10)
  {
    for (n in 1:9)
    {
      for (m in c(1, 10, 2 * n * 100))
      {
        x <- matrix(runif(n * 100), nc = n)
        x <- kd_sort(x)
        y <- runif(n)
        z1 <- kd_nearest_neighbors(x, y, m)
        z2 <- r_nns(x, y, m)
        expect_equal(kd_sort(z1), kd_sort(z2))
      }
    }
  }
})
