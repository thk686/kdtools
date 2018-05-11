library(kdtools)
context("Nearest neighbor")

r_nn <- function(x, y) {
  which.min(vapply(seq_len(nrow(x)),
                   function(i) { dist(rbind(x[i, ], y)) },
                   FUN.VALUE = double(1)))
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

r_nns <- function(x, y, n) {
  i = vapply(seq_len(nrow(x)),
             function(i) { dist(rbind(x[i, ], y)) },
             FUN.VALUE = double(1))
  x[which(rank(i) <= n),, drop = FALSE]
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
