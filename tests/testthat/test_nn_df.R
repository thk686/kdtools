library(kdtools)
context("Nearest neighbor data.frame")

r_nn <- function(x, y) {
  which.min(vapply(seq_len(nrow(x)),
                   function(i) { dist(rbind(x[i, ], y)) },
                   FUN.VALUE = double(1)))
}

mk_ties <- function(nc) {
  x <- double()
  for (i in 1:nc)
    x <- cbind(x, sample(1:5))
  i <- sample(1:5, 100, replace = TRUE)
  return(as.data.frame(x[i,]))
}

pair_dist <- function(a, b) sqrt(sum((a - b)^2))

# test_that("nearest neighbor works", {
#   for (ignore in 1:10)
#   {
#     for (n in 1:9)
#     {
#       x <- as.data.frame(matrix(runif(n * 100), ncol = n))
#       x <- kd_sort(x)
#       y <- runif(n)
#       i <- kd_nearest_neighbor(x, y)
#       j <- r_nn(x, y)
#       expect_equal(i, j)
#     }
#     for (n in 1:9)
#     {
#       x <- mk_ties(n)
#       x <- kd_sort(x)
#       y <- runif(n, 1, 5)
#       i <- kd_nearest_neighbor(x, y)
#       j <- r_nn(x, y)
#       expect_equal(pair_dist(y, x[i,]),
#                    pair_dist(y, x[j,]))
#     }
#   }
# })

r_nns <- function(x, y, n) {
  i = vapply(seq_len(nrow(x)),
             function(i) { dist(rbind(x[i, ], y)) },
             FUN.VALUE = double(1))
  x[which(rank(i, ties.method = "first") <= n),, drop = FALSE]
}

test_that("nearest neighbors works", {
  for (ignore in 1:10)
  {
    for (n in 1:9)
    {
      for (m in c(1, 10, 2 * n * 100))
      {
        x <- as.data.frame(matrix(runif(n * 100), ncol = n))
        x <- kd_sort(x)
        y <- runif(n)
        z1 <- kd_nearest_neighbors(x, y, m)
        z2 <- r_nns(x, y, m)
        expect_equal(kd_sort(z1), kd_sort(z2))
      }
    }
    for (n in 1:9)
    {
      for (m in c(1, 10, 2 * n * 100))
      {
        x <- mk_ties(n)
        x <- kd_sort(x)
        y <- runif(n)
        z1 <- kd_sort(kd_nearest_neighbors(x, y, m))
        z2 <- kd_sort(r_nns(x, y, m))
        row.names(z1) <- NULL
        row.names(z2) <- NULL
        expect_equal(z1, z2)
      }
    }
  }
})

r_nns_i <- function(x, y, n) {
  i = vapply(seq_len(nrow(x)),
             function(i) { dist(rbind(x[i, ], y)) },
             FUN.VALUE = double(1))
  which(rank(i, ties.method = "first") <= n)
}

test_that("nearest neighbors indices works", {
  for (ignore in 1:10)
  {
    for (n in 1:9)
    {
      for (m in c(1, 10, 2 * n * 100))
      {
        x <- as.data.frame(matrix(runif(n * 100), ncol = n))
        x <- kd_sort(x)
        y <- runif(n)
        z1 <- kd_nn_indices(x, y, m)
        z2 <- r_nns_i(x, y, m)
        expect_equal(sort(z1), sort(z2))
      }
    }
  }
})


