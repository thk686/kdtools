library(kdtools)
context("Nearest neighbor data.frame")

reps <- 3
nci <- seq(1, 9, 2)

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
#   for (ignore in 1:reps)
#   {
#     for (n in nci)
#     {
#       x <- as.data.frame(matrix(runif(n * 100), ncol = n))
#       x <- kd_sort(x)
#       y <- runif(n)
#       i <- kd_nearest_neighbor(x, y)
#       j <- r_nn(x, y)
#       expect_equal(i, j)
#     }
#     for (n in nci)
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
  for (ignore in 1:reps)
  {
    for (n in nci)
    {
      for (m in c(1, 10, 50))
      {
        x <- as.data.frame(matrix(runif(n * 100), ncol = n))
        x <- kd_sort(x)
        y <- runif(n)
        z1 <- kd_nearest_neighbors(x, y, m)
        z2 <- r_nns(x, y, m)
        expect_equal(kd_sort(z1), kd_sort(z2))
      }
    }
    for (n in nci)
    {
      for (m in c(1, 10, 50))
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
  for (ignore in 1:reps)
  {
    for (n in nci)
    {
      for (m in c(1, 10, 50))
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

r_nns_i_dist <- function(x, y, n) {
  i <- vapply(seq_len(nrow(x)),
              function(i) { dist(rbind(x[i, ], y)) },
              FUN.VALUE = double(1))
  j <- which(rank(i, ties.method = "first") <= n)
  res <- data.frame(index = j, distance = i[j])
  res <- res[order(res$distance), ]
  rownames(res) <- 1:nrow(res)
  return(res)
}

test_that("nearest neighbors distances works", {
  for (ignore in 1:reps)
  {
    for (n in nci)
    {
      for (m in c(1, 10, 50))
      {
        x <- as.data.frame(matrix(runif(n * 100), ncol = n))
        x <- kd_sort(x)
        y <- runif(n)
        z1 <- kd_nn_indices(x, y, m, distances = TRUE)
        z2 <- r_nns_i_dist(x, y, m)
        expect_equal(z1, z2)
      }
    }
  }
})

test_that("approximate nearest neighbors distances works", {
  for (ignore in 1:reps)
  {
    for (n in nci)
    {
      for (a in c(0.05, 0.5, 2))
      {
        x <- as.data.frame(matrix(runif(n * 100), ncol = n))
        x <- kd_sort(x)
        y <- runif(n)
        z1 <- kd_nn_indices(x, y, 5, a = a, distances = TRUE)
        z2 <- r_nns_i_dist(x, y, 5)
        expect_true(all(!(z1$distance > (1 + a) * z2$distance)))
      }
    }
  }
})

r_nns_i_dist_1 <- function(x, y, n, p) {
  i <- vapply(seq_len(nrow(x)),
              function(i) { dist(rbind(x[i, ], y), method = "minkowski", p = 1) },
              FUN.VALUE = double(1))
  j <- which(rank(i, ties.method = "first") <= n)
  res <- data.frame(index = j, distance = i[j])
  res <- res[order(res$distance), ]
  rownames(res) <- 1:nrow(res)
  return(res)
}

test_that("nearest neighbors manhattan distance works", {
  for (ignore in 1:reps)
  {
    for (n in nci)
    {
      for (m in c(1, 10, 50))
      {
        x <- as.data.frame(matrix(runif(n * 100), ncol = n))
        x <- kd_sort(x)
        y <- runif(n)
        z1 <- kd_nn_indices(x, y, m, p = 1, distances = TRUE)
        z2 <- r_nns_i_dist_1(x, y, m)
        expect_equal(z1, z2)
      }
    }
  }
})



