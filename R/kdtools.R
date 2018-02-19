kd_sort <- function(x, ...) UseMethod("kd_sort")

kd_sort.matrix <- function(x, parallel = FALSE) {
  y <- matrix_to_tuples(x)
  kd_sort_(y, inplace = TRUE, parallel = parallel)
  return(tuples_to_matrix(y))
}

kd_sort.arrayvec <- function(x, inplace = FALSE, parallel = FALSE) {
  return(kd_sort_(x, inplace = inplace, parallel = parallel))
}

lex_sort <- function(x, ...) UseMethod("lex_sort")

lex_sort.matrix <- function(x) {
  y <- matrix_to_tuples(x)
  lex_sort_(y, inplace = TRUE)
  return(tuples_to_matrix(y))
}

lex_sort.arrayvec <- function(x, inplace = FALSE) {
  return(lex_sort_(x, inplace = inplace))
}

kd_lower_bound <- function(x, v) UseMethod("kd_lower_bound")

kd_lower_bound.matrix <- function(x, v) {
  y <- matrix_to_tuples(x)
  return(kd_lower_bound_(y, v))
}

kd_lower_bound.arrayvec <- function(x, v) {
  return(kd_lower_bound_(x, v))
}

kd_upper_bound <- function(x, v) UseMethod("kd_upper_bound")

kd_upper_bound.matrix <- function(x, v) {
  y <- matrix_to_tuples(x)
  return(kd_upper_bound_(y, v))
}

kd_upper_bound.arrayvec <- function(x, v) {
  return(kd_upper_bound_(x, v))
}

kd_range_query <- function(x, l, u) UseMethod("kd_range_query")

kd_range_query.matrix <- function(x, l, u) {
  y <- matrix_to_tuples(x)
  z <- kd_range_query_(y, l, u)
  return(tuples_to_matrix(z))
}

kd_range_query.arrayvec <- function(x, l, u) {
  return(kd_range_query_(x, l, u))
}

kd_nearest_neighbor <- function(x, v) UseMethod("kd_nearest_neighbor")

kd_nearest_neighbor.matrix <- function(x, v) {
  y <- matrix_to_tuples(x)
  return(kd_nearest_neighbor_(y, v))
}

kd_nearest_neighbor.arrayvec <- function(x, v) {
  return(kd_nearest_neighbor_(x, v))
}

kd_approx_nn <- function(x, v, eps) UseMethod("kd_approx_nn")

kd_approx_nn.matrix <- function(x, v, eps) {
  y <- matrix_to_tuples(x)
  return(kd_approx_nn_(y, v, eps))
}

kd_approx_nn.arrayvec <- function(x, v, eps) {
  return(kd_approx_nn_(x, v, eps))
}

kd_binary_search <- function(x, v) UseMethod("kd_binary_search")

kd_binary_search.matrix <- function(x, v) {
  y <- matrix_to_tuples(x)
  return(kd_binary_search_(y, v))
}

kd_binary_search.arrayvec <- function(x, v) {
  return(kd_binary_search_(x, v))
}

kd_nearest_neighbors <- function(x, v, n) UseMethod("kd_nearest_neighbors")

kd_nearest_neighbors.matrix <- function(x, v, n) {
  y <- matrix_to_tuples(x)
  z <- kd_nearest_neighbors_(y, v, n)
  return(tuples_to_matrix(z))
}

kd_nearest_neighbors.arrayvec <- function(x, v, n) {
  return(kd_nearest_neighbors_(x, v, n))
}
