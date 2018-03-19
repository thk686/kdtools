#' @importFrom Rcpp evalCpp
#' @useDynLib kdtools, .registration = TRUE
NULL

#' Sort multidimensional data
#' @param x a matrix or arrayvec object
#' @param ... other arguments
#' @details The algorithm used is a divide-and-conquer quicksort variant that
#'   recursively partions an range of tuples using the median of each successive
#'   dimension. Ties are resolved by cycling over successive dimensions. The
#'   result is an ordering of tuples matching their order if they were inserted
#'   into a kd-tree.
#' @note The matrix version will be slower because of data structure
#'   conversions.
#' @examples
#' x = kd_sort(matrix(runif(200), 100))
#' kd_is_sorted(x)
#' plot(x, type = "o", pch = 19, col = "steelblue", asp = 1)
#'
#' @seealso \code{\link{arrayvec}}
#' @rdname kdsort
#' @export
kd_sort <- function(x, ...) UseMethod("kd_sort")

#' @export
kd_sort.matrix <- function(x, parallel = FALSE, ...) {
  y <- matrix_to_tuples(x)
  kd_sort_(y, inplace = TRUE, parallel = parallel)
  return(tuples_to_matrix(y))
}

#' @export
kd_sort.arrayvec <- function(x, inplace = FALSE, parallel = FALSE, ...) {
  return(kd_sort_(x, inplace = inplace, parallel = parallel))
}

#' @rdname kdsort
#' @export
kd_is_sorted <- function(x) UseMethod("kd_is_sorted")

#' @export
kd_is_sorted.matrix <- function(x) {
  return(kd_is_sorted_(matrix_to_tuples(x)))
}

#' @export
kd_is_sorted.arrayvec <- function(x) {
  return(kd_is_sorted_(x))
}

#' Sort a matrix into lexicographical order
#' @param x a matrix or arrayvec object
#' @param ... other parameters
#' @details Sorts a range of tuples into lexicographical order.
#' @examples
#' x = lex_sort(matrix(runif(200), 100))
#' plot(x, type = "o", pch = 19, col = "steelblue", asp = 1)
#'
#' @rdname lexsort
#' @export
lex_sort <- function(x, ...) UseMethod("lex_sort")

#' @export
lex_sort.matrix <- function(x, ...) {
  y <- matrix_to_tuples(x)
  lex_sort_(y, inplace = TRUE)
  return(tuples_to_matrix(y))
}

#' @export
lex_sort.arrayvec <- function(x, inplace = FALSE, ...) {
  return(lex_sort_(x, inplace = inplace))
}

#' Search sorted data
#' @param x an object sorted by \code{\link{kd_sort}}
#' @param v a vector specifying where to look
#' @examples
#' x = matrix(runif(200), 100)
#' y = matrix_to_tuples(x)
#' kd_sort(y, inplace = TRUE)
#' y[kd_lower_bound(y, c(1/2, 1/2)),]
#' y[kd_upper_bound(y, c(1/2, 1/2)),]
#' kd_binary_search(y, c(1/2, 1/2))
#' kd_range_query(y, c(1/3, 1/3), c(2/3, 2/3))
#'
#' @aliases kd_lower_bound
#' @rdname search
#' @export
kd_lower_bound <- function(x, v) UseMethod("kd_lower_bound")

#' @export
kd_lower_bound.matrix <- function(x, v) {
  y <- matrix_to_tuples(x)
  return(kd_lower_bound_(y, v))
}

#' @export
kd_lower_bound.arrayvec <- function(x, v) {
  return(kd_lower_bound_(x, v))
}

#' @rdname search
#' @export
kd_upper_bound <- function(x, v) UseMethod("kd_upper_bound")

#' @export
kd_upper_bound.matrix <- function(x, v) {
  y <- matrix_to_tuples(x)
  return(kd_upper_bound_(y, v))
}

#' @export
kd_upper_bound.arrayvec <- function(x, v) {
  return(kd_upper_bound_(x, v))
}

#' @param l lower left corner of search region
#' @param u upper right corner of search region
#' @rdname search
#' @export
kd_range_query <- function(x, l, u) UseMethod("kd_range_query")

#' @export
kd_range_query.matrix <- function(x, l, u) {
  y <- matrix_to_tuples(x)
  z <- kd_range_query_(y, l, u)
  return(tuples_to_matrix(z))
}

#' @export
kd_range_query.arrayvec <- function(x, l, u) {
  return(kd_range_query_(x, l, u))
}

#' @rdname search
#' @export
kd_binary_search <- function(x, v) UseMethod("kd_binary_search")

#' @export
kd_binary_search.matrix <- function(x, v) {
  y <- matrix_to_tuples(x)
  return(kd_binary_search_(y, v))
}

#' @export
kd_binary_search.arrayvec <- function(x, v) {
  return(kd_binary_search_(x, v))
}

#' Find nearest neighbors
#' @param x an object sorted by \code{\link{kd_sort}}
#' @param v a vector specifying where to look
#' @param n the number of neighbors to return
#'
#' @examples
#' x = matrix(runif(200), 100)
#' y = matrix_to_tuples(x)
#' kd_sort(y, inplace = TRUE)
#' y[kd_nearest_neighbor(y, c(1/2, 1/2)),]
#' kd_nearest_neighbors(y, c(1/2, 1/2), 3)
#'
#' @rdname nneighb
#' @export
kd_nearest_neighbors <- function(x, v, n) UseMethod("kd_nearest_neighbors")

#' @export
kd_nearest_neighbors.matrix <- function(x, v, n) {
  y <- matrix_to_tuples(x)
  z <- kd_nearest_neighbors_(y, v, n)
  return(tuples_to_matrix(z))
}

#' @export
kd_nearest_neighbors.arrayvec <- function(x, v, n) {
  return(kd_nearest_neighbors_(x, v, n))
}

#' @rdname nneighb
#' @export
kd_nearest_neighbor <- function(x, v) UseMethod("kd_nearest_neighbor")

#' @export
kd_nearest_neighbor.matrix <- function(x, v) {
  y <- matrix_to_tuples(x)
  return(kd_nearest_neighbor_(y, v))
}

#' @export
kd_nearest_neighbor.arrayvec <- function(x, v) {
  return(kd_nearest_neighbor_(x, v))
}
