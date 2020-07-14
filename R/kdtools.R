#' @importFrom Rcpp evalCpp
#' @useDynLib kdtools, .registration = TRUE
NULL

colspec <- function(x, cols) {
  switch(mode(cols),
         "character" = {
           if (!all(cols %in% names(x)))
             warning("Non-existent column name ignored")
           which(names(x) %in% cols)
          },
         "numeric" = {
           (1:ncol(x))[cols]
         },
         "logical" = {
           (1:ncol(x))[cols]
          },
         stop("Invalid column spec"))
}

#' Sort multidimensional data
#' @param x a matrix or arrayvec object
#' @param ... ignored
#' @details The algorithm used is a divide-and-conquer quicksort variant that
#'   recursively partions an range of tuples using the median of each successive
#'   dimension. Ties are resolved by cycling over successive dimensions. The
#'   result is an ordering of tuples matching their order if they were inserted
#'   into a kd-tree.
#'
#'   \code{kd_order} returns permutation vector that will order the rows of the
#'   original matrix, exactly as \code{\link{order}}. If \code{inplace} is true,
#'   then \code{kd_order} will also sort the arrayvec object as a side effect.
#'   This can be more efficient when many subsequent queries are required.
#'
#'   \code{kd_sort} and \code{kd_order} have been extended to work directly on a
#'   data frame. All vector column types are supported (even lists of objects as
#'   long as equality and comparison operators are defined). Additionaly, the
#'   user can specify a sequence of column indices that will be used for
#'   sorting. These can be a subset of columns and given in any order.
#' @note The matrix version will be slower because of data structure
#'   conversions.
#' @examples
#' z <- data.frame(real = runif(10), lgl = runif(10) > 0.5,
#'                 int = as.integer(rpois(10, 2)), char = sample(month.name, 10),
#'                 stringsAsFactors = FALSE)
#' kd_sort(z)
#' x <- matrix(runif(200), 100)
#' y <- kd_sort(x)
#' kd_is_sorted(y)
#' kd_order(x)
#' plot(y, type = "o", pch = 19, col = "steelblue", asp = 1)
#' @seealso \code{\link{arrayvec}}
#' @rdname kdsort
#' @export
kd_sort <- function(x, ...) UseMethod("kd_sort")

#' @param parallel use multiple threads if true
#' @rdname kdsort
#' @export
kd_sort.matrix <- function(x, parallel = TRUE, ...) {
  y <- matrix_to_tuples(x)
  kd_sort_(y, inplace = TRUE, parallel = parallel)
  return(tuples_to_matrix(y))
}

#' @param inplace sort as a side-effect if true
#' @rdname kdsort
#' @export
kd_sort.arrayvec <- function(x, inplace = FALSE, parallel = TRUE, ...) {
  return(kd_sort_(x, inplace = inplace, parallel = parallel))
}

#' @rdname kdsort
#' @export
kd_sort.data.frame <- function(x, cols = 1:ncol(x), parallel = TRUE, ...) {
  return(x[kd_order(x, cols = colspec(x, cols), parallel = parallel),, drop = FALSE])
}

#' @rdname kdsort
#' @export
kd_order <- function(x, ...) UseMethod("kd_order")

#' @rdname kdsort
#' @export
kd_order.matrix <- function(x, parallel = TRUE, ...) {
  y <- matrix_to_tuples(x)
  return(kd_order_(y, inplace = FALSE, parallel = parallel))
}

#' @rdname kdsort
#' @export
kd_order.arrayvec <- function(x, inplace = FALSE, parallel = TRUE, ...) {
  return(kd_order_(x, inplace = inplace, parallel = parallel))
}

#' @param cols integer vector of column indices
#' @rdname kdsort
#' @export
kd_order.data.frame <- function(x, cols = 1:ncol(x), parallel = TRUE, ...) {
  return(kd_order_df(x, colspec(x, cols), parallel = parallel))
}

#' @rdname kdsort
#' @export
kd_is_sorted <- function(x, ...) UseMethod("kd_is_sorted")

#' @export
kd_is_sorted.matrix <- function(x, parallel = FALSE, ...) {
  return(kd_is_sorted_(matrix_to_tuples(x), parallel))
}

#' @export
kd_is_sorted.arrayvec <- function(x, parallel = FALSE, ...) {
  return(kd_is_sorted_(x, parallel))
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
#' @param ... addtional arguments
#' @examples
#' x = matrix(runif(200), 100)
#' y = matrix_to_tuples(x)
#' kd_sort(y, inplace = TRUE)
#' y[kd_lower_bound(y, c(1/2, 1/2)),]
#' y[kd_upper_bound(y, c(1/2, 1/2)),]
#' kd_binary_search(y, c(1/2, 1/2))
#' kd_range_query(y, c(1/3, 1/3), c(2/3, 2/3))
#' kd_rq_indices(y, c(1/3, 1/3), c(2/3, 2/3))
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
kd_range_query <- function(x, l, u, ...) UseMethod("kd_range_query")

#' @rdname search
#' @export
kd_range_query.matrix <- function(x, l, u, ...) {
  y <- matrix_to_tuples(x)
  z <- kd_range_query_(y, l, u)
  return(tuples_to_matrix(z))
}

#' @rdname search
#' @export
kd_range_query.arrayvec <- function(x, l, u, ...) {
  return(kd_range_query_(x, l, u))
}

#' @rdname search
#' @export
kd_range_query.data.frame <- function(x, l, u, cols = 1:ncol(x), ...) {
  return(x[kd_rq_indices(x, l, u, colspec(x, cols)),, drop = FALSE])
}

#' @rdname search
#' @export
kd_rq_indices <- function(x, l, u, ...) UseMethod("kd_rq_indices")

#' @rdname search
#' @export
kd_rq_indices.matrix <- function(x, l, u, ...) {
  y <- matrix_to_tuples(x)
  return(kd_rq_indices_(y, l, u))
}

#' @rdname search
#' @export
kd_rq_indices.arrayvec <- function(x, l, u, ...) {
  return(kd_rq_indices_(x, l, u))
}

#' @rdname search
#' @param cols integer vector of column indices
#' @export
kd_rq_indices.data.frame <- function(x, l, u, cols = 1:ncol(x), ...) {
  return(kd_rq_df(x, colspec(x, cols), l, u))
}

#' @rdname search
#' @export
kd_binary_search <- function(x, v) UseMethod("kd_binary_search")

#' @rdname search
#' @export
kd_binary_search.matrix <- function(x, v) {
  y <- matrix_to_tuples(x)
  return(kd_binary_search_(y, v))
}

#' @rdname search
#' @export
kd_binary_search.arrayvec <- function(x, v) {
  return(kd_binary_search_(x, v))
}

#' Find nearest neighbors
#' @param x an object sorted by \code{\link{kd_sort}}
#' @param v a vector specifying where to look
#' @param n the number of neighbors to return
#' @param ... addtional arguments
#'
#' @examples
#' x = matrix(runif(200), 100)
#' y = matrix_to_tuples(x)
#' kd_sort(y, inplace = TRUE)
#' y[kd_nearest_neighbor(y, c(1/2, 1/2)),]
#' kd_nearest_neighbors(y, c(1/2, 1/2), 3)
#' y[kd_nn_indices(y, c(1/2, 1/2), 5),]
#'
#' @rdname nneighb
#' @export
kd_nearest_neighbors <- function(x, v, n, ...) UseMethod("kd_nearest_neighbors")

#' @rdname nneighb
#' @export
kd_nearest_neighbors.matrix <- function(x, v, n, ...) {
  y <- matrix_to_tuples(x)
  z <- kd_nearest_neighbors_(y, v, n)
  return(tuples_to_matrix(z))
}

#' @rdname nneighb
#' @export
kd_nearest_neighbors.arrayvec <- function(x, v, n, ...) {
  return(kd_nearest_neighbors_(x, v, n))
}

#' @param cols integer indices of columns to use
#' @param w distance weights
#' @rdname nneighb
#' @export
kd_nearest_neighbors.data.frame <- function(x, v, n, cols = 1:ncol(x),
                                            w = rep(1, length(cols)), ...) {
  return(x[kd_nn_indices(x, v, n, colspec(x, cols), w),, drop = FALSE])
}

#' @rdname nneighb
#' @export
kd_nn_indices <- function(x, v, n, ...) UseMethod("kd_nn_indices")

#' @rdname nneighb
#' @export
kd_nn_indices.matrix <- function(x, v, n, ...) {
  y <- matrix_to_tuples(x)
  return(kd_nn_indices_(y, v, n))
}

#' @rdname nneighb
#' @export
kd_nn_indices.arrayvec <- function(x, v, n, ...) {
  return(kd_nn_indices_(x, v, n))
}

#' @rdname nneighb
#' @export
kd_nn_indices.data.frame <- function(x, v, n, cols = 1:ncol(x),
                                     w = rep(1, length(cols)), ...) {
  return(kd_nn_df(x, colspec(x, cols), w, v, n))
}

#' @rdname nneighb
#' @export
kd_nearest_neighbor <- function(x, v) UseMethod("kd_nearest_neighbor")

#' @rdname nneighb
#' @export
kd_nearest_neighbor.matrix <- function(x, v) {
  y <- matrix_to_tuples(x)
  return(kd_nearest_neighbor_(y, v))
}

#' @rdname nneighb
#' @export
kd_nearest_neighbor.arrayvec <- function(x, v) {
  return(kd_nearest_neighbor_(x, v))
}
