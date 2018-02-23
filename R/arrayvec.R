#' Support for C++ vector of arrays
#' @param x an arrayvec object
#' @details Because \code{kdtools} is implemented in C++, it operates natively
#'   on a vector of arrays. An \code{arrayvec} object is a wrapper around a
#'   pointer to a vector of arrays. These functions provide some ability to
#'   manipulate the data as if it were a matrix.
#' @rdname arrayvec
#' @export
dim.arrayvec <- function(x) {
  return(c(x$nrow, x$ncol))
}

#' @rdname arrayvec
#' @export
as.matrix.arrayvec <- function(x) {
  return(tuples_to_matrix(x))
}

#' @rdname arrayvec
#' @export
as.data.frame.arrayvec <- function(x) {
  return(as.data.frame(as.matrix(x)))
}

#' @rdname arrayvec
#' @export
`[.arrayvec` <- function(x, i, j, drop = TRUE) {
  if (missing(i)) i <- 1:nrow(x); rng = range(i - 1)
  tuples_to_matrix_rows(x, rng[1], rng[2])[i - rng[1], j, drop = drop]
}

#' @rdname arrayvec
#' @export
`[[.arrayvec` <- function(x, ...) {
  as.matrix(x)[[...]]
}

#' @rdname arrayvec
#' @export
print.arrayvec <- function(x) {
  if (nrow(x) > 5)
  {
    if (ncol(x) > 5)
    {
      print(x[1:5, 1:5, FALSE])
      cat("(continues for", nrow(x) - 5, "and", ncol(x) - 5, "more rows and columns)\n")
    } else {
      print(x[1:5,, FALSE])
      cat("(continues for", nrow(x) - 5, "more rows)\n")
    }
  }
  else print(as.matrix(x))
}
