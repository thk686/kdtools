#' Support for C++ vector of arrays
#' @param x an arrayvec object
#' @details Because \code{kdtools} is implemented in C++, it operates natively
#'   on a vector of arrays. An \code{arrayvec} object is a wrapper around a
#'   pointer to a vector of arrays. These functions provide some ability to
#'   manipulate the data as if it were a matrix.
#'
#' @return \tabular{ll}{
#' \code{print.arrayvec} \tab the object invisibly \cr
#' \code{dim.arrayvec} \tab the rows and columns \cr
#' \code{as.matrix.arrayvec} \tab a matrix \cr
#' \code{as.data.frame.arrayvec} \tab a data frame \cr
#' \code{`[.arrayvec`} \tab a matrix or vector \cr
#' \code{`[[.arrayvec`} \tab a column vector \cr
#' }
#' @aliases arrayvec
#' @rdname arrayvec
#' @export
print.arrayvec <- function(x, ...) {
  if (nrow(x) > 5) {
    if (ncol(x) > 5) {
      print(x[1:5, 1:5, FALSE])
      cat("(continues for", nrow(x) - 5, "and", ncol(x) - 5, "more rows and columns)\n")
    } else {
      print(x[1:5, , FALSE])
      cat("(continues for", nrow(x) - 5, "more rows)\n")
    }
  }
  else {
    print(as.matrix(x))
  }
  return(invisible(x))
}

#' @rdname arrayvec
#' @export
dim.arrayvec <- function(x) {
  return(c(x$nrow, x$ncol))
}

#' @param ... other parameters
#' @rdname arrayvec
#' @export
as.matrix.arrayvec <- function(x, ...) {
  return(tuples_to_matrix(x))
}

#' @rdname arrayvec
#' @export
as.data.frame.arrayvec <- function(x, ...) {
  return(as.data.frame(as.matrix(x)))
}

#' @param i row
#' @param j column
#' @param drop drop singleton dimensions if true
#' @rdname arrayvec
#' @export
# `[.arrayvec` <- function(x, i, j, drop = TRUE) {
#   if (missing(i)) i <- seq_len(nrow(x))
#   rng <- range(i)
#   tuples_to_matrix_rows(x, rng[1], rng[2])[i - rng[1] + 1, j, drop = drop]
# }
`[.arrayvec` <- function(x, ...) {
  as.matrix(x)[...]
}


#' @rdname arrayvec
#' @export
`[[.arrayvec` <- function(x, ...) {
  as.matrix(x)[[...]]
}
