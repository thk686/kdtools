#' Support for C++ vector of arrays
#' @param x an indexvec object
#' @details Because \code{kdtools} is implemented in C++, it operates natively
#'   on a vector of arrays. An \code{indexvec} object is a wrapper around a
#'   pointer to a vector of arrays. These functions provide some ability to
#'   manipulate the data as if it were a matrix.
#'
#' @aliases indexvec
#' @rdname indexvec
#' @export
print.indexvec <- function(x, ...) {
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
}

#' @rdname indexvec
#' @export
dim.indexvec <- function(x) {
  return(c(x$nrow, x$ncol))
}

#' @param ... other parameters
#' @rdname indexvec
#' @export
as.matrix.indexvec <- function(x, ...) {
  return(indexed_to_matrix(x))
}

#' @rdname indexvec
#' @export
as.data.frame.indexvec <- function(x, ...) {
  return(as.data.frame(as.matrix(x)))
}

#' @param i row
#' @param j column
#' @param drop drop singleton dimensions if true
#' @rdname indexvec
#' @export
`[.indexvec` <- function(x, i, j, drop = TRUE) {
  if (missing(i)) i <- 1:nrow(x)
  rng <- range(i - 1)
  indexed_to_matrix_rows(x, rng[1], rng[2])[i - rng[1], j, drop = drop]
}

#' @rdname indexvec
#' @export
`[[.indexvec` <- function(x, ...) {
  as.matrix(x)[[...]]
}
