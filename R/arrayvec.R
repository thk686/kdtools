
dim.arrayvec <- function(x) {
  return(c(x$nrow, x$ncol))
}

as.matrix.arrayvec <- function(x) {
  return(tuples_to_matrix(x))
}

as.data.frame.arrayvec <- function(x) {
  return(as.data.frame(as.matrix(x)))
}

`[.arrayvec` <- function(x, ...) {
  as.matrix(x)[...]
}

`[[.arrayvec` <- function(x, ...) {
  as.matrix(x)[[...]]
}

print.arrayvec <- function(x) {
  print(as.matrix(x))
}
