get_lhs <- function(x) {
  if (!inherits(x, "formula")) stop("Formula input required")
  if (length(x) == 2) stop("Formula requires a left-hand side")
  return(all.vars(x[[2]]))
}

#' @export
knnmod <- function(formula, data, ...) {
  indep <- colspec(data, formula)
  dep <- colspec(data, get_lhs(formula))
  data <- kd_sort(data, cols = indep, ...)
  class(data) <- c("kdtools_knn", class(data))
  attr(data, "formula") <- formula
  attr(data, "indep") <- indep
  attr(data, "dep") <- dep
  return(data)
}

#' @importMethodsFrom stats predict
#' @export
predict.kdtools_knn <- function(object, newdata, k = 5, ...) {
  dep <- attr(object, "dep")
  indep <- attr(object, "indep")
  if (missing(newdata)) newdata <- object
  for (i in 1:nrow(newdata)) {
    j <- kd_nn_indices(object, newdata[i, indep], n = k,
                       cols = indep, distances = TRUE, ...)
    newdata[i, dep] <- sapply(dep, function(a) {
      weighted.mean(object[j$index, a], 1 / j$distance)
    })
  }
  return(newdata)
}

my_new_function <- function(x) return(x)

