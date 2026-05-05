#' Compute dot products
#'
#' @param a A numeric vector.
#' @param b A numeric vector or matrix.
#'
#' @return If `b` is a vector, the dot product of `a` and `b`.
#'    If `b` is a matrix, a vector containing the dot products of `a` and each row of `b`.
#'
#' @examples
dot <- function(a, b) {
  if (!is.vector(a)) {
    stop("a should be a vector.")
  }
  if (!is.numeric(a) | !is.numeric(b)) {
    stop("a, b should be numeric.")
  }
  if (is.vector(b)) {
    if (length(a) != length(b)) {
      stop("a, b should have the same length.")
    }
    return(sum(a*b))
  }
  b %*% a |> as.vector()
}
