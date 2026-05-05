#' Compute dot products
#'
#' Computes the dot product between a vector and one or several other vectors.
#'
#' @param a A numeric vector.
#' @param b A numeric vector or matrix.
#'
#' @return If `b` is a vector, the dot product of `a` and `b`.
#'    If `b` is a matrix, a vector containing the dot products of `a` and each row of `b`.
#'
#' @examples
#' x <- c(1, 2, 3)
#' y <- c(2, 4, 2)
#' dot(x, y)
#' m <- rbind(x, y)
#' dot(x, m)
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
#' Compute squared norms
#'
#' @param a A numeric vector or matrix.
#'
#' @return The squared Euclidean norm of `a` (if `a` is a vector)
#'    or the squared Euclidean norms of each row of `a` (if `a` is a matrix).
#'
#' @examples
#' x <- c(1, 2, 3)
#' squared_norm(x)
#' m <- matrix(rnorm(40), nrow = 10)
#' squared_norm(m)
squared_norm <- function(a) {
  if (is.vector(a)) return(dot(a, a))
  apply(a, MARGIN = 1, squared_norm)
}
