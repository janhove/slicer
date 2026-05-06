#' Compute dot products
#'
#' Computes the dot product between a vector and one or several other vectors.
#'
#' @noRd
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
#' Computes the squared Euclidean norm(s) of one or more vectors.
#'
#' @noRd
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
#' Compute norms
#'
#' Computes the Euclidean norm(s) of one or more vectors.
#'
#' @noRd
#'
#' @param a A numeric vector or matrix.
#'
#' @return The Euclidean norm of `a` (if `a` is a vector)
#'    or the Euclidean norms of each row of `a` (if `a` is a matrix).
#'
#' @examples
#' norm(c(1, 2, 3))
#' norm(matrix(rnorm(90), nrow = 9))
norm <- function(a) {
  sqrt(squared_norm(a))
}
#' Normalize vectors to unit length
#'
#' @noRd
#'
#' @param a A numeric vector or matrix.
#'
#' @return `a` scaled to have unit norm (if `a` is a vector) or
#'    `a` with each row scaled to have unit length (if `a` is a matrix).
#'
#' @examples
#' x <- c(1, 2, 3)
#' normalize(x)
#' m <- matrix(rnorm(40), nrow = 10)
normalize <- function(a) {
  if (is.vector(a)) {
    if (isTRUE(all.equal(norm(a), 0))) return(NaN)
    return(a / norm(a))
  }
  apply(a, MARGIN = 1, normalize) |> t()
}
#' Find median value in upper triangle
#'
#' @noRd
#'
#' @param M A square matrix.
#'
#' @return The median value in the upper triangle.
#'
#' @examples
#' m <- matrix(rnorm(100), nrow = 10)
#' find_median(m)
#' find_median(m + diag(100, 10))
find_median <- function(M) {
  stats::median(M[upper.tri(M)])
}
#' Sorted projections of a data matrix.
#'
#' @noRd
#'
#' @param M n by d matrix representing an empirical distribution.
#' @param thetas L by d matrix representing projection directions.
#'
#' @return An n by L matrix, the i-th column of which is the sorted
#'         empirical distribution of the projection of M onto the i-th
#'         row of thetas.
#'
#' @examples
#' my_data <- matrix(rnorm(100), ncol = 5)
#' directions <- generate_directions(50, 5)
#' project_and_sort(my_data, directions)
project_and_sort <- function(M, thetas) {
  thetas <- t(thetas)
  apply(M %*% thetas, 2, sort)
}
