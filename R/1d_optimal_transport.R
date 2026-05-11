#' One-Dimensional Wasserstein Distance Using North-West Corner Algorithm
#'
#' Applies the north-west corner algorithm to obtain the p-Wasserstein distance
#' between two one-dimensional empirical distributions.
#'
#' @param x,y Vectors representing one-dimensional empirical distributions.
#' @param presorted Set to `TRUE` if both `x` and `y` are sorted to obtain a speed-up.
#' @param p Order of the Wasserstein distance.
#' @param eps Numerical precision for floating point comparisons.
#'
#' @return The p-Wasserstein distance between `x` and `y`.
#' @export
#'
#' @examples
#' x <- rnorm(10)
#' y <- rnorm(40)
#' nw_corner_distance(x, y)
#' nw_corner_distance(sort(x), sort(y), presorted = TRUE)
nw_corner_distance <- function(x, y, presorted = FALSE, p = 2, eps = sqrt(.Machine$double.eps)) {
  nw_corner_distance_cpp(x, y, presorted = presorted, p = p, eps = eps)
}
