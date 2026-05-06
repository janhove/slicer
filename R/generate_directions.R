#' Generate random projection directions
#'
#' Obtain `L` random vectors on the unit sphere in `d`-dimensional Euclidean space.
#'
#' @param L Number of projection directions.
#' @param d Dimension.
#'
#' @return An `L` by `d` matrix, each row of which contains a unit vector.
#' @export
#'
#' @examples
#' generate_directions(10, 4)
#' plot(generate_directions(25, 2))
generate_directions <- function(L, d) {
  z <- MASS::mvrnorm(L, mu = rep(0, d), Sigma = diag(1, d))
  normalize(z)
}
