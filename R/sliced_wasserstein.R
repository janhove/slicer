#' Sliced Wasserstein Distance Between Two Empirical Distributions
#'
#' @description
#' Estimate the sliced Wasserstein distance between two empirical distributions
#' represented as matrices using Monte Carlo simulation.
#'
#' Optionally, a matrix `thetas` can be supplied. Each row of `thetas` will
#' be interpreted as a projection direction. These rows do not have to be
#' unit vectors, allowing users to transform the projection directions
#' rather than the matrices.
#'
#' @param x,y Matrices representing empirical distributions.
#' @param p Order of sliced Wasserstein distance.
#' @param thetas Optionally, a matrix, each row of which represents a projection direction.
#' @param L If no `thetas` are provided, `L` random projection directions are generated.
#' @param seed Optional random seed.
#'
#' @return Estimated sliced Wasserstein distance of order `p` between `x` and `y`.
#' @export
#'
#' @examples
#' M1 <- matrix(rnorm(50), ncol = 5)
#' M2 <- matrix(rnorm(250), ncol = 5)
#' sliced_wasserstein(M1, M2) # random projection directions
#' cardinal_axes <- diag(1, 5)
#' sliced_wasserstein(M1, M2, thetas = cardinal_axes)
#' first_two_axes <- cardinal_axes[1:2, ]
#' sliced_wasserstein(M1, M2, thetas = first_two_axes)
sliced_wasserstein <- function(x, y, p = 2, thetas = NULL, L = 50, seed = NULL) {
  if (ncol(y) != ncol(x)) {
    stop("x and y should have the same number of columns.")
  }
  d <- ncol(x)

  if (!is.null(thetas)) {
    if (ncol(thetas) != d) {
      stop("thetas should be a matrix with the same number of columns as x, y")
    }
    L <- nrow(thetas)
  } else {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    thetas <- generate_directions(L, d)
  }

  total_cost <- 0
  for (ell in 1:L) {
    theta <- thetas[ell, ]
    x_theta <- dot(theta, x)
    y_theta <- dot(theta, y)
    total_cost <- total_cost + nw_corner_distance(x_theta, y_theta, p = p)^p
  }
  (total_cost / L)^(1/p)
}
