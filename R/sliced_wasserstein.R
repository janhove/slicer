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
#' Pairwise Wasserstein Distances Along Projection Directions Between Several Empirical Distributions
#'
#' This function computes the squared 2-Wasserstein distances between the projections
#' of several empirical directions along specified projection directions.
#' The projection directions can optionally be transformed by means of a linear map.
#' The functions can output the squared distances along each projection direction,
#' or the squared distances averaged over projection directions.
#'
#' @param distributions A list of matrices representing empirical distributions.
#' @param thetas A matrix, each row of which represents a projection direction.
#' @param A Optionally, a matrix used to transform each projection direction.
#' @param verbose If `TRUE`, show progress.
#' @param keep_projections If `TRUE`, the distance matrix for each projection direction is output.
#'                         If `FALSE`, the distance matrices for the different projection directions are averaged.
#'
#' @return A list of squared-distance matrices, one for each projection direction (if `keep_projections = TRUE`);
#'         otherwise, a matrix with the averaged squared distances.
#' @export
#'
#' @examples
#' M1 <- matrix(rnorm(50), ncol = 5)
#' M2 <- matrix(rnorm(150), ncol = 5)
#' M3 <- matrix(rnorm(250), ncol = 5)
#' # Sliced Wasserstein:
#' my_directions <- generate_directions(20, 5)
#' compute_all_distances(list(M1, M2, M3), my_directions, keep_projections = FALSE)
#' # Marginal Wasserstein distances:
#' marginal_wass <- compute_all_distances(list(M1, M2, M3), diag(1, 5), keep_projections = TRUE)
#' marginal_wass[[3]] # along third dimension
compute_all_distances <- function(
    distributions, thetas, A = NULL, verbose = TRUE, keep_projections = TRUE
) {
  d <- ncol(thetas)
  L <- nrow(thetas)
  N <- length(distributions)

  if (!is.null(A)) {
    if (verbose) message("Transforming projection directions...")
    thetas <- thetas %*% A
  }

  if (verbose) message("Projecting distributions...")
  projections <- lapply(distributions, project_and_sort, thetas)

  if (verbose) {
    message("Computing sliced Wasserstein distances...")
    pb <- utils::txtProgressBar(min = 0, max = L, style = 3)
  }

  distance_list <- vector("list", L)
  for (ell in 1:L) {

    distance_matrix <- matrix(0, nrow = N, ncol = N)
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        current_distance <- nw_corner_distance(projections[[i]][, ell], projections[[j]][, ell], presorted = TRUE)^2
        distance_matrix[i, j] <- current_distance
        distance_matrix[j, i] <- current_distance
      }
    }
    distance_list[[ell]] <- distance_matrix
    if (verbose) utils::setTxtProgressBar(pb, ell)
  }
  if (verbose) cat("\n")

  if (keep_projections) return(distance_list)

  M <- matrix(0, nrow = N, ncol = N)
  for (ell in seq_len(L)) {
    M <- M + distance_list[[ell]]
  }
  M / L
}
