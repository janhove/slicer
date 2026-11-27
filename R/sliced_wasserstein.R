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
#' @param thetas Optionally, a matrix, each row of which represents
#'    a projection direction.
#' @param L If no `thetas` are provided, `L` random projection directions are
#'    generated.
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
  if (is.vector(x)) x <- matrix(x, ncol = 1)
  if (is.vector(y)) y <- matrix(y, ncol = 1)

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
#' This function computes Wasserstein distances between the projections
#' of several empirical directions along specified projection directions.
#' The projection directions can optionally be transformed by means of a linear map.
#' The functions can output the distances along each projection direction,
#' or the distances averaged over projection directions.
#'
#' @section Squared Wasserstein distances:
#'
#' The Wasserstein distance of order p = 2 as well as the square root of
#' the Wasserstein distance of order p = 1 are Hilbertian distances between
#' univariate probability distributions.
#' Similarly, the sliced Wasserstein distance of order p = 2 as well as the
#' square root of the sliced Wasserstein distance of order p = 1 are Hilbertian
#' distances between multivariate distributions.
#' We can therefore slot in these distances into the Gaussian radial basis
#' function to obtain a valid kernel.
#' Downstream, the workflow is simplified somewhat if `compute_all_distances()`
#' return not the raw Wasserstein distances but the p-Wasserstein and p-sliced
#' Wasserstein distances _raised to the p-th power_.
#' This way, squared (sliced) 2-Wasserstein distances can directly be fed to,
#' say, `rbf()` with the default parameter setting `squared = TRUE`,
#' as can the raw (sliced) 1-Wasserstein distances.
#'
#' Recall that the p-sliced Wasserstein distance is computed as the p-th root
#' of the average of p-Wasserstein distances raised to the p-th power.
#'
#' @param distributions A list of matrices representing empirical distributions.
#' @param thetas A matrix, each row of which represents a projection direction.
#' @param A Optionally, a matrix used to transform each projection direction.
#' @param p Order of the Wasserstein distances. Defaults to 2.
#' @param verbose If `TRUE`, show progress.
#' @param keep_projections If `TRUE`, the distance matrix for each projection direction is output.
#'    If `FALSE`, the distance matrices for the different projection directions are averaged.
#' @param test_idx Optionally, a vector of indices. If supplied, skip the
#'    distance computations between distribution pairs whose indices both occur
#'    in this vector.
#'
#' @return If `keep_projections = TRUE`, a list of matrices containing
#'  p-Wasserstein distances raised to the p-th power,
#'  one for each projection direction;
#'  otherwise, a matrix with the averages of the p-Wasserstein distances
#'  raised to the p-th power.
#' @export
#'
#' @examples
#' M1 <- matrix(rnorm(50), ncol = 5)
#' M2 <- matrix(rnorm(50), ncol = 5)
#' M3 <- matrix(rnorm(250), ncol = 5)
#'
#' # Squared sliced Wasserstein of order 2
#' my_directions <- generate_directions(20, 5)
#' compute_all_distances(list(M1, M2, M3), my_directions,
#'   keep_projections = FALSE, verbose = FALSE)
#'
#' # To obtain raw sliced Wasserstein of order 2, take square root.
#'
#' # (Raw) sliced Wasserstein of order 1
#' compute_all_distances(list(M1, M2, M3), my_directions, p = 1,
#'   keep_projections = FALSE, verbose = FALSE)
#'
#' # Marginal squared 2-Wasserstein distances:
#' marginal_wass <- compute_all_distances(list(M1, M2, M3), diag(1, 5),
#'   keep_projections = TRUE, verbose = FALSE)
#' marginal_wass[[3]] # along third dimension
#'
#' # Reweight projection directions (shear transformation)
#' A <- diag(c(4, 0.5, 3, 2, 1))
#' shear_wass <- compute_all_distances(list(M1, M2, M3), diag(1, 5), A = A,
#'   keep_projections = TRUE, verbose = FALSE)
#' shear_wass[[3]]
compute_all_distances <- function(
    distributions, thetas = NULL, A = NULL, p = 2, verbose = TRUE,
    keep_projections = TRUE, test_idx = NULL
) {
  if (is.null(thetas)) {
    if (verbose) message("No projection directions (thetas) supplied. Assuming univariate distributions.")
    thetas <- matrix(1, 1, 1)
    keep_projections <- FALSE
  }

  d <- ncol(thetas)
  L <- nrow(thetas)
  N <- length(distributions)

  if (!is.null(A)) {
    if (verbose) message("Transforming...")
    thetas <- thetas %*% A
  }

  if (verbose && d*L > 1) message("Projecting...")
  projections <- lapply(distributions, project_and_sort, thetas)

  if (verbose && d*L > 1) {
    message("Computing distances...")
    pb <- utils::txtProgressBar(min = 0, max = L, style = 3)
  }

  distance_list <- vector("list", L)
  for (ell in 1:L) {

    distance_matrix <- matrix(0, nrow = N, ncol = N)
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        if (i %in% test_idx && j %in% test_idx) {
          distance_matrix[i, j] <- NA
          distance_matrix[j, i] <- NA
          next
        }
        current_distance <- nw_corner_distance(projections[[i]][, ell],
          projections[[j]][, ell], p = p, presorted = TRUE)^p
        distance_matrix[i, j] <- current_distance
        distance_matrix[j, i] <- current_distance
      }
    }
    distance_list[[ell]] <- distance_matrix
    if (verbose && d*L > 1) utils::setTxtProgressBar(pb, ell)
  }
  if (verbose) cat("\n")

  if (keep_projections) return(distance_list)

  M <- matrix(0, nrow = N, ncol = N)
  for (ell in seq_len(L)) {
    M <- M + distance_list[[ell]]
  }
  M / L
}
