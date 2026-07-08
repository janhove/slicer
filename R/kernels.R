#' Gaussian Radial Basis Function Kernel
#'
#' Builds a Gaussian radial basis function kernel matrix from
#' a distance matrix.
#'
#' @param D Matrix with pairwise (possibly squared) distances.
#' @param length_scale Length-scale parameter.
#' @param variance Kernel variance parameter.
#' @param squared If `TRUE` (default), `D` is assumed to contain squared
#'    distances. Set to `FALSE` if this is not the case.
#'
#' @return A kernel matrix with the same dimensions as `D`.
#' @export
#'
#' @examples
#' M <- rbind(c(0, 2, 1),
#'            c(2, 0, 2.5),
#'            c(1, 2.5, 0))
#' rbf(M, squared = FALSE)
rbf <- function(D, length_scale = 1, variance = 1, squared = TRUE) {
  if (squared) return(variance * exp(-D / (2 * length_scale^2)))
  variance * exp(-D^2/(2 * length_scale^2))
}
#' Sum of Gaussian Radial Basis Function Kernels
#'
#' Builds Gaussian radial basis function kernel matrices from a list of
#' distances matrices and computes their sum.
#'
#' @param D_list A list of matrices with pairwise (possibly squared) distances.
#' @param length_scales A vector of length-scale parameters.
#' @param variances A vector of kernel variance parameters.
#' @param squared If `TRUE` (default), each distance matrix is assumed to
#'    contain squared distances. If `FALSE`, each distance matrix is assumed to
#'    contain raw distances.
#'
#' @return A kernel matrix with the same dimensions as the matrices in `D_list`.
#' @export
#'
#' @examples
#' M1 <- rbind(c(0, 2, 1),
#'             c(2, 0, 2.5),
#'             c(1, 2.5, 0))
#' M2 <- rbind(c(0, 1, 4),
#'             c(1, 0, 3),
#'             c(4, 3, 0))
#' M_list <- list(M1, M2)
#' rbf_multiple(M_list, length_scales = c(1, 2), variances = c(2, 1/2), squared = FALSE)
rbf_multiple <- function(D_list, length_scales, variances, squared = TRUE) {
  L <- length(D_list)
  if (length(length_scales) != L) {
    stop("length_scales, D_list don't have same length.")
  }
  if (length(variances) != L) stop("variances, D_list don't have same length.")
  N <- ncol(D_list[[1]])
  outcome <- matrix(0, nrow = N, ncol = N)
  for (ell in 1:L) {
    if (ncol(D_list[[ell]]) != N) {
      stop(paste0("Distance matrix ", ell,
        " doesn't have the same number of columns as distance matrix 1."))
    }
    outcome <- outcome +
      rbf(D_list[[ell]], length_scale = length_scales[[ell]],
          variance = variances[[ell]], squared = squared)
  }
  outcome
}
#' Matérn Kernels
#'
#' Builds a Matérn kernel (with smoothness parameter 0.5, 1.5 or 2.5) from a
#' distance matrix.
#'
#' @param D Matrix with pairwise (possibly squared) distances.
#' @param length_scale Length-scale parameter.
#' @param variance Kernel variance parameter.
#' @param tau Smoothness parameter. needs to be 0.5, 1.5 (default) or 2.5.
#' @param squared If `TRUE` (default), `D` is assumed to contain squared
#'    distances. Set to `FALSE` if this is not the case.
#'
#' @return A kernel matrix with the same dimensions as `D`.
#' @export
#'
#' @examples
#' M <- rbind(c(0, 2, 1),
#'            c(2, 0, 2.5),
#'            c(1, 2.5, 0))
#' matern(M, tau = 0.5, squared = FALSE)
#' matern(M, tau = 1.5, squared = FALSE)
#' matern(M, tau = 2.5, squared = FALSE)
matern <- function(D, length_scale = 1, variance = 1, tau = 1.5,
  squared = TRUE) {
  if (!(tau %in% c(0.5, 1.5, 2.5))) stop("tau must be 0.5, 1.5 or 2.5.")
  if (squared) D <- sqrt(D)
  if (tau == 0.5) {
    return(variance * exp(-D / (2 * length_scale)))
  }
  if (tau == 1.5) {
    return(variance * (1 + sqrt(3) * D / length_scale) *
             exp(-sqrt(3) * D / length_scale))
  }
  variance *
    (1 + sqrt(5) * D / length_scale + 5 * D^2 / (3 * length_scale^2)) *
    exp(-sqrt(5) * D / length_scale)
}
#' Sum of Matérn Kernels
#'
#' Builds Matérn kernel matrices from a list of distances matrices and computes
#' their sum.
#'
#' @param D_list A list of matrices with pairwise (possibly squared) distances.
#' @param length_scales A vector of length-scale parameters.
#' @param variances A vector of kernel variance parameters.
#' @param taus A vector of smoothness parameters.
#'    Values need to be in 0.5, 1.5 and 2.5.
#'    If only a single value is provided, it will be used for all components.
#' @param squared If `TRUE` (default), each distance matrix is assumed to
#'    contain squared distances. If `FALSE`, each distance matrix is assumed to
#'    contain raw distances.
#'
#' @return A kernel matrix with the same dimensions as the matrices in `D_list`.
#' @export
#'
#' @examples
#' M1 <- rbind(c(0, 2, 1),
#'             c(2, 0, 2.5),
#'             c(1, 2.5, 0))
#' M2 <- rbind(c(0, 1, 4),
#'             c(1, 0, 3),
#'             c(4, 3, 0))
#' M_list <- list(M1, M2)
#' matern_multiple(M_list, length_scales = c(1, 2), variances = c(2, 1/2),
#'   taus = c(0.5, 2.5), squared = FALSE)
matern_multiple <- function(D_list, length_scales, variances, taus,
  squared = TRUE) {
  L <- length(D_list)
  if (length(length_scales) != L) {
    stop("length_scales, D_list don't have same length.")
  }
  if (length(variances) != L) stop("variances, D_list don't have same length.")
  if (length(taus) == 1L) taus <- rep(taus, L)
  if (length(taus) != L) stop("taus, D_list don't have same length.")
  N <- ncol(D_list[[1]])
  outcome <- matrix(0, nrow = N, ncol = N)
  for (ell in 1:L) {
    if (ncol(D_list[[ell]]) != N) {
      stop(paste0("Distance matrix ", ell,
        " doesn't have the same number of columns as distance matrix 1."))
    }
    outcome <- outcome +
      matern(D_list[[ell]], length_scale = length_scales[[ell]],
          variance = variances[[ell]], tau = taus[[ell]], squared = squared)
  }
  outcome
}
