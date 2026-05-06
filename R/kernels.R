#' Gaussian radial basis function kernel
#'
#' Builds a Gaussian radial basis function kernel matrix from
#' a distance matrix.
#'
#' @param D A matrix with pairwise (possibly squared) distances.
#' @param length_scale The length-scale parameter.
#' @param variance The kernel variance parameter.
#' @param squared If `TRUE` (default), `D` is assumed to contain squared distances.
#'                Set to `FALSE` if this is not the case.
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
#' Sum of Gaussian radial basis function kernels
#'
#' Builds Gaussian radial basis function kernel matrices from a list of
#' distances matrices and computes their sum.
#'
#' @param D_list A list of matrices with pairwise (possibly squared) distances.
#' @param length_scales A vector of length-scale parameters.
#' @param variances A vector of kernel variance parameters.
#' @param squared If `TRUE` (default), each distance matrix is assumed to contain squared distances.
#'                If `FALSE`, each distance matrix is assumed to contain raw distances.
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
  if (length(length_scales) != L) stop("length_scales, D_list don't have same length.")
  if (length(variances) != L) stop("variances, D_list don't have same length.")
  N <- ncol(D_list[[1]])
  outcome <- matrix(0, nrow = N, ncol = N)
  for (ell in 1:L) {
    if (ncol(D_list[[ell]]) != N) stop(paste0("Distance matrix ", ell, " doesn't have the same number of columns as distance matrix 1."))
    outcome <- outcome +
      rbf(D_list[[ell]], length_scale = length_scales[[ell]], variance = variances[[ell]], squared = squared)
  }
  outcome
}
