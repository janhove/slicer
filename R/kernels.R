#' Gaussian Radial Basis Function Kernel
#'
#' Builds a Gaussian radial basis function kernel matrix from
#' a distance matrix.
#'
#' @param D Matrix with pairwise (possibly squared) distances.
#' @param length_scale Length-scale parameter.
#' @param variance Kernel variance parameter.
#' @param squared If `TRUE`, `D` is assumed to contain squared
#'    distances. The default is `FALSE`.
#'
#' @return A kernel matrix with the same dimensions as `D`.
#' @export
#'
#' @examples
#' M <- rbind(c(0, 2, 1),
#'            c(2, 0, 2.5),
#'            c(1, 2.5, 0))
#' rbf(M)
rbf <- function(D, length_scale = 1, variance = 1, squared = FALSE) {
  if (squared) return(variance * exp(-D / (2 * length_scale^2)))
  variance * exp(-D^2/(2 * length_scale^2))
}
#' Matérn Kernel
#'
#' Builds a Matérn kernel from a distance matrix.
#'
#' @param D Matrix with pairwise (possibly squared) distances.
#' @param length_scale Length-scale parameter.
#' @param variance Kernel variance parameter.
#' @param nu Smoothness parameter. Defaults to 1.5.
#' @param squared If `TRUE`, `D` is assumed to contain squared
#'    distances. Defaults to `FALSE`.
#'
#' @return A kernel matrix with the same dimensions as `D`.
#' @export
#'
#' @examples
#' M <- rbind(c(0, 2, 1),
#'            c(2, 0, 2.5),
#'            c(1, 2.5, 0))
#' matern(M, nu = 0.5)
#' matern(M, nu = 1.5)
#' matern(M, nu = 2.5)
matern <- function(D, length_scale = 1, variance = 1, nu = 1.5,
  squared = FALSE) {
  if (squared) D <- sqrt(D)
  val <- variance * 2^(1 - nu) / gamma(nu) * (sqrt(2 * nu) * D / length_scale)^nu * besselK(sqrt(2 * nu) * D / length_scale, nu)
  if (is.matrix(val)) diag(val) <- variance
  val[D == 0] <- variance
  val
}
#' Power Exponential Kernel
#'
#' Builds a power exponential kernel from a distance matrix.
#'
#' @param D Matrix with pairwise (possibly squared) distances.
#' @param length_scale Length-scale parameter.
#' @param variance Kernel variance parameter.
#' @param expon Exponent. Defaults to 1.
#' @param squared If `TRUE`, `D` is assumed to contain squared
#'    distances. Defaults to `FALSE`.
#'
#' @return A kernel matrix with the same dimensions as `D`.
#' @export
#'
#' @examples
#' M <- rbind(c(0, 2, 1),
#'            c(2, 0, 2.5),
#'            c(1, 2.5, 0))
#' powexp(M, expon = 0.5)
#' powexp(M, expon = 1)
#' powexp(M, expon = 2)
powexp <- function(D, length_scale = 1, variance = 1, expon = 1,
                   squared = FALSE) {
  if (expon > 2) stop("expon must be at most 2.")
  if (expon <= 0) stop("expon must be at least 0.")
  if (squared) D <- sqrt(D)
  variance * exp(- (D / length_scale)^expon)
}
#' Sum of Kernels
#'
#' Builds kernel matrices from a list of distances matrices and computes
#' their sum.
#'
#' @noRd
#'
#' @param D_list A list of matrices with pairwise (possibly squared) distances.
#' @param length_scales A vector of length-scale parameters.
#' @param variances A vector of kernel variance parameters.
#' @param kernels A vector of kernels (`rbf`, `matern`, `powexp`).
#' @param params A vector of additional kernel parameters (smoothness parameter
#'  for Matérn, exponent for power exponential; ignored for RBF).
#' @param squared If `TRUE`, each distance matrix is assumed to
#'    contain squared distances. Defaults to `FALSE`.
#'
#' @return A kernel matrix with the same dimensions as the matrices in `D_list`.
#'
#' @examples
#' M1 <- rbind(c(0, 2, 1),
#'             c(2, 0, 2.5),
#'             c(1, 2.5, 0))
#' M2 <- rbind(c(0, 1, 4),
#'             c(1, 0, 3),
#'             c(4, 3, 0))
#' M3 <- rbind(c(0, 4, 8),
#'             c(4, 0, 1),
#'             c(8, 1, 0))
#' M_list <- list(M1, M2, M3)
#' K_sum <- kernel_sum(M_list, length_scales = c(1, 2, 0.5),
#'                     variances = c(2, 1/2, 0.1),
#'                     kernels = c("powexp", "matern", "rbf"), params = c(1, 2.3, NA))
#'
kernel_sum <- function(D_list, length_scales, variances, kernels, params,
                       squared = FALSE) {
  L <- length(D_list)
  if (length(length_scales) != L) {
    stop("length_scales, D_list don't have same length.")
  }
  if (length(variances) != L) stop("variances, D_list don't have same length.")
  if (length(kernels) == 1L) kernels <- rep(kernels, L)
  if (length(params) == 1L) params <- rep(params, L)
  N <- ncol(D_list[[1]])
  outcome <- matrix(0, nrow = N, ncol = N)
  for (ell in 1:L) {
    if (ncol(D_list[[ell]]) != N) {
      stop(paste0("Distance matrix ", ell,
                  " doesn't have the same number of columns as distance matrix 1."))
    }
    outcome <- outcome +
      my_kernel(D_list[[ell]], length_scale = length_scales[[ell]],
             variance = variances[[ell]], kernel = kernels[[ell]],
             param = params[[ell]], squared = squared)
  }
  outcome
}
