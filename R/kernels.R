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
#' Matérn Kernel
#'
#' Builds a Matérn kernel from a distance matrix.
#'
#' @param D Matrix with pairwise (possibly squared) distances.
#' @param length_scale Length-scale parameter.
#' @param variance Kernel variance parameter.
#' @param nu Smoothness parameter. Defaults to 1.5.
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
#' matern(M, nu = 0.5, squared = FALSE)
#' matern(M, nu = 1.5, squared = FALSE)
#' matern(M, nu = 2.5, squared = FALSE)
matern <- function(D, length_scale = 1, variance = 1, nu = 1.5,
  squared = TRUE) {
  if (squared) D <- sqrt(D)
  val <- variance * 2^(1 - nu) / gamma(nu) * (sqrt(2 * nu) * D / length_scale)^nu * besselK(sqrt(2 * nu) * D / length_scale, nu)
  if (is.matrix(val)) diag(val) <- variance
  val[D == 0] <- variance
  val
}
#' Sum of Matérn Kernels
#'
#' Builds Matérn kernel matrices from a list of distances matrices and computes
#' their sum.
#'
#' @param D_list A list of matrices with pairwise (possibly squared) distances.
#' @param length_scales A vector of length-scale parameters.
#' @param variances A vector of kernel variance parameters.
#' @param nus A vector of smoothness parameters.
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
#'   nus = c(0.5, 2.5), squared = FALSE)
matern_multiple <- function(D_list, length_scales, variances, nus,
  squared = TRUE) {
  L <- length(D_list)
  if (length(length_scales) != L) {
    stop("length_scales, D_list don't have same length.")
  }
  if (length(variances) != L) stop("variances, D_list don't have same length.")
  if (length(nus) == 1L) nus <- rep(nus, L)
  if (length(nus) != L) stop("nus, D_list don't have same length.")
  N <- ncol(D_list[[1]])
  outcome <- matrix(0, nrow = N, ncol = N)
  for (ell in 1:L) {
    if (ncol(D_list[[ell]]) != N) {
      stop(paste0("Distance matrix ", ell,
        " doesn't have the same number of columns as distance matrix 1."))
    }
    outcome <- outcome +
      matern(D_list[[ell]], length_scale = length_scales[[ell]],
          variance = variances[[ell]], nu = nus[[ell]], squared = squared)
  }
  outcome
}
#' Power Exponential Kernel
#'
#' Builds a power exponential kernel from a distance matrix.
#'
#' @param D Matrix with pairwise (possibly squared) distances.
#' @param length_scale Length-scale parameter.
#' @param variance Kernel variance parameter.
#' @param expon Exponent. Defaults to 1.
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
#' powexp(M, expon = 0.5, squared = FALSE)
#' powexp(M, expon = 1, squared = FALSE)
#' powexp(M, expon = 2, squared = FALSE)
powexp <- function(D, length_scale = 1, variance = 1, expon = 1,
                   squared = TRUE) {
  if (expon > 2) stop("expon must be at most 2.")
  if (expon <= 0) stop("expon must be at least 0.")
  if (squared) D <- sqrt(D)
  variance * exp(- (D / length_scale)^expon)
}
#' Sum of Power Exponential Kernels
#'
#' Builds power exponential kernel matrices from a list of distances matrices
#' and computes their sum.
#'
#' @param D_list A list of matrices with pairwise (possibly squared) distances.
#' @param length_scales A vector of length-scale parameters.
#' @param variances A vector of kernel variance parameters.
#' @param expons A vector of exponents.
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
#' powexp_multiple(M_list, length_scales = c(1, 2), variances = c(2, 1/2),
#'   expons = c(0.5, 1.5), squared = FALSE)
powexp_multiple <- function(D_list, length_scales, variances, expons,
                            squared = TRUE) {
  L <- length(D_list)
  if (length(length_scales) != L) {
    stop("length_scales, D_list don't have same length.")
  }
  if (length(variances) != L) stop("variances, D_list don't have same length.")
  if (length(expons) == 1L) nus <- rep(expons, L)
  if (length(expons) != L) stop("nus, D_list don't have same length.")
  N <- ncol(D_list[[1]])
  outcome <- matrix(0, nrow = N, ncol = N)
  for (ell in 1:L) {
    if (ncol(D_list[[ell]]) != N) {
      stop(paste0("Distance matrix ", ell,
                  " doesn't have the same number of columns as distance matrix 1."))
    }
    outcome <- outcome +
      powexp(D_list[[ell]], length_scale = length_scales[[ell]],
             variance = variances[[ell]], expon = expons[[ell]], squared = squared)
  }
  outcome
}
#' Sum of Kernels
#'
#' Builds kernel matrices from a list of distances matrices and computes
#' their sum.
#'
#' @param D_list A list of matrices with pairwise (possibly squared) distances.
#' @param length_scales A vector of length-scale parameters.
#' @param variances A vector of kernel variance parameters.
#' @param kernels A vector of kernels (`rbf`, `matern`, `powexp`).
#' @param params A vector of additional kernel parameters (smoothness parameter
#'  for Matérn, exponent for power exponential; ignored for RBF).
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
#' M3 <- rbind(c(0, 4, 8),
#'             c(4, 0, 1),
#'             c(8, 1, 0))
#' M_list <- list(M1^2, M2^2, M3^2)
#' K_sum <- kernel_sum(M_list, length_scales = c(1, 2, 0.5),
#'                     variances = c(2, 1/2, 0.1),
#'                     kernels = c("powexp", "matern", "rbf"), params = c(1, 2.3, NA))
#'
kernel_sum <- function(D_list, length_scales, variances, kernels, params, squared = TRUE) {
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
             param = params[[ell]])
  }
  outcome
}
