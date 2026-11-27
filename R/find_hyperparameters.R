#' Estimate Gaussian Process Hyperparameters
#'
#' The hyperparameters are estimated by minimising the model's marginal
#' negative log-likelihood using the BFGS algorithm.
#' If you feed multiple distance matrices to this function,
#' a sum of kernels is built, with each summand having a separate
#' length-scale and variances.
#'
#' @noRd
#'
#' @param D Matrix, or list of matrices, with pairwise distances between the
#'    training objects only, or
#' @param y_train Vector with outcomes for the training objects.
#' @param centre Set to `TRUE` (default) when centring training outcomes
#'    around their mean.
#' @param use_gradient If `TRUE`, closed-form expressions for the RBF's
#'    gradient are used. Else, the optimiser uses the finite-differences method.
#'    If the kernel used is not the RBF, this value is currently reset to `FALSE`.
#' @param kernels A vector with the names of the kernel functions to be applied
#'    to each distance matrix.
#'    `"rbf"` for the Gaussian RBF kernel (default),
#'    `"matern"` for the Matérn covariance function, and `"powexp"` for the
#'    power exponential kernel.
#'    If a single kernel is specified but a list of distance matrices is supplied,
#'    the same kernel is used for all distance matrices.
#' @param params A vector with additional kernel parameters
#'    (smoothness for Matérn, exponent for power exponential; ignored for RBF).
#'    If a single kernel is specified but a list of distance matrices is supplied,
#'    the same kernel is used for all distance matrices.
#' @param squared Set to `TRUE` if `D` contains squared distances. Defaults to `FALSE`.
#' @param runs Number of independent attempts to find a minimum.
#' @param verbose If `TRUE`, progress is shown.
#' @param maxit Maximum number of iterations during optimisation.
#'
#' @return A list containing the kernels used, the estimated hyperparameters
#'    (kernel variance and kernel length-scale for each kernel as well as
#'    noise variance) as well as the marginal negative log-likelihood achieved.
find_hyperparameters <- function(
    D, y_train, centre = TRUE, use_gradient = TRUE, kernels = "rbf",
    params = NA, squared = FALSE, runs = 10, verbose = TRUE, maxit = 500
) {
  # Set kernel -----------------------------------------------------------------
  if (any(!(kernels %in% c("rbf", "powexp", "matern")))) {
    stop("There is at least one invalid kernel name.")
  }

  # For later: Generalise use of gradients to other kernels
  use_gradient <- use_gradient && all(kernels == "rbf")

  # Precomputations ------------------------------------------------------------
  if (!is.list(D)) D <- list(D)
  L <- length(D)

  if (centre) y_train <- y_train - mean(y_train)
  n <- length(y_train)
  log_yvar <- log(stats::var(y_train))

  meds <- vapply(D, find_median, FUN.VALUE = 0.0)
  if (squared) meds <- sqrt(meds)

  # Cacheing (initialisation) --------------------------------------------------
  last_log_params <- NULL
  last_K_list <- vector("list", L)
  last_U <- NULL
  last_a <- NULL

  # Auxiliary functions --------------------------------------------------------
  update_cache <- function(log_params) {
    if (identical(log_params, last_log_params)) return(invisible(NULL))
    va_vec <- exp(log_params[1:L])
    ls_vec  <- exp(log_params[(L+1):(2*L)])
    lambda2 <- exp(log_params[2*L + 1])
    K_list <- mapply(
      function(D_ell, va, ls, krnl, param, sqd)
        va * my_kernel(D_ell, variance = 1, length_scale = ls,
                       kernel = krnl, param = param, squared = sqd),
      D, va_vec, ls_vec, kernels, params, squared,
      SIMPLIFY = FALSE
    )
    K <- Reduce(`+`, K_list)
    diag(K) <- diag(K) + lambda2
    last_U <<- tryCatch(chol(K), error = function(e) NULL)
    if (!is.null(last_U)) {
      last_a <<- backsolve(last_U, backsolve(last_U, y_train, transpose = TRUE))
      last_K_list <<- K_list
    }
    last_log_params <<- log_params
  }

  nll <- function(log_params) {
    update_cache(log_params)
    if (is.null(last_U)) return(Inf)
    drop(0.5 * crossprod(y_train, last_a) + sum(log(diag(last_U))) + n/2 * log(2 * pi))
  }

  # Only used when using RBF. May want to generalise in future.
  gradient <- function(log_params) {
    update_cache(log_params)
    if (is.null(last_U)) return(rep(Inf, 2 * L + 1))

    va_vec <- exp(log_params[1:L])
    ls_vec <- exp(log_params[(L+1):(2*L)])
    lambda2 <- exp(log_params[2*L + 1])

    M <- chol2inv(last_U) - tcrossprod(last_a)

    jacobian <- numeric(2*L + 1)

    for (ell in 1:L) {
      DK <- va_vec[[ell]] * last_K_list[[ell]]
      jacobian[[ell]] <- 0.5 * sum(M * DK)
      jacobian[[L + ell]] <- 0.5 * sum(M * DK * D[[ell]]^2) / ls_vec[[ell]]^2
    }
    jacobian[[2*L + 1]] <- 0.5 * lambda2 * sum(diag(M))

    jacobian
  }

  # Hyperparameter search ------------------------------------------------------
  ctrl <- list(reltol = 1e-6, maxit = maxit)
  best <- NULL
  for (i in seq_len(runs)) {
    if (verbose) message(paste0("Hyperparameter search ", i, " of ", runs, "."))
    # Randomly sampled initial values, tailored to scale of data
    init <- c(variance = log_yvar - log(L) + stats::runif(L, log(0.1), log(10)),
              length_scale = log(meds) + stats::runif(L, log(0.1), log(10)),
              lambda2 = log_yvar + stats::runif(1, -10, -2))

    if (use_gradient) {
      fit <- tryCatch(
        stats::optim(init, nll, gr = gradient, method = "BFGS", control = ctrl),
        error = function(e) NULL
      )
    } else {
      fit <- tryCatch(
        stats::optim(init, nll, method = "BFGS", control = ctrl),
        error = function(e) NULL
      )
    }

    if (!is.null(fit) && (is.null(best) || fit$value < best$value)) {
      if (verbose && !is.null(best)) {
        message(paste0("Current optimum improved from ",
                       round(best$value, 4), " to ",
                       round(fit$value,  4), "."))
      }
      if (verbose && is.null(best)) {
        message(paste0("Optimum set at ",
                       round(fit$value,  4), "."))
      }
      best <- fit
    }
  }
  if (is.null(best)) stop("All optimisation runs failed.")

  list(
    kernels = kernels,
    params = params,
    variance = exp(best$par[1:L]),
    length_scale = exp(best$par[(L+1):(2*L)]),
    lambda2 = exp(best$par[2*L + 1]),
    nll = best$value
  )
}
