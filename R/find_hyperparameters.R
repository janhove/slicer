#' Estimate Gaussian Process Hyperparameters When Using Multiple Kernels
#'
#' This functions works exactly like \link{find_gpr_hyperparameters()},
#' but takes a list of squared-distance matrices rather than a single such
#' matrix. You can also specify a vector of kernels rather than a single one.
#'
#' @noRd
find_gpr_hyperparameters_multiple <- function(
    D2_list, y_train, centre = TRUE,
    use_gradient = TRUE, kernels = "rbf", params = NA,
    runs = 10, verbose = TRUE, maxit = 500
) {
  # Set kernel -----------------------------------------------------------------
  if (any(!(kernels %in% c("rbf", "powexp", "matern")))) {
    stop("There is at least one invalid kernel name.")
  }

  if (length(kernels) == 1L) kernels <- rep(kernels, length(D2_list))
  if (length(params) == 1L) params <- rep(params, length(D2_list))

  use_gradient <- use_gradient && all(kernels == "rbf")

  # Precomputations ------------------------------------------------------------
  L <- length(D2_list)
  n <- length(y_train)
  meds <- vapply(D2_list, find_median, FUN.VALUE = 0.0) |> sqrt()
  log_yvar <- log(stats::var(y_train))

  if (centre) y_train <- y_train - mean(y_train)

  # Cacheing -------------------------------------------------------------------
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
      function(K_l, va, ls, krnl, param)
        va * my_kernel(K_l, variance = 1, length_scale = ls,
                       kernel = krnl, param = param),
      D2_list, va_vec, ls_vec, kernels, params,
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

  # Only used when using RBF
  # May want to generalise in future.
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
      jacobian[[L + ell]] <- 0.5 * sum(M * DK * D2_list[[ell]]) / ls_vec[[ell]]^2
    }
    jacobian[[2*L + 1]] <- 0.5 * lambda2 * sum(diag(M))

    jacobian
  }

  # Hyperparameter search ------------------------------------------------------
  ctrl <- list(reltol = 1e-6, maxit = maxit)
  best <- NULL
  for (i in seq_len(runs)) {
    if (verbose) message(paste0("Hyperparameter search ", i, " of ", runs, "."))
    init <- c(variance = log_yvar - log(L) * stats::runif(L, log(0.1), log(10)),
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
      if (verbose && !is.null(best)) message(paste0("Current optimum improved from ",
        round(best$value, 4), " to ",
        round(fit$value,  4), "."))
      if (verbose && is.null(best)) message(paste0("Optimum set at ",
        round(fit$value,  4), "."))
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

#' Estimate Gaussian Process Hyperparameters
#'
#' The hyperparameters are estimated by minimising the model's marginal
#' negative log-likelihood using the BFGS algorithm.
#'
#' @noRd
#'
#' @param D2 Matrix with squared pairwise distances between the
#'    training objects only.
#' @param y_train Vector with outcomes for the training objects.
#' @param centre Set to `TRUE` (default) when centring training outcomes
#'    around their mean.
#' @param use_gradient If `TRUE`, closed-form expressions for the RBF's
#'    gradient are used. Else, the optimiser uses the finite-differences method.
#'    If the kernel used is not the RBF, this value is currently reset to `FALSE`.
#' @param kernel `rbf` for the Gaussian RBF kernel (default);
#'    `matern05`, `matern15` or `matern25` for the
#'    Matérn kernel with smoothness parameter 0.5, 1.5 or 2.5.
#' @param param Additional kernel parameter (smoothness for Matérn, exponent for power exponential; ignored for RBF).
#' @param runs Number of independent attempts to find a minimum.
#' @param verbose If `TRUE`, progress is shown.
#'
#' @return A list containing the kernel used, the three estimated hyperparameters
#'    (kernel variance, kernel length-scale, noise variance) as well
#'    as the marginal negative log-likelihood achieved.
#'
#' @examples
#' N1 <- 25
#' x_train <- seq(-pi, pi, length.out = N1)
#' y_train <- ifelse(x_train == 0, 2*pi, sin(2*pi*x_train) / x_train)
#' D2_train <- outer(x_train, x_train, "-")^2
#' find_gpr_hyperparameters(D2_train, y_train,
#'   use_gradient = FALSE, runs = 20, verbose = FALSE)
#' find_gpr_hyperparameters(D2_train, y_train,
#'   use_gradient = TRUE, runs = 20, verbose = FALSE)
find_gpr_hyperparameters <- function(
    D2, y_train, centre = TRUE, use_gradient = FALSE,
    kernel = "rbf", param = NA,
    runs = 10, verbose = TRUE, maxit = 500
) {
  fit <- find_gpr_hyperparameters_multiple(
    list(D2), y_train, centre = centre, use_gradient = use_gradient,
    kernels = kernel, params = param, runs = runs, verbose = verbose,
    maxit = maxit)
  list(
    kernel = kernel,
    variance = fit$variance,
    length_scale = fit$length_scale,
    lambda2 = fit$lambda2,
    nll = fit$nll
  )
}
