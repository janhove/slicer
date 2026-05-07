#' Predict Outcomes Using Gaussian Process Regression
#'
#' Computes posterior mean predictions from a Gaussian process regression model
#' given precomputed kernel matrices. The implementation follows Algorithm 2.1
#' in Rasmussen and Williams (2006), with an optional centring step and
#' automatic jitter escalation to handle near-singular kernel matrices.
#'
#' @param Kxx Kernel matrix evaluated at the training inputs.
#' @param Kxstar Matrix of kernel values between test and training inputs. Each row corresponds to one test input.
#' @param y_train Vector with training outcomes.
#' @param centre If `TRUE`, the training outcomes are centred around their mean.
#' @param lambda2 Kernel noise variance.
#'
#' @return A vector with posterior mean predictions for the test inputs.
#' @export
#'
#' @examples
#' N1 <- 25
#' N2 <- 10
#' x_train <- seq(-pi, pi, length.out = N1)
#' x_test  <- runif(N2, -pi, pi)
#' distance <- outer(c(x_train, x_test), c(x_train, x_test), "-") |> abs()
#' kernel <- rbf(distance, length_scale = 1, variance = 1)
#' Kxx    <- kernel[1:N1, 1:N1]
#' Kxstar <- kernel[(N1+1):(N1+N2), 1:N1]
#' y_train <- ifelse(x_train == 0, 2*pi, sin(2*pi*x_train) / x_train)
#' curve(sin(2*pi*x)/x, from = -pi, to = pi)
#' points(x_train, y_train)
#' points(x_test, gpr_predict(Kxx, Kxstar, y_train), pch = 16)
gpr_predict <- function(Kxx, Kxstar, y_train, centre = TRUE, lambda2 = 1e-6) {
  if (centre) {
    y_mean <- mean(y_train)
    y_train <- y_train - y_mean
  } else {
    y_mean <- 0
  }
  n <- nrow(Kxx)
  U <- tryCatch(
    chol(Kxx + lambda2 * diag(n)),
    error = function(e) {
      warning("Cholesky decomposition unsuccessful. Trying jitter (lambda2) escalation.")
      for (multiplier in c(10, 100, 1000, 10000)) {
        result <- tryCatch(
          chol(Kxx + (lambda2 * multiplier) * diag(n)),
          error = function(e) NULL
        )
        if (!is.null(result)) return(result)
      }
      stop("Matrix not positive definite even after jitter (lambda2) escalation")
    }
  )
  alpha <- backsolve(U, backsolve(U, y_train, transpose = TRUE))
  as.vector((Kxstar %*% alpha) + y_mean)
}

#' Estimate Gaussian process hyperparameters when using Gaussian RBF kernel
#'
#' The hyperparameters are estimated by minimising the model's marginal
#' negative log-likelihood using the BFGS algorithm.
#'
#' @noRd
#'
#' @param D2 Matrix with squared pairwise distances between the training objects only.
#' @param y_train Vector with outcomes for the training objects.
#' @param centre Set to `TRUE` (default) when centring training outcomes around their mean.
#' @param use_gradient If `TRUE`, closed-form expressions for the RBF's gradient are used.
#'                     Else, the optimiser uses the finite-differences method.
#' @param runs Number of independent attempts to find a minimum.
#' @param verbose If `TRUE`, progress is shown.
#'
#' @return A list containing the three estimated hyperparameters
#'         (kernel variance, kernel length-scale, noise variance) as well
#'         as the marginal negative log-likelihood achieved.
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
    D2, y_train, centre = TRUE, use_gradient = FALSE, runs = 10, verbose = TRUE
) {
  # Precomputations ------------------------------------------------------------
  n <- length(y_train)

  if (centre) {
    y_train <- y_train - mean(y_train)
  }

  med <- find_median(D2[1:n, 1:n]) |> sqrt()

  # Cacheing -------------------------------------------------------------------
  last_log_params <- NULL
  last_K <- NULL
  last_U <- NULL
  last_a <- NULL

  # Auxiliary functions --------------------------------------------------------

  # Look up/set in cache. After the nll() call, gradient() will call the same
  # log_params; we can save some costly computations by storing their results.
  nll <- function(log_params) {
    if (!identical(log_params, last_log_params)) {

      va <- exp(log_params[1])
      ls <- exp(log_params[2])
      lambda2 <- exp(log_params[3])

      K <- rbf(D2, length_scale = ls, variance = 1)
      K_reg <- va * K + lambda2 * diag(n)

      U <- tryCatch(chol(K_reg), error = function(e) NULL)
      if (is.null(U)) return(Inf)

      a <- backsolve(U, backsolve(U, y_train, transpose = TRUE))

      # Update cache
      last_log_params <<- log_params
      last_K <<- K
      last_U <<- U
      last_a <<- a
    }

    0.5 * crossprod(y_train, last_a) + sum(log(diag(last_U))) + n/2 * log(2 * pi)
  }

  gradient <- function(log_params) {
    if (!identical(log_params, last_log_params)) {
      va <- exp(log_params[1])
      ls <- exp(log_params[2])
      lambda2 <- exp(log_params[3])

      K <- rbf(D2, length_scale = ls, variance = 1)
      K_reg <- va * K + lambda2 * diag(n)

      U <- tryCatch(chol(K_reg), error = function(e) NULL)
      if (is.null(U)) return(Inf)

      a <- backsolve(U, backsolve(U, y_train, transpose = TRUE))

      # Update cache
      last_log_params <<- log_params
      last_K <<- K
      last_U <<- U
      last_a <<- a
    }

    va <- exp(log_params[1])
    ls <- exp(log_params[2])
    lambda2 <- exp(log_params[3])
    part_sigma2 <- va * last_K
    part_lengthscale <- va * D2 / ls^2 * last_K
    part_lambda2 <- diag(lambda2, n)

    jacobian <- vector("numeric", length = 3L)
    jacobian[[1]] <- -0.5 * (crossprod(last_a, part_sigma2 %*% last_a) - matrix_trace(backsolve(last_U, backsolve(last_U, part_sigma2, transpose = TRUE))))
    jacobian[[2]] <- -0.5 * (crossprod(last_a, part_lengthscale %*% last_a) - matrix_trace(backsolve(last_U, backsolve(last_U, part_lengthscale, transpose = TRUE))))
    jacobian[[3]] <- -0.5 * (crossprod(last_a, part_lambda2 %*% last_a) - matrix_trace(backsolve(last_U, backsolve(last_U, part_lambda2, transpose = TRUE))))

    jacobian
  }

  # Hyperparameter search ------------------------------------------------------
  best <- NULL
  for (i in seq_len(runs)) {
    if (verbose) message(paste0("Hyperparameter search ", i, " of ", runs, "."))
    init <- c(variance = stats::runif(1, log(0.1), log(10)),
              length_scale = stats::runif(1, log(med) + log(0.1), log(med) + log(10)),
              lambda2 = stats::runif(1, -10, -2))

    if (use_gradient) {
      fit <- tryCatch(
        stats::optim(init, nll, gr = gradient, method = "BFGS"),
        error = function(e) { NULL }
      )
    } else {
      fit <- tryCatch(
        stats::optim(init, nll, method = "BFGS"),
        error = function(e) { NULL }
      )
    }

    if (!is.null(fit) && (is.null(best) || fit$value < best$value)) {
      if (verbose && !is.null(best)) message(paste0("Current optimum improved from ",
        round(best$value, options()$digits), " to ",
        round(fit$value, options()$digits), "."))
      best <- fit
    }
  }

  list(
    variance = exp(best$par[1]),
    length_scale = exp(best$par[2]),
    lambda2 = exp(best$par[3]),
    nll = best$value
  )
}
#' Estimate Gaussian process hyperparameters when using multiple Gaussian RBF kernels
#'
#' This functions works exactly like \link{find_gpr_hyperparameters()},
#' but takes a list of squared-distance matrices rather than a single such
#' matrix.
#'
#' @noRd
find_gpr_hyperparameters_multiple <- function(
    D2_list, y_train, centre = TRUE,
    use_gradient = TRUE, runs = 10, verbose = TRUE
) {
  # Precomputations ------------------------------------------------------------
  L <- length(D2_list)
  n <- length(y_train)
  meds <- vapply(D2_list, find_median, FUN.VALUE = 0.0) |> sqrt()

  if (centre) {
    y_train <- y_train - mean(y_train)
  }

  # Cacheing -------------------------------------------------------------------
  last_log_params <- NULL
  last_K_list <- vector("list", L)
  last_U <- NULL
  last_a <- NULL

  # Auxiliary functions --------------------------------------------------------
  nll <- function(log_params) {
    if (!identical(log_params, last_log_params)) {
      va <- exp(log_params[1:L])
      ls <- exp(log_params[(L+1):(2*L)])
      lambda2 <- exp(log_params[2*L + 1])
      K <- matrix(0, ncol = n, nrow = n)
      for (ell in 1:L) {
        last_K_list[[ell]] <<- rbf(D2_list[[ell]], variance = 1, length_scale = ls[[ell]])
        K <- K + va[[ell]] * last_K_list[[ell]]
      }
      K <- K + lambda2 * diag(n)
      last_U <<- tryCatch(chol(K), error = function(e) NULL)
      if (is.null(last_U)) return(Inf)
      last_a <<- backsolve(last_U, backsolve(last_U, y_train, transpose = TRUE))
      last_log_params <<- log_params
    }

    0.5 * crossprod(y_train, last_a) + sum(log(diag(last_U))) + n/2 * log(2 * pi)
  }

  gradient <- function(log_params) {
    va_vec <- exp(log_params[1:L])
    ls_vec <- exp(log_params[(L+1):(2*L)])
    lambda2 <- exp(log_params[2*L + 1])

    if (!identical(log_params, last_log_params)) {
      K <- matrix(0, ncol = n, nrow = n)
      for (ell in 1:L) {
        last_K_list[[ell]] <<- rbf(D2_list[[ell]], variance = 1, length_scale = ls_vec[[ell]])
        K <- K + va_vec[[ell]] * last_K_list[[ell]]
      }
      K <- K + lambda2 * diag(n)
      last_U <<- tryCatch(chol(K), error = function(e) NULL)
      if (is.null(last_U)) return(Inf)
      last_a <<- backsolve(last_U, backsolve(last_U, y_train, transpose = TRUE))
      last_log_params <<- log_params
    }

    jacobian <- vector("numeric", 2*L)
    for (ell in 1:L) {
      va <- va_vec[[ell]]
      ls <- ls_vec[[ell]]
      part_sigma2 <- va * last_K_list[[ell]]
      part_lengthscale <- va * D2_list[[ell]] / ls^2 * last_K_list[[ell]]

      jacobian[[ell]] <- -0.5 * (crossprod(last_a, part_sigma2 %*% last_a) - matrix_trace(backsolve(last_U, backsolve(last_U, part_sigma2, transpose = TRUE))))
      jacobian[[L + ell]] <- -0.5 * (crossprod(last_a, part_lengthscale %*% last_a) - matrix_trace(backsolve(last_U, backsolve(last_U, part_lengthscale, transpose = TRUE))))
    }
    part_lambda2 <- diag(lambda2, n)
    jacobian[[2 * L + 1]] <-  -0.5 * (crossprod(last_a, part_lambda2 %*% last_a) - matrix_trace(backsolve(last_U, backsolve(last_U, part_lambda2, transpose = TRUE))))
    jacobian
  }

  # Hyperparameter search ------------------------------------------------------
  best <- NULL
  for (i in seq_len(runs)) {
    if (verbose) message(paste0("Hyperparameter search ", i, " of ", runs, "."))
    init <- c(variance = stats::runif(L, log(0.1), log(10)),
              length_scale = stats::runif(L, log(meds) + log(0.1), log(meds) + log(10)),
              lambda2 = stats::runif(1, -10, -2))

    if (use_gradient) {
      fit <- tryCatch(
        stats::optim(init, nll, gr = gradient, method = "BFGS"),
        error = function(e) NULL
      )
    } else {
      fit <- tryCatch(
        stats::optim(init, nll, method = "BFGS"),
        error = function(e) NULL
      )
    }

    if (!is.null(fit) && (is.null(best) || fit$value < best$value)) {
      if (verbose && !is.null(best)) message(paste0("Current optimum improved from ",
        round(best$value, options()$digits), " to ",
        round(fit$value, options()$digits), "."))
      best <- fit
    }
  }

  list(
    variance = exp(best$par[1:L]),
    length_scale = exp(best$par[(L+1):(2*L)]),
    lambda2 = exp(best$par[2*L + 1]),
    nll = best$value
  )
}
#' Gaussian Process Fitting to Test Data with Tuned Hyperparameters
#'
#' This function generates predictions using a Gaussian process model
#' with a Gaussian RBF kernel. The hyperparameters are tuned using the
#' training data by minimising the model's negative marginal log-likelihood.
#'
#' @param D2 A matrix containing squared pairwise distances between all objects.
#' @param training_idx A vector with the row (and column) numbers corresponding to the training entries in `D2`.
#' @param test_idx A vector with the row (and column) numbers corresponding to the test entries in `D2`.
#' @param y_train A vector with the training outcomes.
#' @param y_test Optionally, a vector with the test outcomes.
#' @param centre If `TRUE`, the training outcomes are centred around their mean before fitting the model. The mean is then added back to the predictions at the end.
#' @param use_gradient If `TRUE`, closed-form expressions for the RBF's gradient are used.
#'                     Else, the optimiser uses the finite-differences method.
#' @param runs Number of independent attempts to find a minimum when optimising the hyperparameters.
#' @param verbose If `TRUE`, progress is shown on the console.
#'
#' @return A list containing the predictions for the test objects,
#'     the root mean squared error (if the true test outcomes are provided),
#'     and the three tuned hyperparameter values.
#' @export
#'
#' @examples
#' N1 <- 25
#' N2 <- 10
#' x_train <- seq(-pi, pi, length.out = N1)
#' x_test  <- runif(N2, -pi, pi)
#' y_train <- x_train * plogis(x_train) * cos(x_train)
#' y_test <- x_test * plogis(x_test) * cos(x_test)
#' D2 <- outer(c(x_train, x_test), c(x_train, x_test), "-")^2
#' fit <- fit_gpr(D2, seq_len(N1), N1 + seq_len(N2), y_train, y_test, runs = 50)
#' curve(x * plogis(x) * cos(x), -pi, pi)
#' points(x_train, y_train, pch = 1)
#' points(x_test, fit$test_predictions, pch = 16)
#' fit$RMSE
fit_gpr <- function(
    D2, training_idx, test_idx, y_train, y_test = NULL,
    centre = TRUE, use_gradient = TRUE, runs = 10, verbose = TRUE) {
  D2_train <- D2[training_idx, training_idx]
  params <- find_gpr_hyperparameters(
    D2_train, y_train, centre = centre, use_gradient = use_gradient,
    runs = runs, verbose = verbose
  )
  variance <- params$variance
  length_scale <- params$length_scale
  lambda2 <- params$lambda2

  my_kernel <- rbf(D2, length_scale = length_scale, variance = variance)
  Kxx <- my_kernel[training_idx, training_idx]
  Kxstar <- my_kernel[test_idx, training_idx]
  predictions <- gpr_predict(Kxx, Kxstar, y_train, centre = centre, lambda2 = lambda2)

  list(
    test_predictions = predictions,
    RMSE = if (is.null(y_test)) NA else sqrt(mean((predictions - y_test)^2)),
    length_scale = length_scale,
    variance = variance,
    lambda2 = lambda2
  )
}
#' Gaussian Process Fitting to Test Data Using Multiple Kernels with Tuned Hyperparameters
#'
#' This is the counterpart to \link{fit_gpr} when multiple squared-distance
#' matrices are available. A kernel is then built as a linear combination
#' of base Gaussian RBF kernels. The hyperparameters are tuned by minimising
#' the model's negative marginal log-likelihood on the training data.
#'
#' @param D2_list A list of matrices containing squared pairwise distances between all objects.
#' @param training_idx A vector with the row (and column) numbers corresponding to the training entries in `D2`.
#' @param test_idx A vector with the row (and column) numbers corresponding to the test entries in `D2`.
#' @param y_train A vector with the training outcomes.
#' @param y_test Optionally, a vector with the test outcomes.
#' @param centre If `TRUE`, the training outcomes are centred around their mean before fitting the model. The mean is then added back to the predictions at the end.
#' @param use_gradient If `TRUE`, closed-form expressions for the RBF's gradient are used.
#'                     Else, the optimiser uses the finite-differences method.
#' @param runs Number of independent attempts to find a minimum when optimising the hyperparameters.
#' @param verbose If `TRUE`, progress is shown on the console.
#'
#' @return A list containing the predictions for the test objects,
#'     the root mean squared error (if the true test outcomes are provided),
#'     and the tuned hyperparameter values.
#' @export
#'
#' @examples
#' N1 <- 25
#' N2 <- 10
#' x_train1 <- runif(N1, -pi, pi)
#' x_train2 <- runif(N1, 0, 1)
#' x_test1  <- runif(N2, -pi, pi)
#' x_test2  <- runif(N2, 0, 1)
#' y_train <- x_train2 * plogis(x_train1) * cos(x_train1)
#' y_test <- x_test2 * plogis(x_test1) * cos(x_test1)
#' D2_1 <- outer(c(x_train1, x_test1), c(x_train1, x_test1), "-")^2
#' D2_2 <- outer(c(x_train2, x_test2), c(x_train2, x_test2), "-")^2
#' fit <- fit_gpr_multiple(list(D2_1, D2_2), seq_len(N1), N1 + seq_len(N2), y_train, y_test, runs = 50)
#' plot(fit$test_predictions, y_test)
#' fit$RMSE
fit_gpr_multiple <- function(
    D2_list, training_idx, test_idx, y_train, y_test,
    centre = TRUE, use_gradient = TRUE, runs = 10, verbose = TRUE) {
  select_entries <- function(M, idx) {
    M[idx, idx]
  }
  D2_train_list <- lapply(D2_list, select_entries, training_idx)
  params <- find_gpr_hyperparameters_multiple(
    D2_train_list, y_train, centre = centre, use_gradient = use_gradient,
    runs = runs, verbose = verbose
  )
  variance <- params$variance
  length_scale <- params$length_scale
  lambda2 <- params$lambda2

  my_kernel <- rbf_multiple(D2_list, length_scales = length_scale, variances = variance)
  Kxx <- my_kernel[training_idx, training_idx]
  Kxstar <- my_kernel[test_idx, training_idx]
  predictions <- gpr_predict(Kxx, Kxstar, y_train, centre = centre, lambda2 = lambda2)

  list(
    test_predictions = predictions,
    RMSE = sqrt(mean((predictions - y_test)^2)),
    length_scale = length_scale,
    variance = variance,
    lambda2 = lambda2
  )
}
