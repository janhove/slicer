#' Predict Outcomes Using Gaussian Process Regression
#'
#' Computes posterior mean predictions from a Gaussian process regression model
#' given precomputed kernel matrices. The implementation follows Algorithm 2.1
#' in Rasmussen and Williams (2006), with an optional centring step and
#' automatic jitter escalation to handle near-singular kernel matrices.
#'
#' @param Kxx Kernel matrix evaluated at the training inputs.
#' @param Kxstar Matrix of kernel values between test and training inputs.
#'  Each row corresponds to one test input.
#' @param y_train Vector with training outcomes.
#' @param centre If `TRUE`, the training outcomes are centred around their mean.
#' @param lambda2 Kernel noise variance.
#' @param Kxstarstar Optional. A square matrix of kernel values among the
#'  test inputs. If `NULL` (default), only posterior means are computed
#'  and returned as a vector. Else both the posterior mean vector and covariance
#'  matrix are returned.
#'
#' @return If `Kxstarstar` is `NULL` (default), a vector with posterior mean
#'  predictions for the test inputs. Else, a list with a vector
#'  `mean` with posterior mean predictions and a matrix `var` with the posterior
#'  covariance for the test inputs.
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
#' Kxstarstar <- kernel[(N1+1):(N1+N2), (N1+1):(N1+N2)]
#' y_train <- ifelse(x_train == 0, 2*pi, sin(2*pi*x_train) / x_train)
#' curve(sin(2*pi*x)/x, from = -pi, to = pi)
#' points(x_train, y_train)
#' # Only means: don't use Kxstarstar
#' points(x_test, gpr_predict(Kxx, Kxstar, y_train), pch = 16)
#' # Add credible intervals around predictions:
#' gpr_fit <- gpr_predict(Kxx, Kxstar, y_train, Kxstarstar = Kxstarstar)
#' segments(x0 = x_test,
#'          y0 = gpr_fit$mean - 2 * sqrt(diag(gpr_fit$var)),
#'          y1 = gpr_fit$mean + 2 * sqrt(diag(gpr_fit$var)), col = "#4DAF4A")
gpr_predict <- function(Kxx, Kxstar, y_train, centre = TRUE, lambda2 = 1e-6, Kxstarstar = NULL) {
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
      warning("Cholesky decomposition unsuccessful. Trying jitter escalation.")
      for (multiplier in c(10, 100, 1000, 10000)) {
        result <- tryCatch(
          chol(Kxx + (lambda2 * multiplier) * diag(n)),
          error = function(e) NULL
        )
        if (!is.null(result)) return(result)
      }
      stop("Matrix not positive definite even after jitter escalation")
    }
  )

  alpha <- backsolve(U, backsolve(U, y_train, transpose = TRUE))
  mean_pred <- as.vector((Kxstar %*% alpha) + y_mean)

  if (is.null(Kxstarstar)) return(mean_pred)

  if (is.vector(Kxstarstar)) {
    Kxstarstar <- matrix(Kxstarstar)
  }
  if (is.vector(Kxstar)) {
    Kxstar <- matrix(Kxstar, nrow = 1)
  }

  if (nrow(Kxstarstar) != ncol(Kxstarstar)) {
    stop("Kxstarstar must be square.")
  }

  V <- backsolve(U, t(Kxstar), transpose = TRUE)
  var_pred <- Kxstarstar - crossprod(V)
  var_pred <- (var_pred + t(var_pred)) / 2 # enforce symmetry, just to be safe
  list(mean = mean_pred, var = var_pred)
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
    kernel = "rbf",
    runs = 10, verbose = TRUE
) {
  # Set kernel -----------------------------------------------------------------
  if (!(kernel %in% c("rbf", "matern05", "matern15", "matern25"))) {
    stop("Invalid kernel name.")
  }

  if (kernel != "rbf") use_gradient <- FALSE

  # Precomputations ------------------------------------------------------------
  n <- length(y_train)

  if (centre) {
    y_train <- y_train - mean(y_train)
  }

  med <- find_median(D2) |> sqrt()

  # Cacheing -------------------------------------------------------------------
  last_log_params <- NULL
  last_unit_K <- NULL
  last_U <- NULL
  last_a <- NULL

  # Auxiliary functions --------------------------------------------------------
  # Look up/set in cache. After the nll() call, gradient() will call the same
  # log_params; we can save some costly computations by storing their results.
  update_cache <- function(log_params) {
    if (identical(log_params, last_log_params)) return(invisible(NULL))
    va <- exp(log_params[1])
    ls  <- exp(log_params[2])
    lambda2 <- exp(log_params[3])
    unit_K <- my_kernel(D2, length_scale = ls, variance = 1, kernel = kernel)
    K <- va * unit_K + lambda2 * diag(n)
    last_U <<- tryCatch(chol(K), error = function(e) NULL)
    if (!is.null(last_U)) {
      last_a <<- backsolve(last_U, backsolve(last_U, y_train, transpose = TRUE))
      last_unit_K <<- unit_K
    }
    last_log_params <<- log_params
  }

  nll <- function(log_params) {
    update_cache(log_params)
    if (is.null(last_U)) return(Inf)
    drop(0.5 * crossprod(y_train, last_a) + sum(log(diag(last_U))) + n/2 * log(2 * pi))
  }

  # Only for RBF kernel at the moment
  gradient <- function(log_params) {
    update_cache(log_params)
    if (is.null(last_U)) return(rep(Inf, 3))

    va <- exp(log_params[1])
    ls <- exp(log_params[2])
    lambda2 <- exp(log_params[3])

    W <- chol2inv(last_U)

    jacobian <- numeric(3)
    DK <- va * last_unit_K
    jacobian[[1]] <- -0.5 * (drop(last_a %*% DK %*% last_a) - sum(W * DK))
    jacobian[[2]] <- -0.5 * (drop(last_a %*% (DK * D2 / ls^2) %*% last_a) -
                               sum(W * DK * D2 / ls^2))
    jacobian[[3]] <- -0.5 * lambda2 * (drop(last_a %*% last_a) - sum(diag(W)))

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
      if (verbose && is.null(best)) message(paste0("Optimum set at ",
        round(fit$value, options()$digits), "."))
      best <- fit
    }
  }

  list(
    kernel = kernel,
    variance = exp(best$par[1]),
    length_scale = exp(best$par[2]),
    lambda2 = exp(best$par[3]),
    nll = best$value
  )
}
#' Estimate Gaussian Process Hyperparameters When Using Multiple Kernels
#'
#' This functions works exactly like \link{find_gpr_hyperparameters()},
#' but takes a list of squared-distance matrices rather than a single such
#' matrix. You can also specify a vector of kernels rather than a single one.
#'
#' @noRd
find_gpr_hyperparameters_multiple <- function(
    D2_list, y_train, centre = TRUE,
    use_gradient = TRUE, kernels = "rbf", runs = 10, verbose = TRUE
) {
  # Set kernel -----------------------------------------------------------------
  if (any(!(kernels %in% c("rbf", "matern05", "matern15", "matern25")))) {
    stop("There is at least one invalid kernel name.")
  }

  if (length(kernels) == 1L) kernels <- rep(kernels, length(D2_list))

  if (any(!(kernels %in% "rbf"))) use_gradient <- FALSE

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
  update_cache <- function(log_params) {
    if (identical(log_params, last_log_params)) return(invisible(NULL))
    va_vec <- exp(log_params[1:L])
    ls_vec  <- exp(log_params[(L+1):(2*L)])
    lambda2 <- exp(log_params[2*L + 1])
    K <- Reduce(`+`, mapply(function(K_l, va, ls, krnl)
      va * my_kernel(K_l, variance = 1, length_scale = ls, kernel = krnl),
      D2_list, va_vec, ls_vec, kernels, SIMPLIFY = FALSE))
    K <- K + lambda2 * diag(n)
    last_U <<- tryCatch(chol(K), error = function(e) NULL)
    if (!is.null(last_U)) {
      last_a <<- backsolve(last_U, backsolve(last_U, y_train, transpose = TRUE))
      last_K_list <<- mapply(function(K_l, ls, krnl)
        my_kernel(K_l, variance = 1, length_scale = ls, kernel = krnl),
        D2_list, ls_vec, kernels, SIMPLIFY = FALSE)
    }
    last_log_params <<- log_params
  }

  nll <- function(log_params) {
    update_cache(log_params)
    if (is.null(last_U)) return(Inf)
    drop(0.5 * crossprod(y_train, last_a) + sum(log(diag(last_U))) + n/2 * log(2 * pi))
  }

  # Only used when using RBF
  gradient <- function(log_params) {
    update_cache(log_params)
    if (is.null(last_U)) return(rep(Inf, 2*L + 1))

    va_vec <- exp(log_params[1:L])
    ls_vec <- exp(log_params[(L+1):(2*L)])
    lambda2 <- exp(log_params[2*L + 1])

    W <- chol2inv(last_U)

    jacobian <- numeric(2*L + 1)
    for (ell in 1:L) {
      DK <- last_K_list[[ell]] * va_vec[[ell]]
      jacobian[[ell]] <- -0.5 * (drop(last_a %*% DK %*% last_a) - sum(W * DK))
      jacobian[[L + ell]] <- -0.5 * (drop(last_a %*% (DK * D2_list[[ell]] / ls_vec[[ell]]^2) %*% last_a) - sum(W * DK * D2_list[[ell]] / ls_vec[[ell]]^2))
    }
    jacobian[[2*L + 1]] <- -0.5 * lambda2 * (drop(last_a %*% last_a) - sum(diag(W)))

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
      if (verbose && is.null(best)) message(paste0("Optimum set at ",
        round(fit$value, options()$digits), "."))
      best <- fit
    }
  }

  list(
    kernels = kernels,
    variance = exp(best$par[1:L]),
    length_scale = exp(best$par[(L+1):(2*L)]),
    lambda2 = exp(best$par[2*L + 1]),
    nll = best$value
  )
}
#' Gaussian Process Fitting to Test Data with Tuned Hyperparameters
#'
#' This function generates predictions using a Gaussian process model.
#' The hyperparameters are tuned using the training data by minimising
#' the model's negative marginal log-likelihood.
#'
#' @noRd
#'
#' @param D2 A matrix containing squared pairwise distances between all objects.
#' @param training_idx A vector with the row (and column) numbers corresponding to the training entries in `D2`.
#' @param test_idx A vector with the row (and column) numbers corresponding to the test entries in `D2`.
#' @param y_train A vector with the training outcomes.
#' @param y_test Optionally, a vector with the test outcomes.
#' @param centre If `TRUE`, the training outcomes are centred around their mean before fitting the model. The mean is then added back to the predictions at the end.
#' @param use_gradient If `TRUE`, closed-form expressions for the RBF's gradient are used.
#'  Else, the optimiser uses the finite-differences method.
#'  For non-RBF kernels, this is reset to `FALSE`.
#' @param kernel A kernel (`rbf`, `matern05`, `matern15`, `matern25`).
#' @param runs Number of independent attempts to find a minimum when optimising the hyperparameters.
#' @param verbose If `TRUE`, progress is shown on the console.
#'
#' @return A list containing the predictions for the test objects as well as
#'     their variance, the root mean squared error
#'     (if the true test outcomes are provided),
#'     and the three tuned hyperparameter values.
#'
#' @examples
#' N1 <- 25
#' N2 <- 10
#' x_train <- seq(-pi, pi, length.out = N1)
#' x_test  <- runif(N2, -pi, pi)
#' y_train <- x_train * plogis(x_train) * cos(x_train)
#' y_test <- x_test * plogis(x_test) * cos(x_test)
#' D2 <- outer(c(x_train, x_test), c(x_train, x_test), "-")^2
#' fit <- fit_gpr_single(D2, seq_len(N1), N1 + seq_len(N2), y_train, y_test, runs = 50)
#' curve(x * plogis(x) * cos(x), -pi, pi)
#' points(x_train, y_train, pch = 1)
#' points(x_test, fit$test_predictions, pch = 16)
#' fit$RMSE
#' fit <- fit_gpr_single(D2, seq_len(N1), N1 + seq_len(N2), y_train, y_test, runs = 10,
#'    kernel = "matern15")
#' fit$RMSE
fit_gpr_single <- function(
    D2, training_idx, test_idx, y_train, y_test = NULL,
    centre = TRUE, use_gradient = TRUE,
    kernel = "rbf",
    runs = 10, verbose = TRUE) {
  # Set kernel -----------------------------------------------------------------
  if (!(kernel %in% c("rbf", "matern05", "matern15", "matern25"))) {
    stop("Invalid kernel name.")
  }

  if (kernel != "rbf") use_gradient <- FALSE

  D2_train <- D2[training_idx, training_idx]
  params <- find_gpr_hyperparameters(
    D2_train, y_train, centre = centre, use_gradient = use_gradient,
    kernel = kernel,
    runs = runs, verbose = verbose
  )
  variance <- params$variance
  length_scale <- params$length_scale
  lambda2 <- params$lambda2

  K <- my_kernel(D2, length_scale = length_scale, variance = variance, kernel = kernel)
  Kxx <- K[training_idx, training_idx]
  Kxstar <- K[test_idx, training_idx]
  Kxstarstar <- K[test_idx, test_idx]
  predictions <- gpr_predict(Kxx, Kxstar, y_train, centre = centre,
    lambda2 = lambda2, Kxstarstar)

  list(
    kernels = kernel,
    test_predictions = predictions$mean,
    test_variance = predictions$var,
    RMSE = if (is.null(y_test)) NA else sqrt(mean((predictions$mean - y_test)^2)),
    length_scale = unname(length_scale),
    scaling_factor = unname(variance),
    noise_variance = unname(lambda2),
    nll = params$nll
  )
}
#' Gaussian Process Fitting to Test Data Using Multiple Kernels with Tuned Hyperparameters
#'
#' This is the counterpart to \link{fit_gpr} when multiple squared-distance
#' matrices are available. A kernel is then built as a linear combination
#' of base kernels. The hyperparameters are tuned by minimising
#' the model's negative marginal log-likelihood on the training data.
#'
#' @noRd
#'
#' @param D2_list A list of matrices containing squared pairwise distances between all objects.
#' @param training_idx A vector with the row (and column) numbers corresponding to the training entries in `D2`.
#' @param test_idx A vector with the row (and column) numbers corresponding to the test entries in `D2`.
#' @param y_train A vector with the training outcomes.
#' @param y_test Optionally, a vector with the test outcomes.
#' @param centre If `TRUE`, the training outcomes are centred around their mean before fitting the model. The mean is then added back to the predictions at the end.
#' @param use_gradient If `TRUE`, closed-form expressions for the RBF's gradient are used.
#'    Else, the optimiser uses the finite-differences method.
#'    Only used when RBF kernels are used.
#' @param kernels A vector of kernels (`rbf`, `matern05`, `matern15`, `matern25`).
#' @param runs Number of independent attempts to find a minimum when optimising the hyperparameters.
#' @param verbose If `TRUE`, progress is shown on the console.
#'
#' @return A list containing the predictions for the test objects,
#'     the root mean squared error (if the true test outcomes are provided),
#'     the tuned hyperparameter values, and the negative log marginal likelihood (nll).
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
#' fit <- fit_gpr_multiple(list(D2_1, D2_2), seq_len(N1), N1 + seq_len(N2),
#'   y_train, y_test, runs = 50, kernels = c("rbf", "matern15"))
#' plot(fit$test_predictions, y_test)
#' fit$RMSE
fit_gpr_multiple <- function(
    D2_list, training_idx, test_idx, y_train, y_test = NULL,
    centre = TRUE, use_gradient = TRUE, kernels = "rbf", runs = 10, verbose = TRUE) {
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

  if (length(kernels) == 1L) kernels <- rep(kernels, length(D2_list))

  K <- kernel_sum(D2_list, length_scales = length_scale, variances = variance, kernels = kernels)
  Kxx <- K[training_idx, training_idx]
  Kxstar <- K[test_idx, training_idx]
  Kxstarstar <- K[test_idx, test_idx]
  predictions <- gpr_predict(Kxx, Kxstar, y_train, centre = centre,
    lambda2 = lambda2, Kxstarstar)

  list(
    kernels = kernels,
    test_predictions = predictions$mean,
    test_variance = predictions$var,
    RMSE = if (is.null(y_test)) NA else sqrt(mean((predictions$mean - y_test)^2)),
    length_scale = unname(length_scale),
    scaling_factor = unname(variance),
    noise_variance = unname(lambda2),
    nll = params$nll
  )
}

#' Gaussian Process Fitting to Test Data with Tuned Hyperparameters
#'
#' This function generates predictions using a Gaussian process model
#' with one or several Gaussian RBF kernels. The hyperparameters are tuned using the
#' training data by minimising the model's negative marginal log-likelihood.
#'
#' @param D2 A matrix, or a list of matrices, containing squared pairwise distances between all objects.
#' @param training_idx A vector with the row (and column) numbers corresponding to the training entries in `D2`.
#' @param test_idx A vector with the row (and column) numbers corresponding to the test entries in `D2`.
#' @param y_train A vector with the training outcomes.
#' @param y_test Optionally, a vector with the test outcomes.
#' @param centre If `TRUE`, the training outcomes are centred around their mean before fitting the model. The mean is then added back to the predictions at the end.
#' @param use_gradient If `TRUE`, closed-form expressions for the RBF's gradient are used.
#'                     Else, the optimiser uses the finite-differences method.
#'                     Ignored when non-RBF kernels are used.
#' @param kernels A vector with kernels (`"rbf", "matern05", "matern15", "matern25"`).
#' @param runs Number of independent attempts to find a minimum when optimising the hyperparameters.
#' @param cores Number of cores used for parallel processing.
#' @param verbose If `TRUE`, progress is shown on the console. Only works if `cores == 1L`.
#'
#' @return A list containing the predictions for the test objects as well as
#'     their variance, the root mean squared error
#'     (if the true test outcomes are provided),
#'     and the three tuned hyperparameter values.
#' @export
#'
#' @examples
#' # Multiple kernels
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
#' # Single core
#' fit <- fit_gpr(list(D2_1, D2_2), seq_len(N1), N1 + seq_len(N2),
#'   y_train, y_test, runs = 50)
#' plot(y_test, fit$test_predictions)
#' abline(a = 0, b = 1, lty = 1)
#' segments(x0 = y_test,
#'   y0 = fit$test_predictions - 2 * sqrt(diag(fit$test_variance)),
#'   y1 = fit$test_predictions + 2 * sqrt(diag(fit$test_variance)), lty = 2)
#' fit
#' # Multiple cores
#' fit <- fit_gpr(list(D2_1, D2_2), seq_len(N1), N1 + seq_len(N2),
#'   y_train, y_test, runs = 50, cores = 2)
#' fit
#'
#' # Single kernel
#' N1 <- 40
#' N2 <- 10
#' x_train <- seq(-pi, pi, length.out = N1)
#' x_test  <- runif(N2, -pi, pi)
#' y_train <- x_train * plogis(x_train) * cos(x_train) + rnorm(N1, sd = 0.5)
#' y_test <- x_test * plogis(x_test) * cos(x_test)
#' D2 <- outer(c(x_train, x_test), c(x_train, x_test), "-")^2
#' fit <- fit_gpr(D2, seq_len(N1), N1 + seq_len(N2), y_train, y_test,
#'   runs = 50L, cores = 2)
#' curve(x * plogis(x) * cos(x), -pi, pi,
#'   ylim = range(
#'     c(y_train, fit$test_predictions + 2 * sqrt(diag(fit$test_variance)),
#'                fit$test_predictions - 2 * sqrt(diag(fit$test_variance)))
#'    )
#'  )
#' points(x_train, y_train, pch = 1)
#' points(x_test, fit$test_predictions, pch = 16)
#' segments(x0 = x_test,
#'   y0 = fit$test_predictions - 2 * sqrt(diag(fit$test_variance)),
#'   y1 = fit$test_predictions + 2 * sqrt(diag(fit$test_variance)))
#' fit
#' fit$RMSE
#'
#' fit <- fit_gpr(D2, seq_len(N1), N1 + seq_len(N2), y_train, y_test,
#'   runs = 50L, cores = 2, kernel = "matern15")
#' curve(x * plogis(x) * cos(x), -pi, pi,
#'   ylim = range(
#'     c(y_train, fit$test_predictions + 2 * sqrt(diag(fit$test_variance)),
#'                fit$test_predictions - 2 * sqrt(diag(fit$test_variance)))
#'    )
#'  )
#' points(x_train, y_train, pch = 1)
#' points(x_test, fit$test_predictions, pch = 16)
#' segments(x0 = x_test,
#'   y0 = fit$test_predictions - 2 * sqrt(diag(fit$test_variance)),
#'   y1 = fit$test_predictions + 2 * sqrt(diag(fit$test_variance)))
#' fit
#' fit$RMSE
fit_gpr <- function(
    D2, training_idx, test_idx, y_train, y_test = NULL,
    centre = TRUE, use_gradient = TRUE, kernels = "rbf",
    runs = 10L, cores = 1L, verbose = TRUE) {
  if (cores <= 1L) {
    if (is.list(D2)) {
      my_fit <- fit_gpr_multiple(D2, training_idx, test_idx, y_train,
        y_test, centre = centre, use_gradient = use_gradient, runs = runs,
        verbose = verbose, kernels = kernels)
      return(my_fit)
    }
    my_fit <- fit_gpr_single(D2, training_idx, test_idx, y_train,
      y_test, centre = centre, use_gradient = use_gradient, runs = runs,
      verbose = verbose, kernel = kernels)
    return(my_fit)
  }

  # Parallel processing
  available_cores <- parallel::detectCores()
  if (cores > available_cores - 1) {
    cores <- available_cores - 1
    warning(paste0("Number of cores clamped to ", cores, "."))
  }
  runs_per_core <- ceiling(runs / cores)
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl))
  parallel::clusterEvalQ(cl, library(slicer))

  results <- parallel::parLapply(
    cl, seq_len(cores),
    function(i, D2, training_idx, test_idx, y_train, y_test,
             centre, use_gradient, runs_per_core, kernels = kernels) {

      if (is.list(D2)) {
        fit_gpr_multiple(
          D2, training_idx, test_idx, y_train, y_test,
          centre = centre, use_gradient = use_gradient,
          runs = runs_per_core, verbose = FALSE, kernels = kernels
        )
      } else {
        fit_gpr_single(
          D2, training_idx, test_idx, y_train, y_test,
          centre = centre, use_gradient = use_gradient,
          runs = runs_per_core, verbose = FALSE, kernel = kernels
        )
      }
    },
    D2 = D2, training_idx = training_idx, test_idx = test_idx,
    y_train = y_train, y_test = y_test,
    centre = centre, use_gradient = use_gradient,
    runs_per_core = runs_per_core, kernels = kernels
  )
  best_idx <- which.min(sapply(results, `[[`, "nll"))
  results[[best_idx]]
}
