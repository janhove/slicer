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
#' @param centre If `TRUE` (defaults), the training outcomes are centred
#'  around their mean.
#' @param lambda2 Kernel noise variance.
#' @param Kxstarstar Optional. A square matrix of kernel values among the
#'  test inputs. If `NULL` (default), only posterior means are computed
#'  and returned as a vector. Else both the posterior mean vector and covariance
#'  matrix are returned.
#'
#' @return If `Kxstarstar` is `NULL` (default), a vector with posterior mean
#'  predictions for the test inputs. Else, a list with a vector
#'  `test_predictions` with posterior mean predictions and
#'  a matrix `test_variance` with the posterior covariance for the test inputs.
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
#'          y0 = gpr_fit$test_predictions - 2 * sqrt(diag(gpr_fit$test_variance)),
#'          y1 = gpr_fit$test_predictions + 2 * sqrt(diag(gpr_fit$test_variance)), col = "#4DAF4A")
gpr_predict <- function(Kxx, Kxstar, y_train, centre = TRUE, lambda2 = 1e-6,
                        Kxstarstar = NULL) {
  # To do: Input validation.

  n <- nrow(Kxx)

  if (is.vector(Kxstar)) Kxstar <- matrix(Kxstar, nrow = 1)

  if (centre) {
    y_mean <- mean(y_train)
    y_train <- y_train - y_mean
  } else {
    y_mean <- 0
  }

  base <- max(lambda2, 1e-10)
  U <- NULL
  used <- lambda2
  for (multiplier in c(1, 10, 100, 1000, 10000)) {
    used <- ifelse(multiplier == 1, lambda2, base * multiplier)
    A <- Kxx
    diag(A) <- diag(A) + used
    U <- tryCatch(chol(A), error = function(e) NULL)
    if (!is.null(U)) break
  }
  if (is.null(U)) {
    stop("Matrix not p.d. even after jitter escalation.")
  }
  if (used != lambda2) {
    warning(sprintf("Cholesky decomposition failed with lambda = %g; used %g instead.",
                    lambda2, used))
  }

  alpha <- backsolve(U, backsolve(U, y_train, transpose = TRUE))
  mean_pred <- drop(Kxstar %*% alpha) + y_mean

  if (is.null(Kxstarstar)) return(mean_pred)

  if (is.vector(Kxstarstar)) Kxstarstar <- matrix(Kxstarstar)
  if (nrow(Kxstarstar) != ncol(Kxstarstar)) stop("Kxstarstar must be square.")

  V <- backsolve(U, t(Kxstar), transpose = TRUE)
  var_pred <- Kxstarstar - crossprod(V)

  list(
    test_predictions = mean_pred,
    test_variance = var_pred,
    noise_variance = unname(used)
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
#' @param D A list of matrices containing pairwise distances between all objects.
#' @param training_idx A vector with the row (and column) numbers corresponding to the training entries in `D`.
#' @param test_idx A vector with the row (and column) numbers corresponding to the test entries in `D`.
#' @param y_train A vector with the training outcomes.
#' @param y_test Optionally, a vector with the test outcomes.
#' @param centre If `TRUE`, the training outcomes are centred around their mean before fitting the model. The mean is then added back to the predictions at the end.
#' @param use_gradient If `TRUE`, closed-form expressions for the RBF's gradient are used.
#'    Else, the optimiser uses the finite-differences method.
#'    Only used when RBF kernels are used.
#' @param kernels A vector of kernels (`rbf`, `matern`, `powexp`).
#' @param params A vector of kernel parameters (smoothness for Matern; exponent for power exponential; ignored for RBF).
#' @param squared If `TRUE`, the distance matrices are assumed to contain squared distances. Defaults to `FALSE`.
#' @param runs Number of independent attempts to find a minimum when optimising the hyperparameters.
#' @param verbose If `TRUE`, progress is shown on the console.
#'
#' @return A list containing the predictions for the test objects,
#'     the root mean squared error (if the true test outcomes are provided),
#'     the tuned hyperparameter values, and the negative log marginal likelihood (nll).
fit_gpr_one_core <- function(
    D, training_idx, test_idx, y_train, y_test = NULL,
    centre = TRUE, use_gradient = TRUE, kernels = "rbf",
    params = NA, squared = FALSE, runs = 10, verbose = TRUE, maxit = 500)
  {
  select_entries <- function(M, idx) {
    M[idx, idx]
  }
  D_train <- lapply(D, select_entries, training_idx)
  tuned_hyperparams <- find_hyperparameters(
    D_train, y_train, centre = centre, use_gradient = use_gradient,
    runs = runs, verbose = verbose, kernels = kernels, params = params,
    maxit = maxit, squared = squared
  )
  variance <- tuned_hyperparams$variance
  length_scale <- tuned_hyperparams$length_scale
  lambda2 <- tuned_hyperparams$lambda2

  if (length(kernels) == 1L) kernels <- rep(kernels, length(D))

  K <- kernel_sum(D, length_scales = length_scale, variances = variance,
                  kernels = kernels, params = params, squared = squared)
  Kxx <- K[training_idx, training_idx]
  Kxstar <- K[test_idx, training_idx]
  Kxstarstar <- K[test_idx, test_idx]
  predictions <- gpr_predict(Kxx, Kxstar, y_train, centre = centre,
    lambda2 = lambda2, Kxstarstar)

  list(
    kernels = kernels,
    params = tuned_hyperparams$params,
    test_predictions = predictions$test_predictions,
    test_variance = predictions$test_variance,
    RMSE = if (is.null(y_test)) NA else sqrt(mean((predictions$test_predictions - y_test)^2)),
    length_scale = unname(length_scale),
    scaling_factor = unname(variance),
    noise_variance = unname(lambda2),
    nll = tuned_hyperparams$nll
  )
}

#' Gaussian Process Fitting to Test Data with Tuned Hyperparameters
#'
#' This function generates predictions using a Gaussian process model
#' with one or several Gaussian RBF kernels. The hyperparameters are tuned using the
#' training data by minimising the model's negative marginal log-likelihood.
#'
#' @param D A matrix, or a list of matrices, containing pairwise distances between all objects.
#' @param training_idx A vector with the row (and column) numbers corresponding to the training entries in `D`.
#' @param test_idx A vector with the row (and column) numbers corresponding to the test entries in `D`.
#' @param y_train A vector with the training outcomes.
#' @param y_test Optionally, a vector with the test outcomes.
#' @param centre If `TRUE`, the training outcomes are centred around their mean before fitting the model. The mean is then added back to the predictions at the end.
#' @param use_gradient If `TRUE`, closed-form expressions for the RBF's gradient are used.
#'                     Else, the optimiser uses the finite-differences method.
#'                     Ignored when non-RBF kernels are used.
#' @param kernels A vector with kernels (`"rbf", "powexp", "matern"`).
#' @param params Kernel parameters (smoothness for Matern; exponent for power exponential; ignored for RBF).
#' @param squared Set to `TRUE` if the distances in `D` are squared. Defaults to `FALSE`.
#' @param runs Number of independent attempts to find a minimum when optimising the hyperparameters.
#' @param cores Number of cores used for parallel processing.
#' @param verbose If `TRUE`, progress is shown on the console. Only works if `cores == 1L`.
#' @param maxit Maximum number of iterations per optimisation run.
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
#' D1 <- outer(c(x_train1, x_test1), c(x_train1, x_test1), "-") |> abs()
#' D2 <- outer(c(x_train2, x_test2), c(x_train2, x_test2), "-") |> abs()
#'
#' # Single core
#' fit <- fit_gpr(list(D1, D2), seq_len(N1), N1 + seq_len(N2),
#'   y_train, y_test, runs = 50)
#' plot(y_test, fit$test_predictions)
#' abline(a = 0, b = 1, lty = 1)
#' segments(x0 = y_test,
#'   y0 = fit$test_predictions - 2 * sqrt(diag(fit$test_variance)),
#'   y1 = fit$test_predictions + 2 * sqrt(diag(fit$test_variance)), lty = 2)
#' fit
#'
#' # Multiple cores
#' fit <- fit_gpr(list(D1, D2), seq_len(N1), N1 + seq_len(N2),
#'  y_train, y_test, runs = 50, cores = 2L)
#' fit
#'
#' # Single kernel
#' N1 <- 40
#' N2 <- 10
#' x_train <- seq(-pi, pi, length.out = N1)
#' x_test  <- runif(N2, -pi, pi)
#' y_train <- x_train * plogis(x_train) * cos(x_train) + rnorm(N1, sd = 0.5)
#' y_test <- x_test * plogis(x_test) * cos(x_test)
#' D <- outer(c(x_train, x_test), c(x_train, x_test), "-") |> abs()
#' fit <- fit_gpr(D, seq_len(N1), N1 + seq_len(N2), y_train, y_test, runs = 50L)
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
    D, training_idx, test_idx, y_train, y_test = NULL,
    centre = TRUE, use_gradient = TRUE, kernels = "rbf", params = NA,
    squared = FALSE, runs = 10L, cores = 1L, verbose = TRUE, maxit = 500)
  {

  # Input cleaning -------------------------------------------------------------
  if (!is.list(D)) D <- list(D)
  L <- length(D)

  if (length(kernels) == 1L) kernels <- rep(kernels, L)
  if (length(params) == 1L)  params  <- rep(params, L)

  if (length(kernels) != L)  stop(paste0(length(kernels), " kernels specified; ", L, " needed."))
  if (length(params) != L)   stop(paste0(length(params), " params specified; ", L, " needed."))

  specify_params <- is.na(params) & kernels == "matern"
  if (sum(specify_params) > 0) {
    params[specify_params] <- 3/2
    warning("Unspecified Matern smoothness parameters set to 3/2.")
  }
  specify_params <- is.na(params) & kernels == "powexp"
  if (sum(specify_params) > 0) {
    params[specify_params] <- 1
    warning("Unspecified power exponents set to 1.")
  }

  # One core -------------------------------------------------------------------
  if (cores <= 1L) {
    my_fit <- fit_gpr_one_core(D, training_idx, test_idx, y_train,
      y_test, centre = centre, use_gradient = use_gradient, runs = runs,
      squared = squared,
      verbose = verbose, kernels = kernels, params = params, maxit = maxit)
    return(my_fit)
  }

  # Parallel processing --------------------------------------------------------
  available_cores <- parallel::detectCores()
  if (cores > available_cores - 1) {
    cores <- available_cores - 1
    warning(paste0("Number of cores clamped to ", cores, "."))
  }
  runs_per_core <- ceiling(runs / cores)
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  # Make set.seed() also affect seeds in clusters.
  parallel::clusterSetRNGStream(cl, sample.int(.Machine$integer.max, 1))
  parallel::clusterEvalQ(cl, library(slicer))

  worker <- function(
    i, D, training_idx, test_idx, y_train, y_test, centre, use_gradient,
    runs_per_core, kernels, params, maxit, squared
  ) {
    tryCatch(
      fit_gpr_one_core(D, training_idx, test_idx, y_train, y_test,
                       centre = centre, use_gradient = use_gradient,
                       runs = runs_per_core, verbose = FALSE,
                       squared = squared,
                       kernels = kernels, params = params, maxit = maxit),
      error = function(e) NULL
    )
  }
  environment(worker) <- environment(fit_gpr)

  results <- parallel::parLapply(
    cl, seq_len(cores), worker,
    D = D, training_idx = training_idx, test_idx = test_idx,
    y_train = y_train, y_test = y_test,
    centre = centre, use_gradient = use_gradient, runs_per_core = runs_per_core,
    kernels = kernels, params = params, maxit = maxit, squared = squared
  )
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0L) stop("All optimisation runs failed on all cores.")
  best_idx <- which.min(sapply(results, `[[`, "nll"))
  results[[best_idx]]
}
