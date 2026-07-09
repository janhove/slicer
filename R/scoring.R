#' Negative Log Predictive Density for Gaussian Process Regression Models
#'
#' Computes the average negative log predictive density (NLPD) of test outcomes
#' under the posterior Gaussian distribution returned by a GPR model. If the
#' posterior covariance is not positive definite, a small jitter term is added
#' to the main diagonal for numerical stability; this jitter term is increased
#' adaptively as necessary.
#'
#' @param fit A Gaussian process model fit obtained by [gpr_predict()] or
#'    [fit_gpr()].
#' @param y_test A vector of the true test outcomes.
#' @param add_noise If `TRUE` (default), `y_test` is assumed to be noisy
#'    and the observation noise is added to the variance in the posterior
#'    covariance.
#' @param nugget Numeric; base jitter added to the covariance diagonal if needed
#'   to ensure positive definiteness. Increased multiplicatively if required.
#'
#' @return A scalar giving the average NLPD. The returned value has an attribute
#'   `"jitter_used"` indicating the final diagonal jitter added (0 if none).
#'
#' @export
#'
#' @examples
#' set.seed(2026-06-25)
#' N1 <- 50
#' N2 <- 20
#' x_train <- seq(-pi, pi, length.out = N1)
#' x_test  <- runif(N2, -pi, pi)
#' y_train <- x_train * plogis(x_train) * cos(x_train) + rnorm(N1, sd = 0.5)
#' y_test <- x_test * plogis(x_test) * cos(x_test)
#' D2 <- outer(c(x_train, x_test), c(x_train, x_test), "-")^2
#' fit <- fit_gpr(D2, seq_len(N1), N1 + seq_len(N2), y_train, y_test,
#'                runs = 20L, cores = 1)
#' # noisy outcome
#' nlpd_gpr(fit, y_test + rnorm(N2, sd = 0.5))
#' # clean outcome
#' nlpd_gpr(fit, y_test, add_noise = FALSE)
nlpd_gpr <- function(fit, y_test, add_noise = TRUE, nugget = 1e-10) {
  Sigma <- fit$test_variance
  if (is.vector(Sigma)) {
    Sigma <- as.matrix(Sigma)
  }
  if (any(is.na(Sigma))) {
    warning("The fit's test variance contains NAs, likely because the distances among the test inputs weren't provided.\nReplacing NAs by 0.")
    Sigma[is.na(Sigma)] <- 0
  }
  if (add_noise) {
    diag(Sigma) <- diag(Sigma) + fit$noise_variance
  }

  I <- diag(nrow(Sigma))
  Sigma_base <- Sigma

  pd_check <- tryCatch(
   chol(Sigma_base),
    error = function(e) {
      warning("Covariance matrix not positive definite. Trying jitter escalation.")

      for (multiplier in c(1, 10, 100, 1000, 10000)) {
        Sigma_try <- Sigma_base + (nugget * multiplier) * I

        result <- tryCatch(
          chol(Sigma_try),
          error = function(e) NULL
        )

        if (!is.null(result)) {
          message(sprintf("Cholesky succeeded with jitter = %g.", nugget * multiplier))
          Sigma <<- Sigma_try
          attr(result, "jitter_used") <- nugget * multiplier
          return(result)
        }
      }

      stop("Sigma not positive definite even after jitter escalation.")
    }
  )
  if (is.null(attr(pd_check, "jitter_used"))) {
    attr(pd_check, "jitter_used") <- 0
  }

  means <- fit$test_predictions

  log_density <- mvtnorm::dmvnorm(y_test, means, Sigma, log = TRUE)
  nlpd <- -log_density / length(y_test)
  attr(nlpd, "jitter_used") <- attr(pd_check, "jitter_used")
  nlpd
}

