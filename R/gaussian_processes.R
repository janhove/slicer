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
