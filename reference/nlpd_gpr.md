# Negative log predictive density for Gaussian process regression models

Computes the average negative log predictive density (NLPD) of test
outcomes under the posterior Gaussian distribution returned by a GPR
model. If the posterior covariance is not positive definite, a small
diagonal "nugget" is added adaptively for numerical stability.

## Usage

``` r
nlpd_gpr(fit, y_test, add_noise = TRUE, nugget = 1e-10)
```

## Arguments

- fit:

  A Gaussian process model fit obtained by
  [`gpr_predict()`](https://janhove.github.io/slicer/reference/gpr_predict.md)
  or
  [`fit_gpr()`](https://janhove.github.io/slicer/reference/fit_gpr.md).

- y_test:

  A vector of the true test outcomes.

- add_noise:

  If `TRUE` (default), `y_test` is assumed to be noisy and the
  observation noise is added to the variance in the posterior
  covariance.

- nugget:

  Numeric; base jitter added to the covariance diagonal if needed to
  ensure positive definiteness. Increased multiplicatively if required.

## Value

A scalar giving the average NLPD. The returned value has an attribute
`"jitter_used"` indicating the final diagonal jitter added (0 if none).

## Examples

``` r
set.seed(2026-06-25)
N1 <- 50
N2 <- 20
x_train <- seq(-pi, pi, length.out = N1)
x_test  <- runif(N2, -pi, pi)
y_train <- x_train * plogis(x_train) * cos(x_train) + rnorm(N1, sd = 0.5)
y_test <- x_test * plogis(x_test) * cos(x_test)
D2 <- outer(c(x_train, x_test), c(x_train, x_test), "-")^2
fit <- fit_gpr(D2, seq_len(N1), N1 + seq_len(N2), y_train, y_test,
               runs = 20L, cores = 1)
#> Hyperparameter search 1 of 20.
#> Optimum set at 65.1048849.
#> Hyperparameter search 2 of 20.
#> Current optimum improved from 65.1048849 to 57.4472943.
#> Hyperparameter search 3 of 20.
#> Hyperparameter search 4 of 20.
#> Hyperparameter search 5 of 20.
#> Hyperparameter search 6 of 20.
#> Hyperparameter search 7 of 20.
#> Hyperparameter search 8 of 20.
#> Hyperparameter search 9 of 20.
#> Hyperparameter search 10 of 20.
#> Hyperparameter search 11 of 20.
#> Current optimum improved from 57.4472943 to 40.4504072.
#> Hyperparameter search 12 of 20.
#> Hyperparameter search 13 of 20.
#> Hyperparameter search 14 of 20.
#> Hyperparameter search 15 of 20.
#> Hyperparameter search 16 of 20.
#> Hyperparameter search 17 of 20.
#> Hyperparameter search 18 of 20.
#> Hyperparameter search 19 of 20.
#> Hyperparameter search 20 of 20.
# noisy outcome
nlpd_gpr(fit, y_test + rnorm(N2, sd = 0.5))
#> [1] 1.130341
#> attr(,"jitter_used")
#> [1] 0
# clean outcome
nlpd_gpr(fit, y_test, add_noise = FALSE)
#> Warning: Covariance matrix not positive definite. Trying jitter escalation.
#> Cholesky succeeded with jitter = 1e-10.
#> [1] -3.758439
#> attr(,"jitter_used")
#> [1] 1e-10
```
