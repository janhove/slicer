# Gaussian Process Fitting to Test Data with Tuned Hyperparameters

This function generates predictions using a Gaussian process model with
one or several Gaussian RBF kernels. The hyperparameters are tuned using
the training data by minimising the model's negative marginal
log-likelihood.

## Usage

``` r
fit_gpr(
  D2,
  training_idx,
  test_idx,
  y_train,
  y_test = NULL,
  centre = TRUE,
  use_gradient = TRUE,
  runs = 10L,
  cores = 1L,
  verbose = TRUE
)
```

## Arguments

- D2:

  A matrix, or a list of matrices, containing squared pairwise distances
  between all objects.

- training_idx:

  A vector with the row (and column) numbers corresponding to the
  training entries in `D2`.

- test_idx:

  A vector with the row (and column) numbers corresponding to the test
  entries in `D2`.

- y_train:

  A vector with the training outcomes.

- y_test:

  Optionally, a vector with the test outcomes.

- centre:

  If `TRUE`, the training outcomes are centred around their mean before
  fitting the model. The mean is then added back to the predictions at
  the end.

- use_gradient:

  If `TRUE`, closed-form expressions for the RBF's gradient are used.
  Else, the optimiser uses the finite-differences method.

- runs:

  Number of independent attempts to find a minimum when optimising the
  hyperparameters.

- cores:

  Number of cores used for parallel processing.

- verbose:

  If `TRUE`, progress is shown on the console. Only works if
  `cores == 1L`.

## Value

A list containing the predictions for the test objects, the root mean
squared error (if the true test outcomes are provided), the tuned
hyperparameter values, and the negative log marginal likelihood (nll).

## Examples

``` r
# Multiple kernels
N1 <- 25
N2 <- 10
x_train1 <- runif(N1, -pi, pi)
x_train2 <- runif(N1, 0, 1)
x_test1  <- runif(N2, -pi, pi)
x_test2  <- runif(N2, 0, 1)
y_train <- x_train2 * plogis(x_train1) * cos(x_train1)
y_test <- x_test2 * plogis(x_test1) * cos(x_test1)
D2_1 <- outer(c(x_train1, x_test1), c(x_train1, x_test1), "-")^2
D2_2 <- outer(c(x_train2, x_test2), c(x_train2, x_test2), "-")^2
# Single core
fit <- fit_gpr(list(D2_1, D2_2), seq_len(N1), N1 + seq_len(N2),
  y_train, y_test, runs = 50)
#> Hyperparameter search 1 of 50.
#> Optimum set at -10.8843568.
#> Hyperparameter search 2 of 50.
#> Hyperparameter search 3 of 50.
#> Hyperparameter search 4 of 50.
#> Hyperparameter search 5 of 50.
#> Hyperparameter search 6 of 50.
#> Hyperparameter search 7 of 50.
#> Hyperparameter search 8 of 50.
#> Hyperparameter search 9 of 50.
#> Hyperparameter search 10 of 50.
#> Current optimum improved from -10.8843568 to -10.8843568.
#> Hyperparameter search 11 of 50.
#> Hyperparameter search 12 of 50.
#> Hyperparameter search 13 of 50.
#> Hyperparameter search 14 of 50.
#> Hyperparameter search 15 of 50.
#> Hyperparameter search 16 of 50.
#> Hyperparameter search 17 of 50.
#> Hyperparameter search 18 of 50.
#> Hyperparameter search 19 of 50.
#> Hyperparameter search 20 of 50.
#> Hyperparameter search 21 of 50.
#> Hyperparameter search 22 of 50.
#> Hyperparameter search 23 of 50.
#> Hyperparameter search 24 of 50.
#> Hyperparameter search 25 of 50.
#> Hyperparameter search 26 of 50.
#> Hyperparameter search 27 of 50.
#> Hyperparameter search 28 of 50.
#> Hyperparameter search 29 of 50.
#> Hyperparameter search 30 of 50.
#> Hyperparameter search 31 of 50.
#> Hyperparameter search 32 of 50.
#> Hyperparameter search 33 of 50.
#> Hyperparameter search 34 of 50.
#> Hyperparameter search 35 of 50.
#> Hyperparameter search 36 of 50.
#> Hyperparameter search 37 of 50.
#> Hyperparameter search 38 of 50.
#> Hyperparameter search 39 of 50.
#> Hyperparameter search 40 of 50.
#> Hyperparameter search 41 of 50.
#> Hyperparameter search 42 of 50.
#> Hyperparameter search 43 of 50.
#> Hyperparameter search 44 of 50.
#> Hyperparameter search 45 of 50.
#> Hyperparameter search 46 of 50.
#> Hyperparameter search 47 of 50.
#> Hyperparameter search 48 of 50.
#> Hyperparameter search 49 of 50.
#> Hyperparameter search 50 of 50.
plot(fit$test_predictions, y_test)

fit
#> $test_predictions
#>  [1]  0.193023922  0.197859298 -0.281270778 -0.001806574 -0.050603571
#>  [6]  0.200220328  0.218295659  0.151833107  0.086117642 -0.110344091
#> 
#> $RMSE
#> [1] 0.06029104
#> 
#> $length_scale
#> length_scale1 length_scale2 
#>     1.6227648     0.7937291 
#> 
#> $variance
#>   variance1   variance2 
#> 0.055104948 0.005048196 
#> 
#> $lambda2
#>    lambda2 
#> 0.01457336 
#> 
#> $nll
#> [1] -10.88436
#> 
# Multiple cores
fit <- fit_gpr(list(D2_1, D2_2), seq_len(N1), N1 + seq_len(N2),
  y_train, y_test, runs = 50, cores = 2)
fit
#> $test_predictions
#>  [1]  0.193023866  0.197859317 -0.281270694 -0.001806513 -0.050603628
#>  [6]  0.200220282  0.218295666  0.151833111  0.086117670 -0.110344115
#> 
#> $RMSE
#> [1] 0.06029106
#> 
#> $length_scale
#> length_scale1 length_scale2 
#>     1.6227640     0.7937265 
#> 
#> $variance
#>   variance1   variance2 
#> 0.055104818 0.005048195 
#> 
#> $lambda2
#>    lambda2 
#> 0.01457336 
#> 
#> $nll
#> [1] -10.88436
#> 

# Single kernel
N1 <- 25
N2 <- 10
x_train <- seq(-pi, pi, length.out = N1)
x_test  <- runif(N2, -pi, pi)
y_train <- x_train * plogis(x_train) * cos(x_train)
y_test <- x_test * plogis(x_test) * cos(x_test)
D2 <- outer(c(x_train, x_test), c(x_train, x_test), "-")^2
fit <- fit_gpr(D2, seq_len(N1), N1 + seq_len(N2), y_train, y_test,
  runs = 50L, cores = 2)
curve(x * plogis(x) * cos(x), -pi, pi)
points(x_train, y_train, pch = 1)
points(x_test, fit$test_predictions, pch = 16)

fit$RMSE
#> [1] 1.818487e-07
```
