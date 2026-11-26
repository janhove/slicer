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
  kernels = "rbf",
  params = NA,
  runs = 10L,
  cores = 1L,
  verbose = TRUE,
  maxit = 500
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
  Else, the optimiser uses the finite-differences method. Ignored when
  non-RBF kernels are used.

- kernels:

  A vector with kernels (`"rbf", "powexp", "matern"`).

- params:

  Kernel parameters (smoothness for Matern; exponent for power
  exponential; ignored for RBF).

- runs:

  Number of independent attempts to find a minimum when optimising the
  hyperparameters.

- cores:

  Number of cores used for parallel processing.

- verbose:

  If `TRUE`, progress is shown on the console. Only works if
  `cores == 1L`.

- maxit:

  Maximum number of iterations per optimisation run.

## Value

A list containing the predictions for the test objects as well as their
variance, the root mean squared error (if the true test outcomes are
provided), and the three tuned hyperparameter values.

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
#> Optimum set at -9.0169.
#> Hyperparameter search 2 of 50.
#> Hyperparameter search 3 of 50.
#> Hyperparameter search 4 of 50.
#> Hyperparameter search 5 of 50.
#> Current optimum improved from -9.0169 to -9.017.
#> Hyperparameter search 6 of 50.
#> Hyperparameter search 7 of 50.
#> Hyperparameter search 8 of 50.
#> Hyperparameter search 9 of 50.
#> Hyperparameter search 10 of 50.
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
#> Current optimum improved from -9.017 to -9.017.
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
plot(y_test, fit$test_predictions)
abline(a = 0, b = 1, lty = 1)
segments(x0 = y_test,
  y0 = fit$test_predictions - 2 * sqrt(diag(fit$test_variance)),
  y1 = fit$test_predictions + 2 * sqrt(diag(fit$test_variance)), lty = 2)

fit
#> $kernels
#> [1] "rbf" "rbf"
#> 
#> $params
#> [1] NA NA
#> 
#> $test_predictions
#>  [1]  0.103998074 -0.132056795 -0.003187666  0.001400324  0.005432131
#>  [6]  0.006513607 -0.231347151 -0.120350279 -0.010917716 -0.041973861
#> 
#> $test_variance
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,]  7.165663e-03 -3.035588e-04  1.834325e-03  0.0017402334  1.789813e-03
#>  [2,] -3.035588e-04  6.804513e-03  1.928300e-04 -0.0002674140 -1.831808e-03
#>  [3,]  1.834325e-03  1.928300e-04  4.650215e-03  0.0016076773  2.604366e-03
#>  [4,]  1.740233e-03 -2.674140e-04  1.607677e-03  0.0090087348 -1.679742e-04
#>  [5,]  1.789813e-03 -1.831808e-03  2.604366e-03 -0.0001679742  6.180415e-03
#>  [6,]  8.675107e-05  9.699278e-04  2.495755e-04  0.0002130526 -7.436115e-05
#>  [7,] -2.907428e-07 -1.743166e-03  3.553956e-04  0.0008373809  3.675886e-03
#>  [8,] -3.050445e-04  6.669875e-03  1.258966e-04 -0.0003503136 -1.754112e-03
#>  [9,] -1.677170e-03 -1.341922e-05 -3.755381e-04  0.0011672397 -1.242417e-03
#> [10,] -9.632011e-04  2.091705e-03  9.225876e-05  0.0019781729 -1.924133e-03
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,]  8.675107e-05 -2.907428e-07 -3.050445e-04 -1.677170e-03 -9.632011e-04
#>  [2,]  9.699278e-04 -1.743166e-03  6.669875e-03 -1.341922e-05  2.091705e-03
#>  [3,]  2.495755e-04  3.553956e-04  1.258966e-04 -3.755381e-04  9.225876e-05
#>  [4,]  2.130526e-04  8.373809e-04 -3.503136e-04  1.167240e-03  1.978173e-03
#>  [5,] -7.436115e-05  3.675886e-03 -1.754112e-03 -1.242417e-03 -1.924133e-03
#>  [6,]  5.688554e-03  1.778000e-03  1.052291e-03 -5.916770e-05  2.147873e-03
#>  [7,]  1.778000e-03  1.045605e-02 -1.449804e-03 -3.646119e-04  1.605238e-03
#>  [8,]  1.052291e-03 -1.449804e-03  6.667878e-03 -1.541095e-06  2.171421e-03
#>  [9,] -5.916770e-05 -3.646119e-04 -1.541095e-06  6.233065e-03  7.565016e-04
#> [10,]  2.147873e-03  1.605238e-03  2.171421e-03  7.565016e-04  1.022383e-02
#> 
#> $RMSE
#> [1] 0.1799872
#> 
#> $length_scale
#> [1] 1.24561832 0.05832573
#> 
#> $scaling_factor
#> [1] 0.03437302 0.01704833
#> 
#> $noise_variance
#> [1] 0.0100571
#> 
#> $nll
#> [1] -9.017025
#> 
# Multiple cores
fit <- fit_gpr(list(D2_1, D2_2), seq_len(N1), N1 + seq_len(N2),
 y_train, y_test, runs = 50, cores = 2L)
fit
#> $kernels
#> [1] "rbf" "rbf"
#> 
#> $params
#> [1] NA NA
#> 
#> $test_predictions
#>  [1]  0.103857248 -0.132056121 -0.003165938  0.001262344  0.005393760
#>  [6]  0.006479525 -0.231334167 -0.120325946 -0.010855628 -0.041950517
#> 
#> $test_variance
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,]  7.163489e-03 -3.062189e-04  1.831586e-03  0.0017326799  1.789466e-03
#>  [2,] -3.062189e-04  6.804347e-03  1.939820e-04 -0.0002638328 -1.832054e-03
#>  [3,]  1.831586e-03  1.939820e-04  4.648404e-03  0.0016057137  2.602978e-03
#>  [4,]  1.732680e-03 -2.638328e-04  1.605714e-03  0.0090110984 -1.658209e-04
#>  [5,]  1.789466e-03 -1.832054e-03  2.602978e-03 -0.0001658209  6.179096e-03
#>  [6,]  8.628396e-05  9.691977e-04  2.485787e-04  0.0002124405 -7.394114e-05
#>  [7,]  6.880622e-07 -1.740248e-03  3.486675e-04  0.0008387809  3.670280e-03
#>  [8,] -3.086953e-04  6.669714e-03  1.270534e-04 -0.0003463312 -1.754032e-03
#>  [9,] -1.675990e-03 -1.038831e-05 -3.754165e-04  0.0011677703 -1.241485e-03
#> [10,] -9.617710e-04  2.091902e-03  9.152251e-05  0.0019680415 -1.923263e-03
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,]  8.628396e-05  6.880622e-07 -0.0003086953 -1.675990e-03 -9.617710e-04
#>  [2,]  9.691977e-04 -1.740248e-03  0.0066697143 -1.038831e-05  2.091902e-03
#>  [3,]  2.485787e-04  3.486675e-04  0.0001270534 -3.754165e-04  9.152251e-05
#>  [4,]  2.124405e-04  8.387809e-04 -0.0003463312  1.167770e-03  1.968042e-03
#>  [5,] -7.394114e-05  3.670280e-03 -0.0017540317 -1.241485e-03 -1.923263e-03
#>  [6,]  5.683405e-03  1.776450e-03  0.0010514553 -5.789775e-05  2.145467e-03
#>  [7,]  1.776450e-03  1.045716e-02 -0.0014471735 -3.631984e-04  1.602203e-03
#>  [8,]  1.051455e-03 -1.447174e-03  0.0066681249  1.382900e-06  2.171104e-03
#>  [9,] -5.789775e-05 -3.631984e-04  0.0000013829  6.230243e-03  7.565067e-04
#> [10,]  2.145467e-03  1.602203e-03  0.0021711037  7.565067e-04  1.021377e-02
#> 
#> $RMSE
#> [1] 0.1799812
#> 
#> $length_scale
#> [1] 1.24557241 0.05820329
#> 
#> $scaling_factor
#> [1] 0.03428494 0.01700997
#> 
#> $noise_variance
#> [1] 0.0100493
#> 
#> $nll
#> [1] -9.017049
#> 

# Single kernel
N1 <- 40
N2 <- 10
x_train <- seq(-pi, pi, length.out = N1)
x_test  <- runif(N2, -pi, pi)
y_train <- x_train * plogis(x_train) * cos(x_train) + rnorm(N1, sd = 0.5)
y_test <- x_test * plogis(x_test) * cos(x_test)
D2 <- outer(c(x_train, x_test), c(x_train, x_test), "-")^2
fit <- fit_gpr(D2, seq_len(N1), N1 + seq_len(N2), y_train, y_test, runs = 50L)
#> Hyperparameter search 1 of 50.
#> Optimum set at 1062.563.
#> Hyperparameter search 2 of 50.
#> Current optimum improved from 1062.563 to 56.7493.
#> Hyperparameter search 3 of 50.
#> Current optimum improved from 56.7493 to 34.308.
#> Hyperparameter search 4 of 50.
#> Current optimum improved from 34.308 to 34.3079.
#> Hyperparameter search 5 of 50.
#> Hyperparameter search 6 of 50.
#> Hyperparameter search 7 of 50.
#> Hyperparameter search 8 of 50.
#> Hyperparameter search 9 of 50.
#> Hyperparameter search 10 of 50.
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
#> Current optimum improved from 34.3079 to 34.3079.
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
curve(x * plogis(x) * cos(x), -pi, pi,
  ylim = range(
    c(y_train, fit$test_predictions + 2 * sqrt(diag(fit$test_variance)),
               fit$test_predictions - 2 * sqrt(diag(fit$test_variance)))
   )
 )
points(x_train, y_train, pch = 1)
points(x_test, fit$test_predictions, pch = 16)
segments(x0 = x_test,
  y0 = fit$test_predictions - 2 * sqrt(diag(fit$test_variance)),
  y1 = fit$test_predictions + 2 * sqrt(diag(fit$test_variance)))

fit
#> $kernels
#> [1] "rbf"
#> 
#> $params
#> [1] NA
#> 
#> $test_predictions
#>  [1]  0.12674555  0.04401565  0.20231917 -2.46794471  0.16431597  0.10242354
#>  [7] -0.18537081 -0.24208339  0.06769748  0.16698463
#> 
#> $test_variance
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,]  0.0281285726  0.0177677974  0.0099659296  1.194184e-04  2.488318e-02
#>  [2,]  0.0177677974  0.0271209258 -0.0068130354  1.947958e-04  7.630755e-03
#>  [3,]  0.0099659296 -0.0068130354  0.0638297488 -1.127558e-04  2.912014e-02
#>  [4,]  0.0001194184  0.0001947958 -0.0001127558  3.185629e-02 -1.191066e-05
#>  [5,]  0.0248831808  0.0076307550  0.0291201446 -1.191066e-05  3.021271e-02
#>  [6,]  0.0269533237  0.0226592576  0.0016738649  1.896100e-04  1.991933e-02
#>  [7,] -0.0038297486  0.0069168022  0.0009365585 -5.542256e-04 -3.806256e-03
#>  [8,]  0.0008284878 -0.0031454584  0.0001192469  1.245963e-03  1.461837e-03
#>  [9,]  0.0221047950  0.0265981926 -0.0050795128  2.261793e-04  1.227209e-02
#> [10,]  0.0243655726  0.0068038556  0.0308316157 -2.091805e-05  3.043376e-02
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,]  0.0269533237 -0.0038297486  0.0008284878  0.0221047950  2.436557e-02
#>  [2,]  0.0226592576  0.0069168022 -0.0031454584  0.0265981926  6.803856e-03
#>  [3,]  0.0016738649  0.0009365585  0.0001192469 -0.0050795128  3.083162e-02
#>  [4,]  0.0001896100 -0.0005542256  0.0012459627  0.0002261793 -2.091805e-05
#>  [5,]  0.0199193253 -0.0038062562  0.0014618368  0.0122720907  3.043376e-02
#>  [6,]  0.0281237995 -0.0019736157 -0.0002016434  0.0259126387  1.914148e-02
#>  [7,] -0.0019736157  0.0262216551  0.0038514525  0.0028471325 -3.665881e-03
#>  [8,] -0.0002016434  0.0038514525  0.0260416934 -0.0020454192  1.451789e-03
#>  [9,]  0.0259126387  0.0028471325 -0.0020454192  0.0276329466  1.138556e-02
#> [10,]  0.0191414805 -0.0036658809  0.0014517887  0.0113855610  3.072530e-02
#> 
#> $RMSE
#> [1] 0.07643061
#> 
#> $length_scale
#> [1] 1.151216
#> 
#> $scaling_factor
#> [1] 1.463956
#> 
#> $noise_variance
#> [1] 0.1894824
#> 
#> $nll
#> [1] 34.3079
#> 
fit$RMSE
#> [1] 0.07643061
```
