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
  Else, the optimiser uses the finite-differences method. Ignored when
  non-RBF kernels are used.

- kernels:

  A vector with kernels (`"rbf", "matern05", "matern15", "matern25"`).

- runs:

  Number of independent attempts to find a minimum when optimising the
  hyperparameters.

- cores:

  Number of cores used for parallel processing.

- verbose:

  If `TRUE`, progress is shown on the console. Only works if
  `cores == 1L`.

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
#> Optimum set at -6.7313219.
#> Hyperparameter search 2 of 50.
#> Hyperparameter search 3 of 50.
#> Current optimum improved from -6.7313219 to -9.0170496.
#> Hyperparameter search 4 of 50.
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
#> $test_predictions
#>  [1]  0.103819650 -0.132069036 -0.003162020  0.001227841  0.005412415
#>  [6]  0.006466601 -0.231317504 -0.120333455 -0.010845450 -0.041972268
#> 
#> $test_variance
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,]  7.163161e-03 -3.074587e-04  0.0018305542  0.0017310743  1.789838e-03
#>  [2,] -3.074587e-04  6.804186e-03  0.0001942065 -0.0002627083 -1.832662e-03
#>  [3,]  1.830554e-03  1.942065e-04  0.0046481272  0.0016052214  2.603473e-03
#>  [4,]  1.731074e-03 -2.627083e-04  0.0016052214  0.0090127746 -1.649392e-04
#>  [5,]  1.789838e-03 -1.832662e-03  0.0026034731 -0.0001649392  6.179768e-03
#>  [6,]  8.597651e-05  9.687542e-04  0.0002482526  0.0002122387 -7.382978e-05
#>  [7,]  6.731911e-07 -1.739089e-03  0.0003465327  0.0008385609  3.669843e-03
#>  [8,] -3.101901e-04  6.669570e-03  0.0001272532 -0.0003451489 -1.754525e-03
#>  [9,] -1.675449e-03 -9.226206e-06 -0.0003751851  0.0011681249 -1.241379e-03
#> [10,] -9.622552e-04  2.091820e-03  0.0000906921  0.0019653028 -1.923098e-03
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,]  8.597651e-05  6.731911e-07 -3.101901e-04 -1.675449e-03 -0.0009622552
#>  [2,]  9.687542e-04 -1.739089e-03  6.669570e-03 -9.226206e-06  0.0020918200
#>  [3,]  2.482526e-04  3.465327e-04  1.272532e-04 -3.751851e-04  0.0000906921
#>  [4,]  2.122387e-04  8.385609e-04 -3.451489e-04  1.168125e-03  0.0019653028
#>  [5,] -7.382978e-05  3.669843e-03 -1.754525e-03 -1.241379e-03 -0.0019230985
#>  [6,]  5.682890e-03  1.776530e-03  1.051009e-03 -5.747128e-05  0.0021449764
#>  [7,]  1.776530e-03  1.045821e-02 -1.446089e-03 -3.629493e-04  0.0016013236
#>  [8,]  1.051009e-03 -1.446089e-03  6.668123e-03  2.502042e-06  0.0021708825
#>  [9,] -5.747128e-05 -3.629493e-04  2.502042e-06  6.230363e-03  0.0007568228
#> [10,]  2.144976e-03  1.601324e-03  2.170883e-03  7.568228e-04  0.0102118042
#> 
#> $RMSE
#> [1] 0.1799737
#> 
#> $length_scale
#> [1] 1.24597125 0.05817039
#> 
#> $scaling_factor
#> [1] 0.03429806 0.01700076
#> 
#> $noise_variance
#> [1] 0.01004947
#> 
#> $nll
#> [1] -9.01705
#> 
# Multiple cores
fit <- fit_gpr(list(D2_1, D2_2), seq_len(N1), N1 + seq_len(N2),
  y_train, y_test, runs = 50, cores = 2)
fit
#> $kernels
#> [1] "rbf" "rbf"
#> 
#> $test_predictions
#>  [1]  0.103819163 -0.132069190 -0.003161879  0.001227344  0.005412343
#>  [6]  0.006466821 -0.231317513 -0.120333493 -0.010845136 -0.041972108
#> 
#> $test_variance
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,]  7.163170e-03 -3.074687e-04  1.830548e-03  0.0017310460  1.789841e-03
#>  [2,] -3.074687e-04  6.804206e-03  1.942149e-04 -0.0002626940 -1.832670e-03
#>  [3,]  1.830548e-03  1.942149e-04  4.648129e-03  0.0016052188  2.603469e-03
#>  [4,]  1.731046e-03 -2.626940e-04  1.605219e-03  0.0090128093 -1.649325e-04
#>  [5,]  1.789841e-03 -1.832670e-03  2.603469e-03 -0.0001649325  6.179777e-03
#>  [6,]  8.597537e-05  9.687521e-04  2.482507e-04  0.0002122388 -7.382774e-05
#>  [7,]  6.773605e-07 -1.739091e-03  3.465054e-04  0.0008385733  3.669831e-03
#>  [8,] -3.102044e-04  6.669590e-03  1.272612e-04 -0.0003451333 -1.754531e-03
#>  [9,] -1.675453e-03 -9.214794e-06 -3.751874e-04  0.0011681287 -1.241381e-03
#> [10,] -9.622526e-04  2.091828e-03  9.069209e-05  0.0019652690 -1.923104e-03
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,]  8.597537e-05  6.773605e-07 -3.102044e-04 -1.675453e-03 -9.622526e-04
#>  [2,]  9.687521e-04 -1.739091e-03  6.669590e-03 -9.214794e-06  2.091828e-03
#>  [3,]  2.482507e-04  3.465054e-04  1.272612e-04 -3.751874e-04  9.069209e-05
#>  [4,]  2.122388e-04  8.385733e-04 -3.451333e-04  1.168129e-03  1.965269e-03
#>  [5,] -7.382774e-05  3.669831e-03 -1.754531e-03 -1.241381e-03 -1.923104e-03
#>  [6,]  5.682882e-03  1.776530e-03  1.051007e-03 -5.746600e-05  2.144973e-03
#>  [7,]  1.776530e-03  1.045825e-02 -1.446091e-03 -3.629448e-04  1.601313e-03
#>  [8,]  1.051007e-03 -1.446091e-03  6.668144e-03  2.513272e-06  2.170889e-03
#>  [9,] -5.746600e-05 -3.629448e-04  2.513272e-06  6.230367e-03  7.568269e-04
#> [10,]  2.144973e-03  1.601313e-03  2.170889e-03  7.568269e-04  1.021179e-02
#> 
#> $RMSE
#> [1] 0.1799737
#> 
#> $length_scale
#> [1] 1.24597238 0.05816994
#> 
#> $scaling_factor
#> [1] 0.03429813 0.01700070
#> 
#> $noise_variance
#> [1] 0.01004944
#> 
#> $nll
#> [1] -9.01705
#> 

# Single kernel
N1 <- 40
N2 <- 10
x_train <- seq(-pi, pi, length.out = N1)
x_test  <- runif(N2, -pi, pi)
y_train <- x_train * plogis(x_train) * cos(x_train) + rnorm(N1, sd = 0.5)
y_test <- x_test * plogis(x_test) * cos(x_test)
D2 <- outer(c(x_train, x_test), c(x_train, x_test), "-")^2
fit <- fit_gpr(D2, seq_len(N1), N1 + seq_len(N2), y_train, y_test,
  runs = 50L, cores = 2)
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
#> $test_predictions
#>  [1]  0.12674565  0.04401594  0.20231944 -2.46794517  0.16431599  0.10242373
#>  [7] -0.18537117 -0.24208302  0.06769776  0.16698464
#> 
#> $test_variance
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,]  0.0281285179  0.0177677796  0.0099659146  1.194223e-04  2.488314e-02
#>  [2,]  0.0177677796  0.0271208674 -0.0068130371  1.947949e-04  7.630760e-03
#>  [3,]  0.0099659146 -0.0068130371  0.0638297220 -1.127583e-04  2.912011e-02
#>  [4,]  0.0001194223  0.0001947949 -0.0001127583  3.185623e-02 -1.190842e-05
#>  [5,]  0.0248831360  0.0076307603  0.0291201111 -1.190842e-05  3.021265e-02
#>  [6,]  0.0269532734  0.0226592221  0.0016738594  1.896135e-04  1.991930e-02
#>  [7,] -0.0038297364  0.0069168135  0.0009365492 -5.542289e-04 -3.806260e-03
#>  [8,]  0.0008284749 -0.0031454711  0.0001192588  1.245971e-03  1.461838e-03
#>  [9,]  0.0221047622  0.0265981387 -0.0050795132  2.261805e-04  1.227208e-02
#> [10,]  0.0243655293  0.0068038619  0.0308315811 -2.091602e-05  3.043371e-02
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,]  0.0269532734 -0.0038297364  0.0008284749  0.0221047622  2.436553e-02
#>  [2,]  0.0226592221  0.0069168135 -0.0031454711  0.0265981387  6.803862e-03
#>  [3,]  0.0016738594  0.0009365492  0.0001192588 -0.0050795132  3.083158e-02
#>  [4,]  0.0001896135 -0.0005542289  0.0012459706  0.0002261805 -2.091602e-05
#>  [5,]  0.0199192951 -0.0038062603  0.0014618376  0.0122720836  3.043371e-02
#>  [6,]  0.0281237448 -0.0019735963 -0.0002016619  0.0259125915  1.914145e-02
#>  [7,] -0.0019735963  0.0262215907  0.0038514747  0.0028471516 -3.665886e-03
#>  [8,] -0.0002016619  0.0038514747  0.0260416271 -0.0020454373  1.451791e-03
#>  [9,]  0.0259125915  0.0028471516 -0.0020454373  0.0276328902  1.138556e-02
#> [10,]  0.0191414522 -0.0036658860  0.0014517906  0.0113855555  3.072525e-02
#> 
#> $RMSE
#> [1] 0.07643054
#> 
#> $length_scale
#> [1] 1.151219
#> 
#> $scaling_factor
#> [1] 1.463973
#> 
#> $noise_variance
#> [1] 0.1894821
#> 
#> $nll
#> [1] 34.3079
#> 
fit$RMSE
#> [1] 0.07643054

fit <- fit_gpr(D2, seq_len(N1), N1 + seq_len(N2), y_train, y_test,
  runs = 50L, cores = 2, kernel = "matern15")
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
#> [1] "matern15"
#> 
#> $test_predictions
#>  [1]  0.11465887  0.04580080  0.22597408 -2.42672817  0.15906636  0.10122310
#>  [7] -0.14341704 -0.22417308  0.07728402  0.16326950
#> 
#> $test_variance
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,]  3.601562e-02  1.263452e-02  6.007389e-03  6.105070e-08  2.638459e-02
#>  [2,]  1.263452e-02  3.578988e-02 -2.890401e-03 -3.159343e-07  3.201891e-03
#>  [3,]  6.007389e-03 -2.890401e-03  6.863614e-02  2.188487e-08  2.472915e-02
#>  [4,]  6.105070e-08 -3.159343e-07  2.188487e-08  3.791404e-02  5.889949e-08
#>  [5,]  2.638459e-02  3.201891e-03  2.472915e-02  5.889949e-08  3.693619e-02
#>  [6,]  3.172168e-02  2.093001e-02  3.442747e-04  2.311123e-08  1.724662e-02
#>  [7,] -1.281776e-03  2.137100e-03 -1.297633e-04 -1.440008e-06 -8.372175e-04
#>  [8,]  9.686162e-06 -4.911868e-04  3.581740e-05  4.254383e-05  4.852542e-05
#>  [9,]  1.976326e-02  3.250506e-02 -2.653576e-03 -1.303920e-07  7.270835e-03
#> [10,]  2.523322e-02  2.693544e-03  2.671326e-02  5.719230e-08  3.698641e-02
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,]  3.172168e-02 -1.281776e-03  9.686162e-06  1.976326e-02  2.523322e-02
#>  [2,]  2.093001e-02  2.137100e-03 -4.911868e-04  3.250506e-02  2.693544e-03
#>  [3,]  3.442747e-04 -1.297633e-04  3.581740e-05 -2.653576e-03  2.671326e-02
#>  [4,]  2.311123e-08 -1.440008e-06  4.254383e-05 -1.303920e-07  5.719230e-08
#>  [5,]  1.724662e-02 -8.372175e-04  4.852542e-05  7.270835e-03  3.698641e-02
#>  [6,]  3.600664e-02 -1.202788e-03 -6.672627e-05  2.867467e-02  1.612914e-02
#>  [7,] -1.202788e-03  3.567173e-02  3.050395e-04  9.394330e-05 -7.944846e-04
#>  [8,] -6.672627e-05  3.050395e-04  3.567342e-02 -2.800326e-04  4.893046e-05
#>  [9,]  2.867467e-02  9.394330e-05 -2.800326e-04  3.587027e-02  6.514692e-03
#> [10,]  1.612914e-02 -7.944846e-04  4.893046e-05  6.514692e-03  3.722977e-02
#> 
#> $RMSE
#> [1] 0.07128682
#> 
#> $length_scale
#> [1] 2.17199
#> 
#> $scaling_factor
#> [1] 2.03781
#> 
#> $noise_variance
#> [1] 0.1888017
#> 
#> $nll
#> [1] 35.37379
#> 
fit$RMSE
#> [1] 0.07128682
```
