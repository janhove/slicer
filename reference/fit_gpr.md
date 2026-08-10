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
#>  [7,]  6.731916e-07 -1.739089e-03  0.0003465327  0.0008385609  3.669843e-03
#>  [8,] -3.101901e-04  6.669570e-03  0.0001272532 -0.0003451489 -1.754525e-03
#>  [9,] -1.675449e-03 -9.226205e-06 -0.0003751851  0.0011681249 -1.241379e-03
#> [10,] -9.622552e-04  2.091820e-03  0.0000906921  0.0019653028 -1.923098e-03
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,]  8.597651e-05  6.731916e-07 -3.101901e-04 -1.675449e-03 -0.0009622552
#>  [2,]  9.687542e-04 -1.739089e-03  6.669570e-03 -9.226205e-06  0.0020918200
#>  [3,]  2.482526e-04  3.465327e-04  1.272532e-04 -3.751851e-04  0.0000906921
#>  [4,]  2.122387e-04  8.385609e-04 -3.451489e-04  1.168125e-03  0.0019653028
#>  [5,] -7.382978e-05  3.669843e-03 -1.754525e-03 -1.241379e-03 -0.0019230985
#>  [6,]  5.682890e-03  1.776530e-03  1.051009e-03 -5.747128e-05  0.0021449764
#>  [7,]  1.776530e-03  1.045821e-02 -1.446089e-03 -3.629493e-04  0.0016013236
#>  [8,]  1.051009e-03 -1.446089e-03  6.668123e-03  2.502043e-06  0.0021708825
#>  [9,] -5.747128e-05 -3.629493e-04  2.502043e-06  6.230363e-03  0.0007568228
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
#>  [1]  0.103819619 -0.132068978 -0.003162014  0.001227821  0.005412202
#>  [6]  0.006466820 -0.231317773 -0.120333386 -0.010845338 -0.041972070
#> 
#> $test_variance
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,]  7.163167e-03 -3.074578e-04  1.830556e-03  0.0017310721  1.789838e-03
#>  [2,] -3.074578e-04  6.804196e-03  1.942081e-04 -0.0002627090 -1.832663e-03
#>  [3,]  1.830556e-03  1.942081e-04  4.648129e-03  0.0016052226  2.603470e-03
#>  [4,]  1.731072e-03 -2.627090e-04  1.605223e-03  0.0090127827 -1.649407e-04
#>  [5,]  1.789838e-03 -1.832663e-03  2.603470e-03 -0.0001649407  6.179770e-03
#>  [6,]  8.597742e-05  9.687541e-04  2.482534e-04  0.0002122396 -7.382909e-05
#>  [7,]  6.758993e-07 -1.739098e-03  3.465333e-04  0.0008385667  3.669844e-03
#>  [8,] -3.101895e-04  6.669580e-03  1.272548e-04 -0.0003451494 -1.754526e-03
#>  [9,] -1.675454e-03 -9.227159e-06 -3.751870e-04  0.0011681243 -1.241381e-03
#> [10,] -9.622534e-04  2.091824e-03  9.069556e-05  0.0019653044 -1.923102e-03
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,]  8.597742e-05  6.758993e-07 -3.101895e-04 -1.675454e-03 -9.622534e-04
#>  [2,]  9.687541e-04 -1.739098e-03  6.669580e-03 -9.227159e-06  2.091824e-03
#>  [3,]  2.482534e-04  3.465333e-04  1.272548e-04 -3.751870e-04  9.069556e-05
#>  [4,]  2.122396e-04  8.385667e-04 -3.451494e-04  1.168124e-03  1.965304e-03
#>  [5,] -7.382909e-05  3.669844e-03 -1.754526e-03 -1.241381e-03 -1.923102e-03
#>  [6,]  5.682891e-03  1.776530e-03  1.051009e-03 -5.747144e-05  2.144977e-03
#>  [7,]  1.776530e-03  1.045823e-02 -1.446097e-03 -3.629505e-04  1.601322e-03
#>  [8,]  1.051009e-03 -1.446097e-03  6.668132e-03  2.501218e-06  2.170886e-03
#>  [9,] -5.747144e-05 -3.629505e-04  2.501218e-06  6.230367e-03  7.568235e-04
#> [10,]  2.144977e-03  1.601322e-03  2.170886e-03  7.568235e-04  1.021181e-02
#> 
#> $RMSE
#> [1] 0.1799738
#> 
#> $length_scale
#> [1] 1.24596979 0.05817041
#> 
#> $scaling_factor
#> [1] 0.03429800 0.01700081
#> 
#> $noise_variance
#> [1] 0.01004946
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
#>  [1]  0.12674552  0.04401568  0.20231940 -2.46794511  0.16431599  0.10242351
#>  [7] -0.18537077 -0.24208347  0.06769748  0.16698465
#> 
#> $test_variance
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,]  0.0281285620  0.0177677879  0.0099659022  1.194187e-04  0.0248831674
#>  [2,]  0.0177677879  0.0271209101 -0.0068130424  1.947965e-04  0.0076307483
#>  [3,]  0.0099659022 -0.0068130424  0.0638297547 -1.127562e-04  0.0291201185
#>  [4,]  0.0001194187  0.0001947965 -0.0001127562  3.185627e-02 -0.0000119108
#>  [5,]  0.0248831674  0.0076307483  0.0291201185 -1.191080e-05  0.0302126898
#>  [6,]  0.0269533135  0.0226592456  0.0016738422  1.896106e-04  0.0199193141
#>  [7,] -0.0038297517  0.0069167953  0.0009365659 -5.542267e-04 -0.0038062568
#>  [8,]  0.0008284902 -0.0031454580  0.0001192440  1.245964e-03  0.0014618386
#>  [9,]  0.0221047850  0.0265981779 -0.0050795261  2.261800e-04  0.0122720823
#> [10,]  0.0243655588  0.0068038490  0.0308315904 -2.091822e-05  0.0304337462
#>                [,6]          [,7]          [,8]         [,9]         [,10]
#>  [1,]  0.0269533135 -0.0038297517  0.0008284902  0.022104785  2.436556e-02
#>  [2,]  0.0226592456  0.0069167953 -0.0031454580  0.026598178  6.803849e-03
#>  [3,]  0.0016738422  0.0009365659  0.0001192440 -0.005079526  3.083159e-02
#>  [4,]  0.0001896106 -0.0005542267  0.0012459638  0.000226180 -2.091822e-05
#>  [5,]  0.0199193141 -0.0038062568  0.0014618386  0.012272082  3.043375e-02
#>  [6,]  0.0281237886 -0.0019736198 -0.0002016415  0.025912627  1.914147e-02
#>  [7,] -0.0019736198  0.0262216388  0.0038514481  0.002847127 -3.665881e-03
#>  [8,] -0.0002016415  0.0038514481  0.0260416764 -0.002045418  1.451790e-03
#>  [9,]  0.0259126268  0.0028471269 -0.0020454182  0.027632933  1.138555e-02
#> [10,]  0.0191414690 -0.0036658811  0.0014517904  0.011385553  3.072528e-02
#> 
#> $RMSE
#> [1] 0.07643068
#> 
#> $length_scale
#> [1] 1.151217
#> 
#> $scaling_factor
#> [1] 1.463965
#> 
#> $noise_variance
#> [1] 0.1894823
#> 
#> $nll
#> [1] 34.3079
#> 
fit$RMSE
#> [1] 0.07643068

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
#>  [1]  0.11465893  0.04580072  0.22597410 -2.42672806  0.15906642  0.10122310
#>  [7] -0.14341717 -0.22417297  0.07728394  0.16326955
#> 
#> $test_variance
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,]  3.601558e-02  1.263454e-02  6.007413e-03  6.105020e-08  2.638459e-02
#>  [2,]  1.263454e-02  3.578984e-02 -2.890403e-03 -3.159425e-07  3.201906e-03
#>  [3,]  6.007413e-03 -2.890403e-03  6.863611e-02  2.188564e-08  2.472917e-02
#>  [4,]  6.105020e-08 -3.159425e-07  2.188564e-08  3.791401e-02  5.890009e-08
#>  [5,]  2.638459e-02  3.201906e-03  2.472917e-02  5.890009e-08  3.693616e-02
#>  [6,]  3.172165e-02  2.093002e-02  3.442904e-04  2.310923e-08  1.724662e-02
#>  [7,] -1.281776e-03  2.137116e-03 -1.297661e-04 -1.439988e-06 -8.372208e-04
#>  [8,]  9.685372e-06 -4.911900e-04  3.581784e-05  4.254404e-05  4.852542e-05
#>  [9,]  1.976326e-02  3.250504e-02 -2.653573e-03 -1.303973e-07  7.270851e-03
#> [10,]  2.523321e-02  2.693558e-03  2.671328e-02  5.719293e-08  3.698638e-02
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,]  3.172165e-02 -1.281776e-03  9.685372e-06  1.976326e-02  2.523321e-02
#>  [2,]  2.093002e-02  2.137116e-03 -4.911900e-04  3.250504e-02  2.693558e-03
#>  [3,]  3.442904e-04 -1.297661e-04  3.581784e-05 -2.653573e-03  2.671328e-02
#>  [4,]  2.310923e-08 -1.439988e-06  4.254404e-05 -1.303973e-07  5.719293e-08
#>  [5,]  1.724662e-02 -8.372208e-04  4.852542e-05  7.270851e-03  3.698638e-02
#>  [6,]  3.600661e-02 -1.202785e-03 -6.672778e-05  2.867466e-02  1.612915e-02
#>  [7,] -1.202785e-03  3.567169e-02  3.050510e-04  9.395413e-05 -7.944880e-04
#>  [8,] -6.672778e-05  3.050510e-04  3.567338e-02 -2.800352e-04  4.893049e-05
#>  [9,]  2.867466e-02  9.395413e-05 -2.800352e-04  3.587024e-02  6.514708e-03
#> [10,]  1.612915e-02 -7.944880e-04  4.893049e-05  6.514708e-03  3.722974e-02
#> 
#> $RMSE
#> [1] 0.0712868
#> 
#> $length_scale
#> [1] 2.171996
#> 
#> $scaling_factor
#> [1] 2.037816
#> 
#> $noise_variance
#> [1] 0.1888017
#> 
#> $nll
#> [1] 35.37379
#> 
fit$RMSE
#> [1] 0.0712868
```
