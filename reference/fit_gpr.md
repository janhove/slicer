# Gaussian Process Fitting to Test Data with Tuned Hyperparameters

This function generates predictions using a Gaussian process model with
one or several Gaussian RBF kernels. The hyperparameters are tuned using
the training data by minimising the model's negative marginal
log-likelihood.

## Usage

``` r
fit_gpr(
  D,
  training_idx,
  test_idx,
  y_train,
  y_test = NULL,
  centre = TRUE,
  use_gradient = TRUE,
  kernels = "rbf",
  params = NA,
  squared = FALSE,
  runs = 10L,
  cores = 1L,
  verbose = TRUE,
  maxit = 500
)
```

## Arguments

- D:

  A matrix, or a list of matrices, containing pairwise distances between
  all objects.

- training_idx:

  A vector with the row (and column) numbers corresponding to the
  training entries in `D`.

- test_idx:

  A vector with the row (and column) numbers corresponding to the test
  entries in `D`.

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

- squared:

  Set to `TRUE` if the distances in `D` are squared. Defaults to
  `FALSE`.

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
D1 <- outer(c(x_train1, x_test1), c(x_train1, x_test1), "-") |> abs()
D2 <- outer(c(x_train2, x_test2), c(x_train2, x_test2), "-") |> abs()

# Single core
fit <- fit_gpr(list(D1, D2), seq_len(N1), N1 + seq_len(N2),
  y_train, y_test, runs = 50)
#> Hyperparameter search 1 of 50.
#> Optimum set at -9.017.
#> Hyperparameter search 2 of 50.
#> Hyperparameter search 3 of 50.
#> Hyperparameter search 4 of 50.
#> Hyperparameter search 5 of 50.
#> Current optimum improved from -9.017 to -9.017.
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
#>  [1]  0.103997531 -0.132056414 -0.003188060  0.001400180  0.005431839
#>  [6]  0.006513740 -0.231348566 -0.120349995 -0.010917488 -0.041974402
#> 
#> $test_variance
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,]  7.165613e-03 -0.0003035898  1.834287e-03  0.0017402327  1.789816e-03
#>  [2,] -3.035898e-04  0.0068044566  1.928180e-04 -0.0002674102 -1.831816e-03
#>  [3,]  1.834287e-03  0.0001928180  4.650179e-03  0.0016076576  2.604381e-03
#>  [4,]  1.740233e-03 -0.0002674102  1.607658e-03  0.0090087138 -1.679565e-04
#>  [5,]  1.789816e-03 -0.0018318158  2.604381e-03 -0.0001679565  6.180408e-03
#>  [6,]  8.674228e-05  0.0009698991  2.495675e-04  0.0002130453 -7.435612e-05
#>  [7,] -2.872418e-07 -0.0017431459  3.553835e-04  0.0008373512  3.675912e-03
#>  [8,] -3.050739e-04  0.0066698214  1.258847e-04 -0.0003503113 -1.754120e-03
#>  [9,] -1.677148e-03 -0.0000134011 -3.755255e-04  0.0011672339 -1.242415e-03
#> [10,] -9.632311e-04  0.0020916799  9.222944e-05  0.0019781523 -1.924121e-03
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,]  8.674228e-05 -2.872418e-07 -3.050739e-04 -1.677148e-03 -9.632311e-04
#>  [2,]  9.698991e-04 -1.743146e-03  6.669821e-03 -1.340110e-05  2.091680e-03
#>  [3,]  2.495675e-04  3.553835e-04  1.258847e-04 -3.755255e-04  9.222944e-05
#>  [4,]  2.130453e-04  8.373512e-04 -3.503113e-04  1.167234e-03  1.978152e-03
#>  [5,] -7.435612e-05  3.675912e-03 -1.754120e-03 -1.242415e-03 -1.924121e-03
#>  [6,]  5.688540e-03  1.777998e-03  1.052262e-03 -5.916460e-05  2.147847e-03
#>  [7,]  1.777998e-03  1.045600e-02 -1.449787e-03 -3.646279e-04  1.605202e-03
#>  [8,]  1.052262e-03 -1.449787e-03  6.667826e-03 -1.523998e-06  2.171395e-03
#>  [9,] -5.916460e-05 -3.646279e-04 -1.523998e-06  6.233058e-03  7.565077e-04
#> [10,]  2.147847e-03  1.605202e-03  2.171395e-03  7.565077e-04  1.022379e-02
#> 
#> $RMSE
#> [1] 0.1799876
#> 
#> $length_scale
#> [1] 1.24562886 0.05832583
#> 
#> $scaling_factor
#> [1] 0.03437295 0.01704841
#> 
#> $noise_variance
#> [1] 0.01005711
#> 
#> $nll
#> [1] -9.017025
#> 

# Multiple cores
fit <- fit_gpr(list(D1, D2), seq_len(N1), N1 + seq_len(N2),
 y_train, y_test, runs = 50, cores = 2L)
fit
#> $kernels
#> [1] "rbf" "rbf"
#> 
#> $params
#> [1] NA NA
#> 
#> $test_predictions
#>  [1]  0.103858546 -0.132057111 -0.003164795  0.001262679  0.005393361
#>  [6]  0.006481531 -0.231332599 -0.120326573 -0.010855295 -0.041947791
#> 
#> $test_variance
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,]  7.163606e-03 -3.061337e-04  1.831678e-03  0.0017326507  1.789451e-03
#>  [2,] -3.061337e-04  6.804520e-03  1.940299e-04 -0.0002638420 -1.832041e-03
#>  [3,]  1.831678e-03  1.940299e-04  4.648470e-03  0.0016057624  2.602895e-03
#>  [4,]  1.732651e-03 -2.638420e-04  1.605762e-03  0.0090111597 -1.658757e-04
#>  [5,]  1.789451e-03 -1.832041e-03  2.602895e-03 -0.0001658757  6.179089e-03
#>  [6,]  8.631170e-05  9.692573e-04  2.486055e-04  0.0002124678 -7.394697e-05
#>  [7,]  6.976350e-07 -1.740366e-03  3.486936e-04  0.0008388970  3.670203e-03
#>  [8,] -3.086168e-04  6.669879e-03  1.271007e-04 -0.0003463358 -1.754020e-03
#>  [9,] -1.676075e-03 -1.043906e-05 -3.754616e-04  0.0011677708 -1.241500e-03
#> [10,] -9.616845e-04  2.091981e-03  9.162169e-05  0.0019680927 -1.923317e-03
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,]  8.631170e-05  6.976350e-07 -3.086168e-04 -1.676075e-03 -9.616845e-04
#>  [2,]  9.692573e-04 -1.740366e-03  6.669879e-03 -1.043906e-05  2.091981e-03
#>  [3,]  2.486055e-04  3.486936e-04  1.271007e-04 -3.754616e-04  9.162169e-05
#>  [4,]  2.124678e-04  8.388970e-04 -3.463358e-04  1.167771e-03  1.968093e-03
#>  [5,] -7.394697e-05  3.670203e-03 -1.754020e-03 -1.241500e-03 -1.923317e-03
#>  [6,]  5.683407e-03  1.776451e-03  1.051514e-03 -5.790466e-05  2.145520e-03
#>  [7,]  1.776451e-03  1.045733e-02 -1.447281e-03 -3.631661e-04  1.602263e-03
#>  [8,]  1.051514e-03 -1.447281e-03  6.668282e-03  1.335957e-06  2.171184e-03
#>  [9,] -5.790466e-05 -3.631661e-04  1.335957e-06  6.230246e-03  7.564981e-04
#> [10,]  2.145520e-03  1.602263e-03  2.171184e-03  7.564981e-04  1.021388e-02
#> 
#> $RMSE
#> [1] 0.1799808
#> 
#> $length_scale
#> [1] 1.24554139 0.05820321
#> 
#> $scaling_factor
#> [1] 0.03428545 0.01701015
#> 
#> $noise_variance
#> [1] 0.01004912
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
D <- outer(c(x_train, x_test), c(x_train, x_test), "-") |> abs()
fit <- fit_gpr(D, seq_len(N1), N1 + seq_len(N2), y_train, y_test, runs = 50L)
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
