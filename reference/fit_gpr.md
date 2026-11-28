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
#> Optimum set at -6.7314.
#> Hyperparameter search 2 of 50.
#> Hyperparameter search 3 of 50.
#> Hyperparameter search 4 of 50.
#> Current optimum improved from -6.7314 to -6.7314.
#> Hyperparameter search 5 of 50.
#> Current optimum improved from -6.7314 to -9.017.
#> Hyperparameter search 6 of 50.
#> Current optimum improved from -9.017 to -9.017.
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
#>  [1]  0.103817306 -0.132066863 -0.003169449  0.001231043  0.005424128
#>  [6]  0.006453146 -0.231320803 -0.120334349 -0.010851893 -0.041995968
#> 
#> $test_variance
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,]  7.161847e-03 -3.079328e-04  1.829867e-03  0.0017314074  0.0017898444
#>  [2,] -3.079328e-04  6.802404e-03  1.938487e-04 -0.0002627019 -0.0018326293
#>  [3,]  1.829867e-03  1.938487e-04  4.647367e-03  0.0016048454  0.0026039162
#>  [4,]  1.731407e-03 -2.627019e-04  1.604845e-03  0.0090115629 -0.0001646103
#>  [5,]  1.789844e-03 -1.832629e-03  2.603916e-03 -0.0001646103  0.0061794055
#>  [6,]  8.578205e-05  9.683021e-04  2.480890e-04  0.0002120663 -0.0000738112
#>  [7,]  5.046756e-07 -1.738177e-03  3.464922e-04  0.0008376046  0.0036702723
#>  [8,] -3.105903e-04  6.667859e-03  1.269006e-04 -0.0003451830 -0.0017545017
#>  [9,] -1.674750e-03 -8.959101e-06 -3.748398e-04  0.0011680375 -0.0012412023
#> [10,] -9.628729e-04  2.091099e-03  8.999166e-05  0.0019651153 -0.0019226019
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,]  8.578205e-05  5.046756e-07 -3.105903e-04 -1.674750e-03 -9.628729e-04
#>  [2,]  9.683021e-04 -1.738177e-03  6.667859e-03 -8.959101e-06  2.091099e-03
#>  [3,]  2.480890e-04  3.464922e-04  1.269006e-04 -3.748398e-04  8.999166e-05
#>  [4,]  2.120663e-04  8.376046e-04 -3.451830e-04  1.168037e-03  1.965115e-03
#>  [5,] -7.381120e-05  3.670272e-03 -1.754502e-03 -1.241202e-03 -1.922602e-03
#>  [6,]  5.682625e-03  1.776491e-03  1.050563e-03 -5.744861e-05  2.144571e-03
#>  [7,]  1.776491e-03  1.045620e-02 -1.445258e-03 -3.631485e-04  1.600935e-03
#>  [8,]  1.050563e-03 -1.445258e-03  6.666453e-03  2.745929e-06  2.170163e-03
#>  [9,] -5.744861e-05 -3.631485e-04  2.745929e-06  6.229939e-03  7.568474e-04
#> [10,]  2.144571e-03  1.600935e-03  2.170163e-03  7.568474e-04  1.021048e-02
#> 
#> $RMSE
#> [1] 0.1799733
#> 
#> $length_scale
#> [1] 1.24624279 0.05817444
#> 
#> $scaling_factor
#> [1] 0.03430252 0.01699882
#> 
#> $noise_variance
#> [1] 0.0100503
#> 
#> $nll
#> [1] -9.01705
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
#>  [1]  0.103817978 -0.132068780 -0.003161914  0.001226479  0.005410689
#>  [6]  0.006468691 -0.231320261 -0.120332927 -0.010843836 -0.041970976
#> 
#> $test_variance
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,]  7.163210e-03 -3.074981e-04  1.830532e-03  0.0017310012  1.789857e-03
#>  [2,] -3.074981e-04  6.804287e-03  1.942298e-04 -0.0002626781 -1.832703e-03
#>  [3,]  1.830532e-03  1.942298e-04  4.648136e-03  0.0016052156  2.603466e-03
#>  [4,]  1.731001e-03 -2.626781e-04  1.605216e-03  0.0090129285 -1.649218e-04
#>  [5,]  1.789857e-03 -1.832703e-03  2.603466e-03 -0.0001649218  6.179819e-03
#>  [6,]  8.597413e-05  9.687326e-04  2.482492e-04  0.0002122410 -7.381729e-05
#>  [7,]  7.028099e-07 -1.739145e-03  3.464667e-04  0.0008386102  3.669855e-03
#>  [8,] -3.102401e-04  6.669668e-03  1.272753e-04 -0.0003451163 -1.754562e-03
#>  [9,] -1.675478e-03 -9.193740e-06 -3.751954e-04  0.0011681298 -1.241399e-03
#> [10,] -9.622629e-04  2.091853e-03  9.069511e-05  0.0019652289 -1.923136e-03
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,]  8.597413e-05  7.028099e-07 -3.102401e-04 -1.675478e-03 -9.622629e-04
#>  [2,]  9.687326e-04 -1.739145e-03  6.669668e-03 -9.193740e-06  2.091853e-03
#>  [3,]  2.482492e-04  3.464667e-04  1.272753e-04 -3.751954e-04  9.069511e-05
#>  [4,]  2.122410e-04  8.386102e-04 -3.451163e-04  1.168130e-03  1.965229e-03
#>  [5,] -7.381729e-05  3.669855e-03 -1.754562e-03 -1.241399e-03 -1.923136e-03
#>  [6,]  5.682890e-03  1.776537e-03  1.050988e-03 -5.745813e-05  2.144961e-03
#>  [7,]  1.776537e-03  1.045840e-02 -1.446143e-03 -3.629595e-04  1.601270e-03
#>  [8,]  1.050988e-03 -1.446143e-03  6.668222e-03  2.534307e-06  2.170910e-03
#>  [9,] -5.745813e-05 -3.629595e-04  2.534307e-06  6.230415e-03  7.568441e-04
#> [10,]  2.144961e-03  1.601270e-03  2.170910e-03  7.568441e-04  1.021186e-02
#> 
#> $RMSE
#> [1] 0.1799746
#> 
#> $length_scale
#> [1] 1.24597202 0.05816952
#> 
#> $scaling_factor
#> [1] 0.03429798 0.01700111
#> 
#> $noise_variance
#> [1] 0.01004939
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
D <- outer(c(x_train, x_test), c(x_train, x_test), "-") |> abs()
fit <- fit_gpr(D, seq_len(N1), N1 + seq_len(N2), y_train, y_test, runs = 50L)
#> Hyperparameter search 1 of 50.
#> Optimum set at 56.7493.
#> Hyperparameter search 2 of 50.
#> Hyperparameter search 3 of 50.
#> Current optimum improved from 56.7493 to 34.3079.
#> Hyperparameter search 4 of 50.
#> Current optimum improved from 34.3079 to 34.3079.
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
#>  [1]  0.12674548  0.04401561  0.20231942 -2.46794517  0.16431600  0.10242345
#>  [7] -0.18537067 -0.24208360  0.06769741  0.16698466
#> 
#> $test_variance
#>                [,1]         [,2]          [,3]          [,4]          [,5]
#>  [1,]  0.0281285682  0.017767786  0.0099658930  1.194179e-04  2.488317e-02
#>  [2,]  0.0177677863  0.027120915 -0.0068130443  1.947970e-04  7.630743e-03
#>  [3,]  0.0099658930 -0.006813044  0.0638297571 -1.127557e-04  2.912011e-02
#>  [4,]  0.0001194179  0.000194797 -0.0001127557  3.185627e-02 -1.191143e-05
#>  [5,]  0.0248831702  0.007630743  0.0291201125 -1.191143e-05  3.021269e-02
#>  [6,]  0.0269533189  0.022659247  0.0016738335  1.896099e-04  1.991931e-02
#>  [7,] -0.0038297557  0.006916789  0.0009365713 -5.542263e-04 -3.806256e-03
#>  [8,]  0.0008284944 -0.003145454  0.0001192397  1.245962e-03  1.461839e-03
#>  [9,]  0.0221047866  0.026598182 -0.0050795312  2.261800e-04  1.227208e-02
#> [10,]  0.0243655612  0.006803844  0.0308315847 -2.091881e-05  3.043375e-02
#>                [,6]          [,7]          [,8]         [,9]         [,10]
#>  [1,]  0.0269533189 -0.0038297557  0.0008284944  0.022104787  2.436556e-02
#>  [2,]  0.0226592470  0.0069167886 -0.0031454543  0.026598182  6.803844e-03
#>  [3,]  0.0016738335  0.0009365713  0.0001192397 -0.005079531  3.083158e-02
#>  [4,]  0.0001896099 -0.0005542263  0.0012459622  0.000226180 -2.091881e-05
#>  [5,]  0.0199193147 -0.0038062556  0.0014618390  0.012272079  3.043375e-02
#>  [6,]  0.0281237947 -0.0019736263 -0.0002016360  0.025912631  1.914147e-02
#>  [7,] -0.0019736263  0.0262216452  0.0038514401  0.002847119 -3.665880e-03
#>  [8,] -0.0002016360  0.0038514401  0.0260416831 -0.002045413  1.451791e-03
#>  [9,]  0.0259126308  0.0028471193 -0.0020454130  0.027632938  1.138555e-02
#> [10,]  0.0191414691 -0.0036658795  0.0014517905  0.011385549  3.072529e-02
#> 
#> $RMSE
#> [1] 0.07643072
#> 
#> $length_scale
#> [1] 1.151216
#> 
#> $scaling_factor
#> [1] 1.463964
#> 
#> $noise_variance
#> [1] 0.1894823
#> 
#> $nll
#> [1] 34.3079
#> 
fit$RMSE
#> [1] 0.07643072
```
