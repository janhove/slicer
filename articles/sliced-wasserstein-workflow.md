# (Sliced) Wasserstein - Gaussian process regression workflow

## Goal

This vignette illustrates the main use of the `slicer` package by means
of a toy example. Empirical distributions are generated from bivariate
Gaussian mixtures with four isotropic components laid out symmetrically
around the origin. These distributions differ from one another in that
they are rotated counterclockwise by an angle $`\omega \in [0, \pi/2)`$.
The goal of the analysis is to retrieve the specific angles by which the
distributions have been rotated.

## Functions

Let’s load `slicer`.

``` r

library(slicer)
```

To generate the input distributions, we define `generate_data()`:

``` r

generate_data <- function(n = 200, a = 5, angle = 0, Sigma = diag(1, 2)) {
  mus <- cbind(c(a, 0), c(0, a), c(-a, 0), c(0, -a))
  R <- rbind(
    c(cos(angle), -sin(angle)),
    c(sin(angle), cos(angle))
  )
  cluster <- sample(1:4, size = n, replace = TRUE)
  d <- matrix(0, nrow = n, ncol = 2)
  for (i in 1:4) {
    d[which(cluster == i), ] <- MASS::mvrnorm(
      sum(cluster == i), mu = mus[, i], Sigma = Sigma
    )
  }
  d %*% t(R)
}
```

The plots below show two empirical distributions, once using a rotation
angle of $`\omega = \pi/16`$ and once using one of $`\omega = 3\pi/8`$.

``` r

set.seed(2026) # for reproducibility

op <- par(no.readonly = TRUE, bty = "n")
par(mfrow = c(1, 2))
generate_data(angle = pi/16) |> 
  plot(xlab = "first dimension", ylab = "second dimension", asp = 1)
generate_data(angle = 3*pi/8) |> 
  plot(xlab = "first dimension", ylab = "second dimension", asp = 1)
```

![](sliced-wasserstein-workflow_files/figure-html/unnamed-chunk-3-1.png)

``` r

par(op)
```

## Data generation

We generate 60 matrices (40 for training, 20 for testing) containing
empirical distributions with the rotation angles sampled uniformly from
$`[0, \pi/4)`$. Each empirical distribution is based on 200 data points.

``` r

sample_size <- 200
N_train <- 40
N_test  <- 20
angle_range <- c(0, pi/4)

angles <- runif(N_train + N_test, angle_range[1], angle_range[2])

distributions <- vector("list", N_train + N_test)
for (i in seq_len(N_train + N_test)) {
  distributions[[i]] <- generate_data(angle = angles[i])
}
```

## Distances

We can use
[`compute_all_distances()`](https://janhove.github.io/slicer/reference/compute_all_distances.md)
to compute pairwise distances between the input distributions. If we
supply the optional vector `test_idx` with the test set indices, the
pairwise distance computations between test objects (which are not
needed for obtaining predictions) are skipped. To compute **sliced
Wasserstein distances**, we need to first generate
$`\boldsymbol{\theta_1}, \dots, \boldsymbol{\theta_L}`$ sampled
uniformly at random from the unit sphere;
[`generate_directions()`](https://janhove.github.io/slicer/reference/generate_directions.md)
takes care of this. We’ll use $`L = 25`$ and set $`d = 2`$ since we’re
working in two dimensions. If we’re only interested in the estimated
sliced Wasserstein distances, we can set `keep_projections` to `FALSE`.

``` r

thetas <- generate_directions(L = 25, d = 2)
sw_distances <- compute_all_distances(distributions, thetas, verbose = FALSE,
    keep_projections = FALSE, test_idx = N_train + seq_len(N_test))
```

The $`60 \times 60`$ matrix `sw_distances` contains the *squared*
(estimated) sliced Wasserstein distances among the forty training
objects and between the training objects and the twenty test objects.
The pairwise distances among the test objects are set to `NA`:

``` r

sw_distances # output not shown
```

We can also compute (normal) **Wasserstein distances** between the
distributions when they are projected along certain dimensions. For
instances, to compute the pairwise Wasserstein distances along the first
margin and along the second margin, we set the projection directions to
$`\boldsymbol{e_1} = (1, 0)^{\top}`$ and
$`\boldsymbol{e_2} = (0, 1)^{\top}`$. We also set
`keep_projections = TRUE`, which will cause the function to output a
list of two matrices with *squared* Wasserstein distances: one for each
margin.

``` r

marginal_distances <- compute_all_distances(distributions, diag(1, 2), 
    verbose = FALSE, keep_projections = TRUE, 
    test_idx = N_train + seq_len(N_test))
str(marginal_distances)
#> List of 2
#>  $ : num [1:60, 1:60] 0 1.274 0.931 0.548 0.649 ...
#>  $ : num [1:60, 1:60] 0 2.177 0.121 0.445 0.599 ...
```

## Gaussian process models with tuned hyperparameters

The function
[`fit_gpr()`](https://janhove.github.io/slicer/reference/fit_gpr.md)
takes a single matrix with squared pairwise distances and uses it as
input to a Gaussian process regression model with (by default) a
Gaussian RBF kernel:
``` math
k_{\textrm{RBF}}(x, x') = s^2\exp\left(-\frac{d^2(x, x')}{2t^2}\right)
```
, where $`s^2 > 0`$ is a scaling factor, $`t > 0`$ is the length-scale,
and $`d(x, x')`$ is the distance between the inputs $`x`$ and $`x'`$.

The model’s and the kernel’s hyperparameters are tuned using the
training data by minimising the negative marginal log-likelihood.

``` r

sw_fit <- fit_gpr(sw_distances,
  training_idx = seq_len(N_train),
  test_idx = N_train + seq_len(N_test),
  y_train = angles[seq_len(N_train)], 
  y_test = angles[N_train + seq_len(N_test)],
  verbose = TRUE)
#> Hyperparameter search 1 of 10.
#> Optimum set at -58.815.
#> Hyperparameter search 2 of 10.
#> Current optimum improved from -58.815 to -58.8157.
#> Hyperparameter search 3 of 10.
#> Hyperparameter search 4 of 10.
#> Current optimum improved from -58.8157 to -58.8172.
#> Hyperparameter search 5 of 10.
#> Current optimum improved from -58.8172 to -58.8176.
#> Hyperparameter search 6 of 10.
#> Hyperparameter search 7 of 10.
#> Hyperparameter search 8 of 10.
#> Hyperparameter search 9 of 10.
#> Hyperparameter search 10 of 10.
```

The output includes predictions for the test data, the covariannce for
the predictions (if the true test outcomes were provided), and the
estimated hyperparameters.

``` r

str(sw_fit)
#> List of 9
#>  $ kernels         : chr "rbf"
#>  $ params          : logi NA
#>  $ test_predictions: num [1:20] 0.241 0.7091 0.2402 0.5074 0.0604 ...
#>  $ test_variance   : num [1:20, 1:20] 0.00134 NA NA NA NA ...
#>  $ RMSE            : num 0.038
#>  $ length_scale    : num 1.57
#>  $ scaling_factor  : num 0.0758
#>  $ noise_variance  : num 2.06e-07
#>  $ nll             : num -58.8
plot(angles[N_train + seq_len(N_test)], sw_fit$test_predictions,
     xlab = "true test outcomes", ylab = "predicted test outcomes", asp = 1)
```

![](sliced-wasserstein-workflow_files/figure-html/unnamed-chunk-9-1.png)

Since we didn’t compute the distances between the test inputs, only the
predictions variances are reported (the main diagonal of
`test_variance`):

``` r

diag(sw_fit$test_variance)
#>  [1] 0.0013366046 0.0021747134 0.0009410526 0.0024842650 0.0009932470
#>  [6] 0.0012133368 0.0026981347 0.0017291192 0.0009257147 0.0010941824
#> [11] 0.0014113696 0.0014526947 0.0008645598 0.0009060875 0.0030615161
#> [16] 0.0008919558 0.0017045107 0.0039858230 0.0025413739 0.0023768660
```

Had we also computed the distances between the test inputs,
`test_variance` would have also included the covariances.

The function
[`fit_gpr()`](https://janhove.github.io/slicer/reference/fit_gpr.md) can
also be used when multiple matrices with squared distances are provided.
Now, estimated length-scale and (kernel) variance hyperparameters are
provided for the Gaussian RBF corresponding to each squared distance
matrix.

``` r

marginal_fit <- fit_gpr(marginal_distances,
  training_idx = seq_len(N_train),
  test_idx = N_train + seq_len(N_test),
  y_train = angles[seq_len(N_train)], 
  y_test = angles[N_train + seq_len(N_test)],
  verbose = TRUE)
#> Hyperparameter search 1 of 10.
#> Optimum set at -61.5134.
#> Hyperparameter search 2 of 10.
#> Hyperparameter search 3 of 10.
#> Hyperparameter search 4 of 10.
#> Hyperparameter search 5 of 10.
#> Hyperparameter search 6 of 10.
#> Hyperparameter search 7 of 10.
#> Hyperparameter search 8 of 10.
#> Hyperparameter search 9 of 10.
#> Hyperparameter search 10 of 10.
str(marginal_fit)
#> List of 9
#>  $ kernels         : chr [1:2] "rbf" "rbf"
#>  $ params          : logi [1:2] NA NA
#>  $ test_predictions: num [1:20] 0.2394 0.6678 0.2428 0.5619 0.0918 ...
#>  $ test_variance   : num [1:20, 1:20] 0.000587 NA NA NA NA ...
#>  $ RMSE            : num 0.0522
#>  $ length_scale    : num [1:2] 3.89 6.04
#>  $ scaling_factor  : num [1:2] 0.345 0.144
#>  $ noise_variance  : num 0.000448
#>  $ nll             : num -61.5
plot(angles[N_train + seq_len(N_test)], marginal_fit$test_predictions,
     xlab = "true test outcomes", ylab = "predicted test outcomes", asp = 1)
```

![](sliced-wasserstein-workflow_files/figure-html/unnamed-chunk-11-1.png)

We can combine the marginal and sliced Wasserstein distances into a list
with three distance matrices, too:

``` r

total_fit <- fit_gpr(list(sw_distances, marginal_distances[[1]], marginal_distances[[2]]),
  training_idx = seq_len(N_train),
  test_idx = N_train + seq_len(N_test),
  y_train = angles[seq_len(N_train)], 
  y_test = angles[N_train + seq_len(N_test)],
  verbose = TRUE)
#> Hyperparameter search 1 of 10.
#> Optimum set at -63.0742.
#> Hyperparameter search 2 of 10.
#> Hyperparameter search 3 of 10.
#> Hyperparameter search 4 of 10.
#> Hyperparameter search 5 of 10.
#> Current optimum improved from -63.0742 to -63.0791.
#> Hyperparameter search 6 of 10.
#> Hyperparameter search 7 of 10.
#> Hyperparameter search 8 of 10.
#> Hyperparameter search 9 of 10.
#> Hyperparameter search 10 of 10.
str(total_fit)
#> List of 9
#>  $ kernels         : chr [1:3] "rbf" "rbf" "rbf"
#>  $ params          : logi [1:3] NA NA NA
#>  $ test_predictions: num [1:20] 0.2474 0.6891 0.2471 0.5495 0.0852 ...
#>  $ test_variance   : num [1:20, 1:20] 0.000804 NA NA NA NA ...
#>  $ RMSE            : num 0.0471
#>  $ length_scale    : num [1:3] 0.419 3.834 5.509
#>  $ scaling_factor  : num [1:3] 0.00114 0.24911 0.126
#>  $ noise_variance  : num 1.68e-07
#>  $ nll             : num -63.1
plot(angles[N_train + seq_len(N_test)], total_fit$test_predictions,
     xlab = "true test outcomes", ylab = "predicted test outcomes", asp = 1)
```

![](sliced-wasserstein-workflow_files/figure-html/unnamed-chunk-12-1.png)

Parallel processing can be enabled using the `cores` parameter:

``` r

total_fit <- fit_gpr(list(sw_distances, marginal_distances[[1]], marginal_distances[[2]]),
  training_idx = seq_len(N_train),
  test_idx = N_train + seq_len(N_test),
  y_train = angles[seq_len(N_train)],
  y_test = angles[N_train + seq_len(N_test)],
  runs = 50L, cores = 2L)
str(total_fit)
#> List of 9
#>  $ kernels         : chr [1:3] "rbf" "rbf" "rbf"
#>  $ params          : logi [1:3] NA NA NA
#>  $ test_predictions: num [1:20] 0.2474 0.6896 0.2471 0.5491 0.0849 ...
#>  $ test_variance   : num [1:20, 1:20] 0.00081 NA NA NA NA ...
#>  $ RMSE            : num 0.047
#>  $ length_scale    : num [1:3] 0.421 3.843 5.515
#>  $ scaling_factor  : num [1:3] 0.00117 0.24792 0.12613
#>  $ noise_variance  : num 3.12e-07
#>  $ nll             : num -63.1
```

Instead of the Gaussian RBF kernel, the Matérn kernels with smoothness
parameter $`\nu`$ can be used. We can use a different kernel for each
distance matrix; the total kernel used is then the conical combination
of these kernels. For instance, the following fit uses the Matérn kernel
with $`\nu = 2.5`$ for the sliced Wasserstein distances, the Matérn
kernel with $`\nu = 0.5`$ for the Wasserstein distances along the first
axis, and the Matérn kernel with $`\nu = 1.5`$ for the Wasserstein
distances alogn the second axis.

``` r

total_fit_matern <- fit_gpr(list(sw_distances, marginal_distances[[1]], marginal_distances[[2]]),
  training_idx = seq_len(N_train),
  test_idx = N_train + seq_len(N_test),
  y_train = angles[seq_len(N_train)],
  y_test = angles[N_train + seq_len(N_test)],
  kernels = c("matern", "matern", "matern"),
  params = c(2.5, 0.5, 1.5),
  runs = 50L, cores = 2L)
str(total_fit_matern)
#> List of 9
#>  $ kernels         : chr [1:3] "matern" "matern" "matern"
#>  $ params          : num [1:3] 2.5 0.5 1.5
#>  $ test_predictions: num [1:20] 0.242 0.709 0.24 0.507 0.061 ...
#>  $ test_variance   : num [1:20, 1:20] 0.00133 NA NA NA NA ...
#>  $ RMSE            : num 0.0377
#>  $ length_scale    : num [1:3] 2.71 1.99e+03 6.38e+05
#>  $ scaling_factor  : num [1:3] 1.34e-01 1.55e-06 3.89e-06
#>  $ noise_variance  : num 5.6e-09
#>  $ nll             : num -58.6
```

Power exponential (`"powexp"`) kernels are available, too. They take an
exponent $`\kappa \in (0, 2]`$ as their parameter; here,
$`\kappa_1 = 1.2, \kappa_2 = 1, \kappa_3 = 0.8`$:

``` r

total_fit_powexp <- fit_gpr(list(sw_distances, marginal_distances[[1]], marginal_distances[[2]]),
  training_idx = seq_len(N_train),
  test_idx = N_train + seq_len(N_test),
  y_train = angles[seq_len(N_train)],
  y_test = angles[N_train + seq_len(N_test)],
  kernels = c("powexp", "powexp", "powexp"),
  params = c(1.2, 1, 0.8),
  runs = 50L, cores = 2L)
str(total_fit_powexp)
#> List of 9
#>  $ kernels         : chr [1:3] "powexp" "powexp" "powexp"
#>  $ params          : num [1:3] 1.2 1 0.8
#>  $ test_predictions: num [1:20] 0.2678 0.6818 0.2302 0.5324 0.0807 ...
#>  $ test_variance   : num [1:20, 1:20] 0.00276 NA NA NA NA ...
#>  $ RMSE            : num 0.055
#>  $ length_scale    : num [1:3] 11.6 21.3 430.7
#>  $ scaling_factor  : num [1:3] 6.70e-02 1.23e-01 2.47e-06
#>  $ noise_variance  : num 7.72e-09
#>  $ nll             : num -47.9
```

Different kernels can be applied to different distances matrices:

``` r

total_fit_mixed <- fit_gpr(list(sw_distances, marginal_distances[[1]], marginal_distances[[2]]),
  training_idx = seq_len(N_train),
  test_idx = N_train + seq_len(N_test),
  y_train = angles[seq_len(N_train)],
  y_test = angles[N_train + seq_len(N_test)],
  kernels = c("matern", "powexp", "rbf"),
  params = c(2.8, 0.9, NA),
  runs = 50L, cores = 2L)
str(total_fit_mixed)
#> List of 9
#>  $ kernels         : chr [1:3] "matern" "powexp" "rbf"
#>  $ params          : num [1:3] 2.8 0.9 NA
#>  $ test_predictions: num [1:20] 0.235 0.7046 0.2401 0.5036 0.0611 ...
#>  $ test_variance   : num [1:20, 1:20] 0.00135 NA NA NA NA ...
#>  $ RMSE            : num 0.0374
#>  $ length_scale    : num [1:3] 2.41 6057.61 3.41
#>  $ scaling_factor  : num [1:3] 1.06e-01 2.16e-07 2.27e-02
#>  $ noise_variance  : num 1.03e-07
#>  $ nll             : num -58.7
```

Finally, a GPR’s **negative log predictive density** (NLPD) for some
test data can be obtained like so:

``` r

nlpd_gpr(total_fit_mixed, angles[N_train + seq_len(N_test)])
#> Warning in nlpd_gpr(total_fit_mixed, angles[N_train + seq_len(N_test)]): The fit's test variance contains NAs, likely because the distances among the test inputs weren't provided.
#> Replacing NAs by 0.
#> [1] -1.947789
#> attr(,"jitter_used")
#> [1] 0
```

The warning message occurs because we didn’t compute the distances among
the test inputs themselves. As a result, the kernel matrices don’t
specify the covariance among the test inputs, and the NLPD assumes that
they are zero. If you remove the `test_idx` parameter setting when
running
[`compute_all_distances()`](https://janhove.github.io/slicer/reference/compute_all_distances.md)
and refit the models, this warning will disappear.
