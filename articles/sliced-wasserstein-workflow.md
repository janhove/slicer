# (Sliced) Wasserstein kernels and Gaussian process regression

## Goal

This vignette illustrates the main use of the `slicer` package by means
of a toy example. Empirical distributions are generated from bivariate
Gaussian mixtures with four isotropic components laid out symmetrically
around the origin. These distributions differ from one another in that
they are rotated counterclockwise by an angle $`\omega \in [0, \pi/2)`$.
The goal of the analysis is to retrieve the specific angles by which the
distributions were rotated.

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

set.seed(20261127)

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
