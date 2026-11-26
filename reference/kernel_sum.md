# Sum of Kernels

Builds kernel matrices from a list of distances matrices and computes
their sum.

## Usage

``` r
kernel_sum(D_list, length_scales, variances, kernels, params, squared = TRUE)
```

## Arguments

- D_list:

  A list of matrices with pairwise (possibly squared) distances.

- length_scales:

  A vector of length-scale parameters.

- variances:

  A vector of kernel variance parameters.

- kernels:

  A vector of kernels (`rbf`, `matern`, `powexp`).

- params:

  A vector of additional kernel parameters (smoothness parameter for
  Matérn, exponent for power exponential; ignored for RBF).

- squared:

  If `TRUE` (default), each distance matrix is assumed to contain
  squared distances. If `FALSE`, each distance matrix is assumed to
  contain raw distances.

## Value

A kernel matrix with the same dimensions as the matrices in `D_list`.

## Examples

``` r
M1 <- rbind(c(0, 2, 1),
            c(2, 0, 2.5),
            c(1, 2.5, 0))
M2 <- rbind(c(0, 1, 4),
            c(1, 0, 3),
            c(4, 3, 0))
M3 <- rbind(c(0, 4, 8),
            c(4, 0, 1),
            c(8, 1, 0))
M_list <- list(M1^2, M2^2, M3^2)
K_sum <- kernel_sum(M_list, length_scales = c(1, 2, 0.5),
                    variances = c(2, 1/2, 0.1),
                    kernels = c("powexp", "matern", "rbf"), params = c(1, 2.3, NA))
```
