# Sum of Kernels

Builds kernel matrices from a list of distances matrices and computes
their sum.

## Usage

``` r
kernel_sum(D2_list, length_scales, variances, kernels)
```

## Arguments

- D2_list:

  A list of matrices with squared pairwise distances.

- length_scales:

  A vector of length-scale parameters.

- variances:

  A vector of kernel variance parameters.

- kernels:

  A vector of kernels (`rbf`, `matern05`, `matern15`, `matern25`).

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
M_list <- list(M1^2, M2^2)
kernel_sum(M_list, length_scales = c(1, 2), variances = c(2, 1/2),
  kernels = c("rbf", "matern15"))
#>           [,1]      [,2]      [,3]
#> [1,] 2.5000000 0.6631144 1.2829270
#> [2,] 0.6631144 2.5000000 0.2217522
#> [3,] 1.2829270 0.2217522 2.5000000
```
