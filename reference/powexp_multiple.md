# Sum of Power Exponential Kernels

Builds power exponential kernel matrices from a list of distances
matrices and computes their sum.

## Usage

``` r
powexp_multiple(D_list, length_scales, variances, expons, squared = TRUE)
```

## Arguments

- D_list:

  A list of matrices with pairwise (possibly squared) distances.

- length_scales:

  A vector of length-scale parameters.

- variances:

  A vector of kernel variance parameters.

- expons:

  A vector of exponents. If only a single value is provided, it will be
  used for all components.

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
M_list <- list(M1, M2)
powexp_multiple(M_list, length_scales = c(1, 2), variances = c(2, 1/2),
  expons = c(0.5, 1.5), squared = FALSE)
#>           [,1]      [,2]      [,3]
#> [1,] 2.5000000 0.8373277 0.7653118
#> [2,] 0.8373277 2.5000000 0.4911193
#> [3,] 0.7653118 0.4911193 2.5000000
```
