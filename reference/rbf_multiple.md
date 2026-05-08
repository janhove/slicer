# Sum of Gaussian radial basis function kernels

Builds Gaussian radial basis function kernel matrices from a list of
distances matrices and computes their sum.

## Usage

``` r
rbf_multiple(D_list, length_scales, variances, squared = TRUE)
```

## Arguments

- D_list:

  A list of matrices with pairwise (possibly squared) distances.

- length_scales:

  A vector of length-scale parameters.

- variances:

  A vector of kernel variance parameters.

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
rbf_multiple(M_list, length_scales = c(1, 2), variances = c(2, 1/2), squared = FALSE)
#>          [,1]      [,2]      [,3]
#> [1,] 2.500000 0.7119190 1.2807290
#> [2,] 0.711919 2.5000000 0.2502001
#> [3,] 1.280729 0.2502001 2.5000000
```
