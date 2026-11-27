# Gaussian Radial Basis Function Kernel

Builds a Gaussian radial basis function kernel matrix from a distance
matrix.

## Usage

``` r
rbf(D, length_scale = 1, variance = 1, squared = FALSE)
```

## Arguments

- D:

  Matrix with pairwise (possibly squared) distances.

- length_scale:

  Length-scale parameter.

- variance:

  Kernel variance parameter.

- squared:

  If `TRUE`, `D` is assumed to contain squared distances. The default is
  `FALSE`.

## Value

A kernel matrix with the same dimensions as `D`.

## Examples

``` r
M <- rbind(c(0, 2, 1),
           c(2, 0, 2.5),
           c(1, 2.5, 0))
rbf(M)
#>           [,1]       [,2]       [,3]
#> [1,] 1.0000000 0.13533528 0.60653066
#> [2,] 0.1353353 1.00000000 0.04393693
#> [3,] 0.6065307 0.04393693 1.00000000
```
