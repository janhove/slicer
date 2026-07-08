# Matérn Kernels

Builds a Matérn kernel (with smoothness parameter 0.5, 1.5 or 2.5) from
a distance matrix.

## Usage

``` r
matern(D, length_scale = 1, variance = 1, tau = 1.5, squared = TRUE)
```

## Arguments

- D:

  Matrix with pairwise (possibly squared) distances.

- length_scale:

  Length-scale parameter.

- variance:

  Kernel variance parameter.

- tau:

  Smoothness parameter. needs to be 0.5, 1.5 (default) or 2.5.

- squared:

  If `TRUE` (default), `D` is assumed to contain squared distances. Set
  to `FALSE` if this is not the case.

## Value

A kernel matrix with the same dimensions as `D`.

## Examples

``` r
M <- rbind(c(0, 2, 1),
           c(2, 0, 2.5),
           c(1, 2.5, 0))
matern(M, tau = 0.5, squared = FALSE)
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.3678794 0.6065307
#> [2,] 0.3678794 1.0000000 0.2865048
#> [3,] 0.6065307 0.2865048 1.0000000
matern(M, tau = 1.5, squared = FALSE)
#>           [,1]       [,2]       [,3]
#> [1,] 1.0000000 0.13973135 0.48335772
#> [2,] 0.1397314 1.00000000 0.07017579
#> [3,] 0.4833577 0.07017579 1.00000000
matern(M, tau = 2.5, squared = FALSE)
#>           [,1]       [,2]       [,3]
#> [1,] 1.0000000 0.13866022 0.52399411
#> [2,] 0.1386602 1.00000000 0.06351021
#> [3,] 0.5239941 0.06351021 1.00000000
```
