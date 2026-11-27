# Matérn Kernel

Builds a Matérn kernel from a distance matrix.

## Usage

``` r
matern(D, length_scale = 1, variance = 1, nu = 1.5, squared = FALSE)
```

## Arguments

- D:

  Matrix with pairwise (possibly squared) distances.

- length_scale:

  Length-scale parameter.

- variance:

  Kernel variance parameter.

- nu:

  Smoothness parameter. Defaults to 1.5.

- squared:

  If `TRUE`, `D` is assumed to contain squared distances. Defaults to
  `FALSE`.

## Value

A kernel matrix with the same dimensions as `D`.

## Examples

``` r
M <- rbind(c(0, 2, 1),
           c(2, 0, 2.5),
           c(1, 2.5, 0))
matern(M, nu = 0.5)
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.1353353 0.3678794
#> [2,] 0.1353353 1.0000000 0.0820850
#> [3,] 0.3678794 0.0820850 1.0000000
matern(M, nu = 1.5)
#>           [,1]       [,2]       [,3]
#> [1,] 1.0000000 0.13973135 0.48335772
#> [2,] 0.1397314 1.00000000 0.07017579
#> [3,] 0.4833577 0.07017579 1.00000000
matern(M, nu = 2.5)
#>           [,1]       [,2]       [,3]
#> [1,] 1.0000000 0.13866022 0.52399411
#> [2,] 0.1386602 1.00000000 0.06351021
#> [3,] 0.5239941 0.06351021 1.00000000
```
