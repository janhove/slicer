# Power Exponential Kernel

Builds a power exponential kernel from a distance matrix.

## Usage

``` r
powexp(D, length_scale = 1, variance = 1, expon = 1, squared = FALSE)
```

## Arguments

- D:

  Matrix with pairwise (possibly squared) distances.

- length_scale:

  Length-scale parameter.

- variance:

  Kernel variance parameter.

- expon:

  Exponent. Defaults to 1.

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
powexp(M, expon = 0.5)
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.2431167 0.3678794
#> [2,] 0.2431167 1.0000000 0.2057407
#> [3,] 0.3678794 0.2057407 1.0000000
powexp(M, expon = 1)
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.1353353 0.3678794
#> [2,] 0.1353353 1.0000000 0.0820850
#> [3,] 0.3678794 0.0820850 1.0000000
powexp(M, expon = 2)
#>            [,1]        [,2]        [,3]
#> [1,] 1.00000000 0.018315639 0.367879441
#> [2,] 0.01831564 1.000000000 0.001930454
#> [3,] 0.36787944 0.001930454 1.000000000
```
