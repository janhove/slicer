# Sum of Matérn Kernels

Builds Matérn kernel matrices from a list of distances matrices and
computes their sum.

## Usage

``` r
matern_multiple(D_list, length_scales, variances, taus, squared = TRUE)
```

## Arguments

- D_list:

  A list of matrices with pairwise (possibly squared) distances.

- length_scales:

  A vector of length-scale parameters.

- variances:

  A vector of kernel variance parameters.

- taus:

  A vector of smoothness parameters. Values need to be in 0.5, 1.5 and
  2.5. If only a single value is provided, it will be used for all
  components.

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
matern_multiple(M_list, length_scales = c(1, 2), variances = c(2, 1/2),
  taus = c(0.5, 2.5), squared = FALSE)
#>          [,1]      [,2]      [,3]
#> [1,] 2.500000 1.1500835 1.2823914
#> [2,] 1.150083 2.5000000 0.7145912
#> [3,] 1.282391 0.7145912 2.5000000
```
