# Pairwise Wasserstein Distances Along Projection Directions Between Several Empirical Distributions

This function computes the squared 2-Wasserstein distances between the
projections of several empirical directions along specified projection
directions. The projection directions can optionally be transformed by
means of a linear map. The functions can output the squared distances
along each projection direction, or the squared distances averaged over
projection directions.

## Usage

``` r
compute_all_distances(
  distributions,
  thetas,
  A = NULL,
  verbose = TRUE,
  keep_projections = TRUE
)
```

## Arguments

- distributions:

  A list of matrices representing empirical distributions.

- thetas:

  A matrix, each row of which represents a projection direction.

- A:

  Optionally, a matrix used to transform each projection direction.

- verbose:

  If `TRUE`, show progress.

- keep_projections:

  If `TRUE`, the distance matrix for each projection direction is
  output. If `FALSE`, the distance matrices for the different projection
  directions are averaged.

## Value

If `keep_projections = TRUE`, a list of squared-distance matrices, one
for each projection direction; otherwise, a matrix with the averaged
squared distances.

## Examples

``` r
M1 <- matrix(rnorm(50), ncol = 5)
M2 <- matrix(rnorm(150), ncol = 5)
M3 <- matrix(rnorm(250), ncol = 5)
# Sliced Wasserstein:
my_directions <- generate_directions(20, 5)
compute_all_distances(list(M1, M2, M3), my_directions,
  keep_projections = FALSE, verbose = FALSE)
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.2662411 0.2123149
#> [2,] 0.2662411 0.0000000 0.2030276
#> [3,] 0.2123149 0.2030276 0.0000000
# Marginal Wasserstein distances:
marginal_wass <- compute_all_distances(list(M1, M2, M3), diag(1, 5),
  keep_projections = TRUE, verbose = FALSE)
marginal_wass[[3]] # along third dimension
#>           [,1]       [,2]       [,3]
#> [1,] 0.0000000 0.11677560 0.15943208
#> [2,] 0.1167756 0.00000000 0.05865423
#> [3,] 0.1594321 0.05865423 0.00000000
# Reweight projection directions
A <- diag(c(4, 0.5, 3, 2, 1))
shear_wass <- compute_all_distances(list(M1, M2, M3), diag(1, 5), A = A,
  keep_projections = TRUE, verbose = FALSE)
shear_wass[[3]]
#>          [,1]      [,2]      [,3]
#> [1,] 0.000000 1.0509804 1.4348888
#> [2,] 1.050980 0.0000000 0.5278881
#> [3,] 1.434889 0.5278881 0.0000000
```
