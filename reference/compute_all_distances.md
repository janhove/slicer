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
  keep_projections = TRUE,
  test_idx = NULL
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

- test_idx:

  Optionally, a vector of indices. If supplied, skip the distance
  computations between distribution pairs whose indices both occur in
  this vector.

## Value

If `keep_projections = TRUE`, a list of squared-distance matrices, one
for each projection direction; otherwise, a matrix with the averaged
squared distances.

## Examples

``` r
M1 <- matrix(rnorm(50), ncol = 5)
M2 <- matrix(rnorm(50), ncol = 5)
M3 <- matrix(rnorm(250), ncol = 5)
# Sliced Wasserstein:
my_directions <- generate_directions(20, 5)
compute_all_distances(list(M1, M2, M3), my_directions,
  keep_projections = FALSE, verbose = FALSE)
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.3726307 0.2292539
#> [2,] 0.3726307 0.0000000 0.2697058
#> [3,] 0.2292539 0.2697058 0.0000000
# Marginal Wasserstein distances:
marginal_wass <- compute_all_distances(list(M1, M2, M3), diag(1, 5),
  keep_projections = TRUE, verbose = FALSE)
marginal_wass[[3]] # along third dimension
#>            [,1]      [,2]       [,3]
#> [1,] 0.00000000 0.4029302 0.09973663
#> [2,] 0.40293015 0.0000000 0.35042965
#> [3,] 0.09973663 0.3504296 0.00000000
# Reweight projection directions
A <- diag(c(4, 0.5, 3, 2, 1))
shear_wass <- compute_all_distances(list(M1, M2, M3), diag(1, 5), A = A,
  keep_projections = TRUE, verbose = FALSE)
shear_wass[[3]]
#>           [,1]     [,2]      [,3]
#> [1,] 0.0000000 3.626371 0.8976297
#> [2,] 3.6263714 0.000000 3.1538668
#> [3,] 0.8976297 3.153867 0.0000000
```
