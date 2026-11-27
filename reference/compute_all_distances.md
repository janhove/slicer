# Pairwise Wasserstein Distances Along Projection Directions Between Several Empirical Distributions

This function computes Wasserstein distances between the projections of
several empirical directions along specified projection directions. The
projection directions can optionally be transformed by means of a linear
map. The functions can output the distances along each projection
direction, or the distances averaged over projection directions. If the
projection directions are uniformly sampled from the unit sphere (using
[`generate_directions()`](https://janhove.github.io/slicer/reference/generate_directions.md)),
the averaged version estimates the sliced Wasserstein distances.

## Usage

``` r
compute_all_distances(
  distributions,
  thetas = NULL,
  A = NULL,
  p = 2,
  verbose = TRUE,
  average = FALSE,
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

- p:

  Order of the Wasserstein distances. Defaults to 2.

- verbose:

  If `TRUE`, show progress.

- average:

  If `FALSE`, the distance matrix for each projection direction is
  output. If `TRUE`, the distance matrices for the different projection
  directions are averaged. See details below.

- test_idx:

  Optionally, a vector of indices. If supplied, skip the distance
  computations between distribution pairs whose indices both occur in
  this vector.

## Value

If `keep_projections = TRUE`, a list of matrices containing
p-Wasserstein distances raised to the p-th power, one for each
projection direction; otherwise, a matrix with the averages of the
p-Wasserstein distances raised to the p-th power.

## Averaging

When using `average = TRUE`, the distances are raised to the p-th power,
averaged for each distribution pair, and then the p-th root is computed.

## Examples

``` r
M1 <- matrix(rnorm(50), ncol = 5)
M2 <- matrix(rnorm(50), ncol = 5)
M3 <- matrix(rnorm(250), ncol = 5)

# Sliced Wasserstein of order 2
my_directions <- generate_directions(20, 5)
compute_all_distances(list(M1, M2, M3), my_directions,
  average = TRUE, verbose = FALSE)
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.6104348 0.4788046
#> [2,] 0.6104348 0.0000000 0.5193320
#> [3,] 0.4788046 0.5193320 0.0000000

# Sliced Wasserstein of order 1
compute_all_distances(list(M1, M2, M3), my_directions, p = 1,
  average = TRUE, verbose = FALSE)
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.4833153 0.3750074
#> [2,] 0.4833153 0.0000000 0.4101071
#> [3,] 0.3750074 0.4101071 0.0000000

# Marginal squared 2-Wasserstein distances:
marginal_wass <- compute_all_distances(list(M1, M2, M3), diag(1, 5),
  average = FALSE, verbose = FALSE)
marginal_wass[[3]] # along third dimension
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.6347678 0.3158111
#> [2,] 0.6347678 0.0000000 0.5919710
#> [3,] 0.3158111 0.5919710 0.0000000

# Reweight projection directions (shear transformation)
A <- diag(c(4, 0.5, 3, 2, 1))
shear_wass <- compute_all_distances(list(M1, M2, M3), diag(1, 5), A = A,
  average = FALSE, verbose = FALSE)
shear_wass[[3]]
#>           [,1]     [,2]      [,3]
#> [1,] 0.0000000 1.904303 0.9474332
#> [2,] 1.9043034 0.000000 1.7759130
#> [3,] 0.9474332 1.775913 0.0000000
```
