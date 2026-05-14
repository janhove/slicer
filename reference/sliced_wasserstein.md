# Sliced Wasserstein Distance Between Two Empirical Distributions

Estimate the sliced Wasserstein distance between two empirical
distributions represented as matrices using Monte Carlo simulation.

Optionally, a matrix `thetas` can be supplied. Each row of `thetas` will
be interpreted as a projection direction. These rows do not have to be
unit vectors, allowing users to transform the projection directions
rather than the matrices.

## Usage

``` r
sliced_wasserstein(x, y, p = 2, thetas = NULL, L = 50, seed = NULL)
```

## Arguments

- x, y:

  Matrices representing empirical distributions.

- p:

  Order of sliced Wasserstein distance.

- thetas:

  Optionally, a matrix, each row of which represents a projection
  direction.

- L:

  If no `thetas` are provided, `L` random projection directions are
  generated.

- seed:

  Optional random seed.

## Value

Estimated sliced Wasserstein distance of order `p` between `x` and `y`.

## Examples

``` r
M1 <- matrix(rnorm(50), ncol = 5)
M2 <- matrix(rnorm(250), ncol = 5)
sliced_wasserstein(M1, M2) # random projection directions
#> [1] 0.6139966
cardinal_axes <- diag(1, 5)
sliced_wasserstein(M1, M2, thetas = cardinal_axes)
#> [1] 0.6674493
first_two_axes <- cardinal_axes[1:2, ]
sliced_wasserstein(M1, M2, thetas = first_two_axes)
#> [1] 0.6105983
```
