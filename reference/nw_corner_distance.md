# One-Dimensional Wasserstein Distance Using North-West Corner Algorithm

Applies the north-west corner algorithm to obtain the p-Wasserstein
distance between two one-dimensional empirical distributions.

## Usage

``` r
nw_corner_distance(
  x,
  y,
  presorted = FALSE,
  p = 2,
  eps = sqrt(.Machine$double.eps)
)
```

## Arguments

- x, y:

  Vectors representing one-dimensional empirical distributions.

- presorted:

  Set to `TRUE` if both `x` and `y` are sorted to obtain a speed-up.

- p:

  Order of the Wasserstein distance.

- eps:

  Numerical precision for floating point comparisons.

## Value

The p-Wasserstein distance between `x` and `y`.

## Examples

``` r
x <- rnorm(10)
y <- rnorm(40)
nw_corner_distance(x, y)
#> [1] 0.37397
nw_corner_distance(sort(x), sort(y), presorted = TRUE)
#> [1] 0.37397
```
