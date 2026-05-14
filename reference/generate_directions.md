# Generate random projection directions

Obtain `L` random vectors on the unit sphere in `d`-dimensional
Euclidean space.

## Usage

``` r
generate_directions(L, d)
```

## Arguments

- L:

  Number of projection directions.

- d:

  Dimension.

## Value

An `L` by `d` matrix, each row of which contains a unit vector.

## Examples

``` r
generate_directions(10, 4)
#>              [,1]       [,2]        [,3]         [,4]
#>  [1,] -0.11383381  0.4323414 -0.07970842  0.890937348
#>  [2,] -0.02487093  0.5449624  0.14185975 -0.825998352
#>  [3,] -0.90399859  0.1152910  0.28803495 -0.294160509
#>  [4,]  0.35569561 -0.9212104 -0.15670579 -0.017182480
#>  [5,]  0.04206901 -0.4218964  0.20647188 -0.881818015
#>  [6,]  0.06063430 -0.5024280 -0.17696780  0.844139774
#>  [7,] -0.45031604 -0.3337077 -0.39787643 -0.726325686
#>  [8,]  0.46584574 -0.1046329 -0.84791751  0.230381415
#>  [9,] -0.33685618  0.7113108  0.21448007  0.578414301
#> [10,]  0.31178657 -0.4823502  0.81855934 -0.009383175
plot(generate_directions(25, 2))
```
