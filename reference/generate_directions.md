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
#>              [,1]        [,2]        [,3]        [,4]
#>  [1,]  0.78950938  0.10494027  0.45009761  0.40382500
#>  [2,]  0.57696299  0.66656186 -0.17245929  0.43939367
#>  [3,]  0.63657217  0.66923123 -0.24363613  0.29588319
#>  [4,] -0.21335308 -0.55370112  0.49662507  0.63345013
#>  [5,] -0.44940206  0.52530858  0.17299691  0.70154169
#>  [6,] -0.48806772  0.03635732 -0.79908094  0.34919578
#>  [7,]  0.39892418 -0.57282875 -0.07182232  0.71243826
#>  [8,]  0.26701676 -0.51094099 -0.81380055  0.07328049
#>  [9,] -0.82875595 -0.09060162  0.09226407  0.54446512
#> [10,]  0.03118836 -0.89215118  0.02843406 -0.44976112
plot(generate_directions(25, 2))
```
