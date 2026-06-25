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
#>              [,1]       [,2]        [,3]        [,4]
#>  [1,]  0.86655316  0.1667153 -0.37197610  0.28796775
#>  [2,]  0.92526042 -0.2520103 -0.14534940 -0.24342863
#>  [3,] -0.53683919  0.6050859  0.08095148 -0.58234147
#>  [4,]  0.06805144 -0.3282733  0.57597267  0.74556095
#>  [5,] -0.24419160 -0.5198069  0.78482336  0.23285952
#>  [6,] -0.25996781 -0.4977248 -0.72171839 -0.40473366
#>  [7,] -0.29779626 -0.8404266  0.17249429  0.41862426
#>  [8,]  0.58119687 -0.6419864 -0.30473829  0.39648225
#>  [9,] -0.21775591 -0.1270875 -0.61284329  0.74890202
#> [10,] -0.82992581  0.5525779 -0.02754795 -0.07156759
plot(generate_directions(25, 2))
```
