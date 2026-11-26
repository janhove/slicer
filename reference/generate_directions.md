# Generate Random Projection Directions

Obtain `L` vectors sampled uniformly at random from the unit sphere in
`d`-dimensional Euclidean space.

## Usage

``` r
generate_directions(L, d)
```

## Arguments

- L:

  Number of projection directions.

- d:

  Dimensionality of the ambient space.

## Value

An `L` by `d` matrix, each row of which contains a unit vector.

## Examples

``` r
generate_directions(10, 4)
#>             [,1]        [,2]        [,3]        [,4]
#>  [1,]  0.5390429  0.09525974 -0.78934144 -0.27802593
#>  [2,]  0.7060359 -0.46395968 -0.34381833 -0.40993126
#>  [3,]  0.1941825  0.46961090 -0.10547860 -0.85477075
#>  [4,] -0.6296354 -0.32081232  0.67743238  0.20426477
#>  [5,] -0.0241143  0.31170963 -0.48223114  0.81835734
#>  [6,]  0.6310542  0.19053111 -0.19146696  0.72719247
#>  [7,]  0.4724374  0.11056033 -0.03986437  0.87349306
#>  [8,]  0.4324194 -0.76390225 -0.47517810  0.06060161
#>  [9,]  0.5038601  0.14898804  0.30647199 -0.79372697
#> [10,] -0.1192164 -0.52583209  0.08356404 -0.83803647
plot(generate_directions(25, 2))
```
