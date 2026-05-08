# slicer: Sliced Wasserstein Distances for Regression

`slicer` provides functions for computing sliced Wasserstein distances
as well as one-dimensional Wasserstein distances between empirical
distributions represented as matrices. The distances can be fed to a
kernel function (currently, only the Gaussian radial basis function
kernel is supported), which can then be used in a Gaussian process
regression model.

## Installation

You can install the development version of slicer from
[GitHub](https://github.com/) with:

``` r

# install.packages("pak")
pak::pak("janhove/slicer")
```

## Documentation

See <https://janhove.github.io/slicer> for documentation.

## Use

Refer to the vignette for a brief tutorial.
