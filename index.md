# slicer: Sliced Wasserstein Distances for Regression

`slicer` provides functions for computing sliced Wasserstein distances
as well as one-dimensional Wasserstein distances between empirical
distributions represented as matrices. These distances can be fed to a
Gaussian process regression model, which currently accepts the Gaussian
radial basis function kernel as well as the Matérn kernels with
smoothness parameters $`\nu = 0.5, 1.5, 2.5`$.

## Installation

You can install the development version of slicer from
[GitHub](https://github.com/) with:

``` r

# install.packages("pak")
pak::pak("janhove/slicer")
```

## Use

Please refer to the vignette for a brief tutorial.
