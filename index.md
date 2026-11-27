# slicer: Sliced Wasserstein Distances for Regression

`slicer` provides functions for computing sliced Wasserstein distances
as well as one-dimensional Wasserstein distances between empirical
distributions represented as matrices. These distances can be converted
into a kernel (Gaussian radial basis function, Matérn and power
exponential) and then fed to a Gaussian process regression model.

## Installation

You can install the development version of `slicer` from
[GitHub](https://github.com/) with:

``` r

# install.packages("pak")
pak::pak("janhove/slicer")
```

## Use

Please refer to the vignettes.

## To do

- Unit testing.
- Vignette explaining (sliced) Wasserstein kernels and their use in
  Gaussian process regression.
