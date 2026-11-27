# How to compute Wasserstein and sliced Wasserstein distances

## Goal

This vignette explains how you can use `slicer` to compute Wasserstein
distances between univariate empirical distributions as well as sliced
Wasserstein distances between multivariate empirical distributions.

## Wasserstein distances between univariate empirical distributions

Let’s first load the `slicer` package and set a random seed.

``` r

library(slicer)
set.seed(20261127)
```

Given two empirical distributions on $`\mathbb{R}`$,

``` math
P = \frac{1}{m}\sum_{i=1}^m \delta_{x_i}
```
and

``` math
Q = \frac{1}{n}\sum_{j=1}^n \delta_{y_j},
```

we want to compute the Wasserstein distance of order
$`p \in [1, \infty)`$ between $`P`$ and $`Q`$. This is defined as

``` math
\mathbf{W}_p(P, Q) := \min_{\substack{X \sim P \\ Y \sim Q}} \left(\mathbb{E}|X - Y|^p\right)^{1/p}.
```

The distributions $`P`$ and $`Q`$ are represented as the vectors (or
matrices) $`x = (x_1, \dots, x_m) \in \mathbb{R}^{m \times 1}`$ and
$`y = (y_1, \dots, y_n) \in \mathbb{R}^{n \times 1}`$, respectively. For
instance,

``` r

P <- c(rnorm(20, 14, 3), rnorm(8, 1, 8))
Q <- rcauchy(14, location = 9, scale = 4)
```

The **north-west corner** algorithm is used to compute exact Wasserstein
distances of order $`p`$. The function
[`nw_corner_distance()`](https://janhove.github.io/slicer/reference/nw_corner_distance.md)
runs it.

``` r

nw_corner_distance(P, Q, p = 2)
#> [1] 12.35423
nw_corner_distance(P, Q, p = 1)
#> [1] 7.530863
```

This algorithm is essentially a sorting operation with some bells and
whistles. When working with large presorted vectors, you can set the
parameter value for `presorted` to `TRUE` for some speed-up.

Typically, we will have several empirical distributions
$`P_1, \dots, P_N`$. The function
[`compute_all_distances()`](https://janhove.github.io/slicer/reference/compute_all_distances.md)
takes a list of empirical distributions represented as vectors (or
matrices). If we don’t specify anything to the contrary, it will assume
that the input distributions are all univariate. By default, Wasserstein
distances of order $`p = 2`$ are computed.

``` r

P1 <- P
P2 <- Q
P3 <- rnorm(41)
P4 <- 2*(rbinom(47, 1, 0.4) - 1.5)
pw_dists <- compute_all_distances(list(P1, P2, P3, P4))
#> No projection directions (thetas) supplied. Assuming univariate distributions.
pw_dists
#>          [,1]     [,2]      [,3]      [,4]
#> [1,]  0.00000 12.35423 12.963304 15.023061
#> [2,] 12.35423  0.00000 19.953468 21.676767
#> [3,] 12.96330 19.95347  0.000000  2.252589
#> [4,] 15.02306 21.67677  2.252589  0.000000
```

Different orders ($`p`$) can be specified, too:

``` r

compute_all_distances(list(P1, P2, P3, P4), p = 1)
#> No projection directions (thetas) supplied. Assuming univariate distributions.
#>           [,1]      [,2]      [,3]      [,4]
#> [1,]  0.000000  7.530863 12.059093 14.045166
#> [2,]  7.530863  0.000000 14.134749 16.180809
#> [3,] 12.059093 14.134749  0.000000  2.094052
#> [4,] 14.045166 16.180809  2.094052  0.000000
```

## Wasserstein distances between projected multivariate distributions

While it is possible to define Wasserstein distances between
$`d`$-dimensional distributions for $`d \geq 2`$, these distances are
expensive to compute. Moreover, they are not Hilbertian. However, what
we can do is define a number of vectors
$`\theta_1, \dots, \theta_L \in \mathbb{R}^d`$ and compute
$`p`$-Wasserstein distances between the projections of the distributions
onto these vectors.

More precisely, let $`P_1, \dots, P_N`$ be $`d`$-dimensional empirical
distributions and let $`\theta_1, \dots, \theta_L`$ be projection
vectors in $`\mathbb{R}^d`$. We now compute

``` math
\mathbf{W}_p({\theta_{\ell}}_{\#}^{\star} P_i, {\theta_{\ell}}_{\#}^{\star} P_j)
```
for all $`\ell = 1, \dots, L, 1 \leq i < j \leq N`$, where
$`v_{\#}^{\star}Q`$ is the pushforward (image) measure of $`Q`$ under
the map $`x \mapsto v^{\top}x`$.

Using `slicer`, we can compute these $`L \cdot N(N-1)/2`$ pairwise
distances as follows. We represent each $`d`$-dimensional empirical
distribution

``` math
P_i = \frac{1}{n_i}\sum_{k = 1}^{n_i} \delta_{x_{i, k}},
```
with $`x_{i,k} \in \mathbb{R}^d`$, as a matrix,

``` math
M_i = 
\begin{pmatrix}
  x^{\top}_{i, 1} \\
  \vdots \\
  x^{\top}_{i, n_i}
\end{pmatrix} \in \mathbb{R}^{n_i \times d}.
```

We put these matrices $`M_1, \dots, M_N`$ in a list.

Then, we define a matrix

``` math
\Theta = \begin{pmatrix} \theta^{\top}_1 \\ \vdots \\ \theta^{\top}_L \end{pmatrix} \in \mathbb{R}^{L \times d}
```
that contains the projection vectors. If we’re just interested in
obtaining Wasserstein distances between the corresponding marginal
distributions, we set $`\Theta = I_d`$.

Then, we feed the list of $`M`$-matrices and $`\Theta`$ to
[`compute_all_distances()`](https://janhove.github.io/slicer/reference/compute_all_distances.md).

By way of example, we generate a number ($`N`$) of samples from
different Gaussian mixtures with a randomly chosen number of components
$`c ~ \textrm{Unif}(1, 2, \dots, C)`$ using the following code.

``` r

rgmm_random <- function(n, C, d) {
  p <- sample.int(C, 1) 
  mu <- matrix(runif(p * d, -5, 5), nrow = p, ncol = d)
  chol_list <- lapply(seq_len(p), function(j) {
    a <- runif(1, 1, 4)
    A <- matrix(runif(d * d, -1, 1), d, d)
    B <- diag(runif(d, 0, 1), nrow = d)
    Sigma <- a * tcrossprod(A) + B
    chol(Sigma)       
  })
  z <- sample.int(p, n, replace = TRUE)
  X <- matrix(0, nrow = n, ncol = d)
  for (j in seq_len(p)) {
    idx <- which(z == j)
    m <- length(idx)
    if (m > 0) {
      Z <- matrix(rnorm(m * d), m, d)
      X[idx, ] <- Z %*% chol_list[[j]] +
        matrix(mu[j, ], m, d, byrow = TRUE)
    }
  }
  list(distribution = X, components = p)
}
rgmm_random_N <- function(N = 100, n = 50, C = 4, d = 8) {
  replicate(N, rgmm_random(n, C, d), simplify = FALSE)
}
```

We generate $`N = 6`$ 8-dimensional empirical distributions by drawing
samples of size $`n = 500`$ from Gaussian mixture models with between 1
and $`C = 4`$ mixture components. (It is not necessary for the
distributions to be supported on the same number of points.)

``` r

rgmms <- rgmm_random_N(N = 6, n = 500, C = 4, d = 8)
distributions <- lapply(rgmms, `[[`, "distribution") 
```

If we want to compute the pairwise Wasserstein distances between the
corresponding margins (i.e., after projecting the distributions onto
each of the cardinal axes) as well as those between the projections
onto, say, the vector

``` math
\theta_9 = \frac{1}{\sqrt{3}}(0, 1, 0, 1, 0, 0, 0, 1)^{\top},
```
then we put the canonical basis vectors $`e_1, \dots, e_8`$ as well as
$`\theta_9`$ into a matrix:

``` r

Theta <- rbind(
  diag(1, 8),
  1/sqrt(3) * c(0, 1, 0, 1, 0, 0, 0, 1)
)
Theta
#>       [,1]      [,2] [,3]      [,4] [,5] [,6] [,7]      [,8]
#>  [1,]    1 0.0000000    0 0.0000000    0    0    0 0.0000000
#>  [2,]    0 1.0000000    0 0.0000000    0    0    0 0.0000000
#>  [3,]    0 0.0000000    1 0.0000000    0    0    0 0.0000000
#>  [4,]    0 0.0000000    0 1.0000000    0    0    0 0.0000000
#>  [5,]    0 0.0000000    0 0.0000000    1    0    0 0.0000000
#>  [6,]    0 0.0000000    0 0.0000000    0    1    0 0.0000000
#>  [7,]    0 0.0000000    0 0.0000000    0    0    1 0.0000000
#>  [8,]    0 0.0000000    0 0.0000000    0    0    0 1.0000000
#>  [9,]    0 0.5773503    0 0.5773503    0    0    0 0.5773503
```

The squared 2-Wasserstein distances between the projections onto each
row of `Theta` can now be computed as follows:

``` r

all_dists <- compute_all_distances(distributions, Theta)
#> Projecting...
#> Computing distances...
#>   |                                                                              |                                                                      |   0%  |                                                                              |========                                                              |  11%  |                                                                              |================                                                      |  22%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================                                       |  44%  |                                                                              |=======================================                               |  56%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================                |  78%  |                                                                              |==============================================================        |  89%  |                                                                              |======================================================================| 100%
```

The object `all_dists` is a list with 9 elements (the number of
projection directions), each of which contains Wasserstein distances of
order $`p = 2`$:

``` r

# Along fourth canonical basis vector
all_dists[[4]]
#>          [,1]     [,2]      [,3]     [,4]      [,5]     [,6]
#> [1,] 0.000000 2.511518 2.3869818 3.511217 2.5519653 1.079536
#> [2,] 2.511518 0.000000 4.1255069 5.956054 4.6548244 3.089585
#> [3,] 2.386982 4.125507 0.0000000 2.837708 0.9174726 1.379150
#> [4,] 3.511217 5.956054 2.8377080 0.000000 2.0281511 3.005094
#> [5,] 2.551965 4.654824 0.9174726 2.028151 0.0000000 1.615499
#> [6,] 1.079536 3.089585 1.3791504 3.005094 1.6154990 0.000000

# Along theta_9
all_dists[[9]]
#>          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
#> [1,] 0.000000 2.135175 1.975174 4.061653 2.578371 1.765189
#> [2,] 2.135175 0.000000 2.176853 5.892668 3.707649 2.853606
#> [3,] 1.975174 2.176853 0.000000 4.348005 1.969027 1.420678
#> [4,] 4.061653 5.892668 4.348005 0.000000 2.737301 3.201242
#> [5,] 2.578371 3.707649 1.969027 2.737301 0.000000 1.043625
#> [6,] 1.765189 2.853606 1.420678 3.201242 1.043625 0.000000
```

## Sliced Wasserstein distances

To estimate sliced Wasserstein distances between $`d`$-variate empirical
probability distributions, we generate independent
$`\theta_1, \dots, \theta_L`$ uniformly at random from the unit sphere
in $`\mathbb{R}^d`$. The function
[`generate_directions()`](https://janhove.github.io/slicer/reference/generate_directions.md)
takes care of this:

``` r

Theta <- generate_directions(L = 50, d = 8)
```

We again use
[`compute_all_distances()`](https://janhove.github.io/slicer/reference/compute_all_distances.md)
but set the `average` parameter to `TRUE`. The function now averages the
$`L`$ Wasserstein distances computed for each pair of distributions like
so:

``` math
\widehat{\mathbf{SW}}_p(P_i, P_j) = \left(\frac{1}{L}\sum_{\ell = 1}^L \mathbf{W}^p_p({\theta_{\ell}}_{\#}^{\star} P_i, {\theta_{\ell}}_{\#}^{\star} P_j)\right)^{1/p}.
```
By default, squared sliced Wasserstein distances of order $`p = 2`$ are
output:

``` r

sw_dists <- compute_all_distances(distributions, Theta, average = TRUE)
#> Projecting...
#> Computing distances...
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   2%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   6%  |                                                                              |======                                                                |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  12%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  16%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  22%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  28%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  32%  |                                                                              |========================                                              |  34%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  38%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  42%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  46%  |                                                                              |==================================                                    |  48%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  58%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  62%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  66%  |                                                                              |================================================                      |  68%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  72%  |                                                                              |====================================================                  |  74%  |                                                                              |=====================================================                 |  76%  |                                                                              |=======================================================               |  78%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  82%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  88%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  92%  |                                                                              |==================================================================    |  94%  |                                                                              |===================================================================   |  96%  |                                                                              |===================================================================== |  98%  |                                                                              |======================================================================| 100%
sw_dists
#>          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
#> [1,] 0.000000 2.225677 2.601511 2.679255 2.085345 3.394862
#> [2,] 2.225677 0.000000 3.108849 3.218715 2.279165 2.964327
#> [3,] 2.601511 3.108849 0.000000 3.201459 2.690495 4.084576
#> [4,] 2.679255 3.218715 3.201459 0.000000 2.424931 4.234922
#> [5,] 2.085345 2.279165 2.690495 2.424931 0.000000 2.595893
#> [6,] 3.394862 2.964327 4.084576 4.234922 2.595893 0.000000
```

Of course, sliced Wasserstein distances of other orders can be computed.
This outputs the sliced Wasserstein distances of order $`p = 1`$:

``` r

sw_dists_1 <- compute_all_distances(distributions, Theta, p = 1, average = TRUE)
#> Projecting...
#> Computing distances...
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   2%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   6%  |                                                                              |======                                                                |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  12%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  16%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  22%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  28%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  32%  |                                                                              |========================                                              |  34%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  38%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  42%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  46%  |                                                                              |==================================                                    |  48%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  58%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  62%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  66%  |                                                                              |================================================                      |  68%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  72%  |                                                                              |====================================================                  |  74%  |                                                                              |=====================================================                 |  76%  |                                                                              |=======================================================               |  78%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  82%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  88%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  92%  |                                                                              |==================================================================    |  94%  |                                                                              |===================================================================   |  96%  |                                                                              |===================================================================== |  98%  |                                                                              |======================================================================| 100%
sw_dists_1
#>          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
#> [1,] 0.000000 1.845777 2.059437 2.154981 1.692770 2.736052
#> [2,] 1.845777 0.000000 2.523656 2.570453 1.860875 2.516199
#> [3,] 2.059437 2.523656 0.000000 2.645451 2.270762 3.513010
#> [4,] 2.154981 2.570453 2.645451 0.000000 1.881794 3.427599
#> [5,] 1.692770 1.860875 2.270762 1.881794 0.000000 2.068663
#> [6,] 2.736052 2.516199 3.513010 3.427599 2.068663 0.000000
```
