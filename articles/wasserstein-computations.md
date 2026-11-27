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
$`p \in [1, \infty)`$. This is defined as

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
that the input distributions are all univariate. By default, it outputs
the pairwise $`p`$-Wasserstein distances *raised to the $`p`$-th_power*.
(The reason is that, in practice, we will often ) Also by default,
$`p = 2`$ is used.

``` r

P1 <- runif(14, -1, 1)
P2 <- rcauchy(39)
P3 <- rnorm(41)
P4 <- 2*(rbinom(47, 1, 0.4) - 1.5)
pw_dists <- compute_all_distances(list(P1, P2, P3, P4))
#> No projection directions (thetas) supplied. Assuming univariate distributions.
pw_dists
#>            [,1]     [,2]       [,3]      [,4]
#> [1,]  0.0000000 13.96197  0.3926126  4.010920
#> [2,] 13.9619674  0.00000 11.3016628 18.016477
#> [3,]  0.3926126 11.30166  0.0000000  4.670274
#> [4,]  4.0109200 18.01648  4.6702741  0.000000
```

The actual $`\mathbf{W}_2`$ distances can now be obtain by taking the
square root:

``` r

sqrt(pw_dists)
#>           [,1]     [,2]      [,3]     [,4]
#> [1,] 0.0000000 3.736572 0.6265881 2.002728
#> [2,] 3.7365716 0.000000 3.3617946 4.244582
#> [3,] 0.6265881 3.361795 0.0000000 2.161082
#> [4,] 2.0027281 4.244582 2.1610817 0.000000
```

Different orders ($`p`$) can be specified, too:

``` r

compute_all_distances(list(P1, P2, P3, P4), p = 1)
#> No projection directions (thetas) supplied. Assuming univariate distributions.
#>           [,1]     [,2]      [,3]     [,4]
#> [1,] 0.0000000 1.412198 0.4530223 1.896643
#> [2,] 1.4121982 0.000000 1.2034373 2.483512
#> [3,] 0.4530223 1.203437 0.0000000 2.075587
#> [4,] 1.8966435 2.483512 2.0755873 0.000000
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
is.list(distributions)
#> [1] TRUE
length(distributions)
#> [1] 6
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
projection directions), each of which contains squared 2-Wasserstein
distances:

``` r

# Along fourth canonical basis vector
all_dists[[4]]
#>           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
#> [1,]  0.000000 14.404932  7.581264 31.425098 46.213927  3.357988
#> [2,] 14.404932  0.000000  3.461265  4.851487  9.977685  4.144342
#> [3,]  7.581264  3.461265  0.000000  9.053734 17.962454  1.483038
#> [4,] 31.425098  4.851487  9.053734  0.000000  1.942366 14.662450
#> [5,] 46.213927  9.977685 17.962454  1.942366  0.000000 24.906594
#> [6,]  3.357988  4.144342  1.483038 14.662450 24.906594  0.000000

# Along theta_9
all_dists[[9]]
#>           [,1]       [,2]       [,3]       [,4]      [,5]       [,6]
#> [1,]  0.000000  3.0972480  6.9837456  5.0764543 30.675007  4.5659659
#> [2,]  3.097248  0.0000000 16.5037894 13.4191350 49.443632  0.2261069
#> [3,]  6.983746 16.5037894  0.0000000  0.2240779  9.147574 19.4126790
#> [4,]  5.076454 13.4191350  0.2240779  0.0000000 11.465228 16.1293612
#> [5,] 30.675007 49.4436324  9.1475739 11.4652281  0.000000 54.4684200
#> [6,]  4.565966  0.2261069 19.4126790 16.1293612 54.468420  0.0000000
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
but set the `keep_projections` parameter to `FALSE`. The function now
works just as before, i.e., it computes $`\mathbf{W}_p^p`$ values
between the projected distributions, but it now averages the $`L`$
values for each distribution pair. By default, squared sliced
Wasserstein distances of order $`p = 2`$ are output:

``` r

sw_dists <- compute_all_distances(distributions, Theta, 
                                  keep_projections = FALSE)
#> Projecting...
#> Computing distances...
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   2%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   6%  |                                                                              |======                                                                |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  12%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  16%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  22%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  28%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  32%  |                                                                              |========================                                              |  34%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  38%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  42%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  46%  |                                                                              |==================================                                    |  48%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  58%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  62%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  66%  |                                                                              |================================================                      |  68%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  72%  |                                                                              |====================================================                  |  74%  |                                                                              |=====================================================                 |  76%  |                                                                              |=======================================================               |  78%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  82%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  88%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  92%  |                                                                              |==================================================================    |  94%  |                                                                              |===================================================================   |  96%  |                                                                              |===================================================================== |  98%  |                                                                              |======================================================================| 100%
sw_dists
#>           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
#> [1,]  0.000000 10.025981 12.103010  9.053233 12.079050 11.420344
#> [2,] 10.025981  0.000000  8.685864  6.618335  8.985870  8.797133
#> [3,] 12.103010  8.685864  0.000000  5.966392  6.935964 12.782313
#> [4,]  9.053233  6.618335  5.966392  0.000000  6.134122 10.016689
#> [5,] 12.079050  8.985870  6.935964  6.134122  0.000000 18.033213
#> [6,] 11.420344  8.797133 12.782313 10.016689 18.033213  0.000000
```

Of course, sliced Wasserstein distances of other orders can be computed.
This outputs the (raw) sliced Wasserstein distances of order $`p = 1`$:

``` r

sw_dists_1 <- compute_all_distances(distributions, Theta, p = 1,
                                    keep_projections = FALSE)
#> Projecting...
#> Computing distances...
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   2%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   6%  |                                                                              |======                                                                |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  12%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  16%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  22%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  28%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  32%  |                                                                              |========================                                              |  34%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  38%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  42%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  46%  |                                                                              |==================================                                    |  48%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  58%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  62%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  66%  |                                                                              |================================================                      |  68%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  72%  |                                                                              |====================================================                  |  74%  |                                                                              |=====================================================                 |  76%  |                                                                              |=======================================================               |  78%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  82%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  88%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  92%  |                                                                              |==================================================================    |  94%  |                                                                              |===================================================================   |  96%  |                                                                              |===================================================================== |  98%  |                                                                              |======================================================================| 100%
sw_dists_1
#>          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
#> [1,] 0.000000 2.538208 2.944338 2.546741 2.749255 2.747485
#> [2,] 2.538208 0.000000 2.400416 2.048687 2.595333 2.480271
#> [3,] 2.944338 2.400416 0.000000 1.980011 2.227120 2.914414
#> [4,] 2.546741 2.048687 1.980011 0.000000 2.015000 2.596229
#> [5,] 2.749255 2.595333 2.227120 2.015000 0.000000 3.553161
#> [6,] 2.747485 2.480271 2.914414 2.596229 3.553161 0.000000
```
