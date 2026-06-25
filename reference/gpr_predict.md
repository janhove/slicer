# Predict Outcomes Using Gaussian Process Regression

Computes posterior mean predictions from a Gaussian process regression
model given precomputed kernel matrices. The implementation follows
Algorithm 2.1 in Rasmussen and Williams (2006), with an optional
centring step and automatic jitter escalation to handle near-singular
kernel matrices.

## Usage

``` r
gpr_predict(
  Kxx,
  Kxstar,
  y_train,
  centre = TRUE,
  lambda2 = 1e-06,
  Kxstarstar = NULL
)
```

## Arguments

- Kxx:

  Kernel matrix evaluated at the training inputs.

- Kxstar:

  Matrix of kernel values between test and training inputs. Each row
  corresponds to one test input.

- y_train:

  Vector with training outcomes.

- centre:

  If `TRUE`, the training outcomes are centred around their mean.

- lambda2:

  Kernel noise variance.

- Kxstarstar:

  Optional. A square matrix of kernel values among the test inputs. If
  `NULL` (default), only posterior means are computed and returned as a
  vector. Else both the posterior mean vector and covariance matrix are
  returned.

## Value

If `Kxstarstar` is `NULL` (default), a vector with posterior mean
predictions for the test inputs. Else, a list with a vector `mean` with
posterior mean predictions and a matrix `var` with the posterior
covariance for the test inputs.

## Examples

``` r
N1 <- 25
N2 <- 10
x_train <- seq(-pi, pi, length.out = N1)
x_test  <- runif(N2, -pi, pi)
distance <- outer(c(x_train, x_test), c(x_train, x_test), "-") |> abs()
kernel <- rbf(distance, length_scale = 1, variance = 1)
Kxx    <- kernel[1:N1, 1:N1]
Kxstar <- kernel[(N1+1):(N1+N2), 1:N1]
Kxstarstar <- kernel[(N1+1):(N1+N2), (N1+1):(N1+N2)]
y_train <- ifelse(x_train == 0, 2*pi, sin(2*pi*x_train) / x_train)
curve(sin(2*pi*x)/x, from = -pi, to = pi)
points(x_train, y_train)
# Only means: don't use Kxstarstar
points(x_test, gpr_predict(Kxx, Kxstar, y_train), pch = 16)
# Add credible intervals around predictions:
gpr_fit <- gpr_predict(Kxx, Kxstar, y_train, Kxstarstar = Kxstarstar)
segments(x0 = x_test,
         y0 = gpr_fit$mean - 2 * sqrt(diag(gpr_fit$var)),
         y1 = gpr_fit$mean + 2 * sqrt(diag(gpr_fit$var)), col = "#4DAF4A")
```
