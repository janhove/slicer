#' One-Dimensional Wasserstein Distance Using North-West Corner Algorithm
#'
#' Applies the north-west corner algorithm to obtain the p-Wasserstein distance
#' between two one-dimensional empirical distributions.
#'
#' @param x,y Vectors representing one-dimensional empirical distributions.
#' @param presorted Set to `TRUE` if both `x` and `y` are sorted to obtain a speed-up.
#' @param p Order of the Wasserstein distance.
#' @param eps Numerical precision for floating point comparisons.
#'
#' @return The p-Wasserstein distance between `x` and `y`.
#' @export
#'
#' @examples
#' x <- rnorm(10)
#' y <- rnorm(40)
#' nw_corner_distance(x, y)
#' nw_corner_distance(sort(x), sort(y), presorted = TRUE)
nw_corner_distance <- function(x, y, presorted = FALSE, p = 2, eps = sqrt(.Machine$double.eps)) {
  if (!presorted) {
    tx <- table(x)
    ty <- table(y)
    alpha <- as.numeric(tx) / sum(tx)
    beta  <- as.numeric(ty) / sum(ty)
    x <- as.numeric(names(tx))
    y <- as.numeric(names(ty))
  } else {
    # Run length encoding is in O(n) instead of table()'s O(n log n)
    rle_x <- rle(x)
    rle_y <- rle(y)
    alpha <- rle_x$lengths / sum(rle_x$lengths)
    beta <- rle_y$lengths / sum(rle_y$lengths)
    x <- rle_x$values
    y <- rle_y$values
  }

  n <- length(x)
  m <- length(y)

  k  <- 0

  total_distance <- 0

  i <- 1; j <- 1
  a <- alpha[1]; b <- beta[1]

  while (i <= n && j <= m) {
    t <- min(a, b)
    k <- k + 1
    dx <- abs(x[i] - y[j])
    total_distance <- total_distance + t * dx^p

    a <- a - t
    b <- b - t

    if (a <= eps) {
      i <- i + 1
      if (i <= n) a <- alpha[i]
    } else if (b <= eps) {
      j <- j + 1
      if (j <= m) b <- beta[j]
    }
  }
  total_distance^(1/p)
}
