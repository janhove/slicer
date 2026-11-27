test_that("summing kernels works", {
  M1 <- rbind(c(0, 2, 1),
              c(2, 0, 2.5),
              c(1, 2.5, 0))
  M2 <- rbind(c(0, 1, 4),
              c(1, 0, 3),
              c(4, 3, 0))
  M3 <- rbind(c(0, 4, 8),
              c(4, 0, 1),
              c(8, 1, 0))
  M_list <- list(M1^2, M2^2, M3^2)
  K_sum <- kernel_sum(M_list, length_scales = c(1, 2, 0.5),
                      variances = c(2, 1/2, 0.1),
                      kernels = c("powexp", "matern", "rbf"),
                      params = c(1, 2.3, NA),
                      squared = TRUE)
  expect_equal(
    K_sum,
    powexp(M1^2, 1, 2, 1, squared = TRUE) +
      matern(M2^2, 2, 1/2, 2.3, squared = TRUE) +
      rbf(M3^2, 0.5, 0.1, squared = TRUE)
  )
})
