test_that("Matérn kernels work", {
  # Matérn with nu = 0.5
  expect_equal(matern(0.8, 1.2, 0.6, 0.5, squared = FALSE),
               0.6*exp(-0.8/1.2))
  expect_equal(matern(0.8, 1.2, 0.6, 0.5, squared = TRUE),
               0.6*exp(-sqrt(0.8)/1.2))
  # Matérn with nu = 1.5
  expect_equal(matern(0.8, 1.2, 0.6, 1.5, squared = FALSE),
               0.6*(1 + sqrt(3)*0.8/1.2)*exp(-sqrt(3)*0.8/1.2))
  expect_equal(matern(0.8, 1.2, 0.6, 1.5, squared = TRUE),
               0.6*(1 + sqrt(3)*sqrt(0.8)/1.2)*exp(-sqrt(3)*sqrt(0.8)/1.2))
  # Matérn with nu = 2.5
  expect_equal(matern(0.8, 1.2, 0.6, 2.5, squared = FALSE),
               0.6*(1 + sqrt(5)*0.8/1.2 + 5*0.8^2 / (3*1.2^2))*exp(-sqrt(5)*0.8/1.2))
  expect_equal(matern(0.8, 1.2, 0.6, 2.5, squared = TRUE),
               0.6*(1 + sqrt(5)*sqrt(0.8)/1.2 + 5*0.8 / (3*1.2^2))*exp(-sqrt(5)*sqrt(0.8)/1.2))
  # If distance = 0, output variance
  expect_equal(matern(0.00, 0.5, 12.1, 2.3), 12.1)
  M <- rbind(c(0, 2, 0),
             c(2, 0, 2.5),
             c(0, 2.5, 0))
  K <- matern(M, nu = 2.3)
  expect_equal(dim(K), dim(M))
  expect_equal(diag(K), rep(1, nrow(M)))
})

test_that("Sum of Matérn kernels works", {
  M1 <- rbind(c(0, 2, 1),
              c(2, 0, 2.5),
              c(1, 2.5, 0))
  M2 <- rbind(c(0, 1, 4),
              c(1, 0, 3),
              c(4, 3, 0))
  M_list <- list(M1, M2)
  Matern_sum <- matern_multiple(M_list, length_scales = c(1, 2), variances = c(2, 1/2),
                                nus = c(0.5, 2.5), squared = FALSE)
  expect_equal(Matern_sum, matern(M1, 1, 2, 0.5, FALSE) + matern(M2, 2, 1/2, 2.5, FALSE))
})
