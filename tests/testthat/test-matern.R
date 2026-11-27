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
