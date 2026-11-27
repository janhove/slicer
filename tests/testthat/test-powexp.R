test_that("Power exponential works", {
  M <- rbind(c(0, 2, 1),
             c(2, 0, 2.5),
             c(1, 2.5, 0))
  # For expon = 1, identical to Matérn 1/2.
  expect_equal(powexp(M, expon = 1, squared = FALSE),
               matern(M, nu = 0.5, squared = FALSE))
  # For expon = 2, identical to RBF with length-scale t/sqrt(2)
  expect_equal(powexp(M, expon = 2, squared = FALSE),
               rbf(M, length_scale = 1/sqrt(2), squared = FALSE))
})
