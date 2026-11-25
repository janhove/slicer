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

test_that("Sum of power exponential works", {
  M1 <- rbind(c(0, 2, 1),
              c(2, 0, 2.5),
              c(1, 2.5, 0))
  M2 <- rbind(c(0, 1, 4),
              c(1, 0, 3),
              c(4, 3, 0))
  M_list <- list(M1, M2)
  Powexp_sum <- powexp_multiple(M_list, length_scales = c(1, 2), variances = c(2, 1/2),
                                expons = c(0.5, 1.5), squared = FALSE)
  expect_equal(Powexp_sum, powexp(M1, 1, 2, 0.5, FALSE) + powexp(M2, 2, 1/2, 1.5, FALSE))
})
