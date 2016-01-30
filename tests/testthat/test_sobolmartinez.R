context("Comparison of sobolmartinez() to original version")

# # Original call:
# library(sensitivity)
# library(boot)
# set.seed(1025)
# n <- 1000
# X1 <- data.frame(matrix(runif(8 * n), nrow = n))
# X2 <- data.frame(matrix(runif(8 * n), nrow = n))
# x_orig <- sobolmartinez(model = sobol.fun, X1, X2, nboot = 0)
# save(x_orig, file = "tests/testthat/sobolmartinez_original.RData")

load("sobolmartinez_original.RData")
library(boot)
set.seed(1025)
n <- 1000
X1 <- data.frame(matrix(runif(8 * n), nrow = n))
X2 <- data.frame(matrix(runif(8 * n), nrow = n))
x <- sobolmartinez(model = sobol.fun, X1, X2, nboot = 0)

test_that("Calls are equal", {
  expect_equal(x$call, x_orig$call)
})

test_that("Random matrices (X1, X2, X) are equal", {
  expect_equal(x$X1, x_orig$X1)
  expect_equal(x$X2, x_orig$X2)
  expect_equal(x$X, x_orig$X)
})

test_that("Response vectors (y) are equal", {
  expect_equal(x$y, x_orig$y)
})

test_that("SA indices (S and T) are equal", {
  expect_equal(x$S, x_orig$S)
  expect_equal(x$T, x_orig$T)
})

test_that("VCEs (V) are equal", {
  expect_equal(x$V, x_orig$V)
})
