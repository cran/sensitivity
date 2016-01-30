context("Comparison of soboljansen() to original version")

# # Original call:
# library(sensitivity)
# library(boot)
# set.seed(720)
# n <- 1000
# X1 <- data.frame(matrix(runif(8 * n), nrow = n))
# X2 <- data.frame(matrix(runif(8 * n), nrow = n))
# x_orig <- soboljansen(model = sobol.fun, X1, X2, nboot = 100)
# save(x_orig, file = "tests/testthat/soboljansen_original.RData")

load("soboljansen_original.RData")
library(boot)
set.seed(720)
n <- 1000
X1 <- data.frame(matrix(runif(8 * n), nrow = n))
X2 <- data.frame(matrix(runif(8 * n), nrow = n))
x <- soboljansen(model = sobol.fun, X1, X2, nboot = 100)

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
