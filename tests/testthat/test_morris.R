context("Comparison of morris() to original version")

# # Original call:
# library(sensitivity)
# set.seed(6923)
# x_oat_orig <- morris(model = ishigami.fun, factors = 3, r = 20,
#                      design = list(type = "oat", levels = 20, grid.jump = 2))
# x_simplex_orig <- morris(model = ishigami.fun, factors = 3, r = 20,
#                          design = list(type = "simplex", scale.factor = 1))
# save(x_oat_orig, x_simplex_orig, file = "tests/testthat/morris_original.RData")

load("morris_original.RData")
set.seed(6923)
x_oat <- morris(model = ishigami.fun, factors = 3, r = 20,
                design = list(type = "oat", levels = 20, grid.jump = 2))
x_simplex <- morris(model = ishigami.fun, factors = 3, r = 20,
                    design = list(type = "simplex", scale.factor = 1))

test_that("Calls are equal", {
  expect_equal(x_oat$call, x_oat_orig$call)
})

test_that("Random matrices (X) are equal", {
  expect_equal(x_oat$X, x_oat_orig$X)
  expect_equal(x_simplex$X, x_simplex_orig$X)
})

test_that("Response vectors (y) are equal", {
  expect_equal(x_oat$y, x_oat_orig$y)
  expect_equal(x_simplex$y, x_simplex_orig$y)
})

test_that("Elementary effects (ee) are equal", {
  expect_equal(x_oat$ee, x_oat_orig$ee)
  expect_equal(x_simplex$ee, x_simplex_orig$ee)
})
