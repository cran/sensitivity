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

test_that("Models with matrix output work correctly", {
  # A model function returning a matrix of two columns:
  sobol.fun_matrix <- function(X){
    res_vector <- sobol.fun(X)
    cbind(res_vector, 2 * res_vector)
  }
  x_matrix <- soboljansen(model = sobol.fun_matrix, X1, X2)
  expect_equal(x_matrix$V[, 1], x_orig$V[, "original"], 
               check.attributes = FALSE)
  expect_equal(x_matrix$S[, 1], x_orig$S[, "original"], 
               check.attributes = FALSE)
  expect_equal(x_matrix$T[, 1], x_orig$T[, "original"], 
               check.attributes = FALSE)
  
  # A model function returning a matrix of only one column:
  sobol.fun_onecol <- function(X){
    res_vector <- sobol.fun(X)
    matrix(res_vector)
  }
  x_onecol <- soboljansen(model = sobol.fun_onecol, X1, X2)
  expect_equal(drop(x_onecol$V), x_orig$V[, "original"], 
               check.attributes = FALSE)
  expect_equal(drop(x_onecol$S), x_orig$S[, "original"], 
               check.attributes = FALSE)
  expect_equal(drop(x_onecol$T), x_orig$T[, "original"], 
               check.attributes = FALSE)
  
  # Trying to perform a bootstrap should fail:
  expect_error(soboljansen(model = sobol.fun_matrix, X1, X2, nboot = 100))
})

test_that("Models with array output work correctly", {
  # A model function returning a 3-dimensional array with dimension 
  # lengths c(., 2, 2):
  sobol.fun_array <- function(X){
    res_vector <- sobol.fun(X)
    res_matrix <- cbind(res_vector, 2 * res_vector)
    array(data = c(res_matrix, 5 * res_matrix), 
          dim = c(length(res_vector), 2, 2))
  }
  x_array <- soboljansen(model = sobol.fun_array, X1, X2)
  expect_equal(x_array$V[, 1, 1], x_orig$V[, "original"], 
               check.attributes = FALSE)
  expect_equal(x_array$S[, 1, 1], x_orig$S[, "original"], 
               check.attributes = FALSE)
  expect_equal(x_array$T[, 1, 1], x_orig$T[, "original"], 
               check.attributes = FALSE)
  
  # A model function returning a 3-dimensional array with dimension 
  # lengths c(., 2, 1):
  sobol.fun_onedim3 <- function(X){
    res_vector <- sobol.fun(X)
    res_matrix <- cbind(res_vector, 2 * res_vector)
    array(data = res_matrix, dim = c(length(res_vector), 2, 1))
  }
  x_onedim3 <- soboljansen(model = sobol.fun_onedim3, X1, X2)
  expect_equal(x_onedim3$V[, 1, 1], x_orig$V[, "original"], 
               check.attributes = FALSE)
  expect_equal(x_onedim3$S[, 1, 1], x_orig$S[, "original"], 
               check.attributes = FALSE)
  expect_equal(x_onedim3$T[, 1, 1], x_orig$T[, "original"], 
               check.attributes = FALSE)
  
  # A model function returning a 3-dimensional array with dimension 
  # lengths c(., 1, 1):
  sobol.fun_onecol_onedim3 <- function(X){
    res_vector <- sobol.fun(X)
    res_matrix <- matrix(res_vector)
    array(data = res_matrix, dim = c(length(res_vector), 1, 1))
  }
  x_onecol_onedim3 <- soboljansen(model = sobol.fun_onecol_onedim3, X1, X2)
  expect_equal(x_onecol_onedim3$V[, 1, 1], x_orig$V[, "original"], 
               check.attributes = FALSE)
  expect_equal(x_onecol_onedim3$S[, 1, 1], x_orig$S[, "original"], 
               check.attributes = FALSE)
  expect_equal(x_onecol_onedim3$T[, 1, 1], x_orig$T[, "original"], 
               check.attributes = FALSE)
  
  # A model function returning an array with more than 3 dimensions should
  # throw an error:
  sobol.fun_dim4 <- function(X){
    res_vector <- sobol.fun(X)
    res_matrix <- cbind(res_vector, 2 * res_vector)
    array(data = rep(c(res_matrix, 5 * res_matrix), 3), 
          dim = c(length(res_vector), 2, 2, 3))
  }
  expect_error(soboljansen(model = sobol.fun_dim4, X1, X2))
  
  # A model function returning a list should throw an error:
  expect_error(soboljansen(model = function(X) list(sobol.fun(X)), X1, X2))
  
  # Trying to perform a bootstrap should fail:
  expect_error(soboljansen(model = sobol.fun_array, X1, X2, nboot = 100))
})
