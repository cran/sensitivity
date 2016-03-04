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

test_that("Models with matrix output work correctly (for OAT design)", {
  # A model function returning a matrix of two columns:
  ishigami.fun_matrix <- function(X){
    res_vector <- ishigami.fun(X)
    cbind(res_vector, 2 * res_vector)
  }
  set.seed(6923)
  x_oat_matrix <- morris(model = ishigami.fun_matrix, factors = 3, r = 20,
                         design = list(type = "oat", levels = 20, 
                                       grid.jump = 2))
  expect_equal(x_oat_matrix$ee[, , 1], x_oat_orig$ee, 
               check.attributes = FALSE)
  
  # A model function returning a matrix of only one column:
  ishigami.fun_onecol <- function(X){
    res_vector <- ishigami.fun(X)
    matrix(res_vector)
  }
  set.seed(6923)
  x_oat_onecol <- morris(model = ishigami.fun_onecol, factors = 3, r = 20,
                         design = list(type = "oat", levels = 20, 
                                       grid.jump = 2))
  expect_equal(x_oat_onecol$ee[, , 1], x_oat_orig$ee, 
               check.attributes = FALSE)
})

test_that("Models with matrix output work correctly (for simplex design)", {
  # A model function returning a matrix of two columns:
  ishigami.fun_matrix <- function(X){
    res_vector <- ishigami.fun(X)
    cbind(res_vector, 2 * res_vector)
  }
  set.seed(6923)
  # A dummy call to morris() with OAT design is made to ensure that the same
  # random numbers as in "x_simplex_orig" are drawn for the simplex design:
  oat_dummy <- morris(model = ishigami.fun, factors = 3, r = 20,
                      design = list(type = "oat", levels = 20, grid.jump = 2))
  x_simplex_matrix <- morris(model = ishigami.fun_matrix, factors = 3, r = 20,
                             design = list(type = "simplex", scale.factor = 1))
  expect_equal(x_simplex_matrix$ee[, , 1], x_simplex_orig$ee, 
               check.attributes = FALSE)
  
  # A model function returning a matrix of only one column:
  ishigami.fun_onecol <- function(X){
    res_vector <- ishigami.fun(X)
    matrix(res_vector)
  }
  set.seed(6923)
  oat_dummy <- morris(model = ishigami.fun, factors = 3, r = 20,
                      design = list(type = "oat", levels = 20, grid.jump = 2))
  x_simplex_onecol <- morris(model = ishigami.fun_onecol, factors = 3, r = 20,
                             design = list(type = "simplex", scale.factor = 1))
  expect_equal(x_simplex_onecol$ee[, , 1], x_simplex_orig$ee, 
               check.attributes = FALSE)
})

test_that("Models with array output work correctly (for OAT design)", {
  # A model function returning a 3-dimensional array with dimension 
  # lengths c(., 2, 2):
  ishigami.fun_array <- function(X){
    res_vector <- ishigami.fun(X)
    res_matrix <- cbind(res_vector, 2 * res_vector)
    array(data = c(res_matrix, 5 * res_matrix), 
          dim = c(length(res_vector), 2, 2))
  }
  set.seed(6923)
  x_oat_array <- morris(model = ishigami.fun_array, factors = 3, r = 20,
                        design = list(type = "oat", levels = 20, 
                                      grid.jump = 2))
  expect_equal(x_oat_array$ee[, , 1, 1], x_oat_orig$ee, 
               check.attributes = FALSE)
  
  # A model function returning a 3-dimensional array with dimension 
  # lengths c(., 2, 1):
  ishigami.fun_onedim3 <- function(X){
    res_vector <- ishigami.fun(X)
    res_matrix <- cbind(res_vector, 2 * res_vector)
    array(data = res_matrix, dim = c(length(res_vector), 2, 1))
  }
  set.seed(6923)
  x_oat_onedim3 <- morris(model = ishigami.fun_onedim3, factors = 3, r = 20,
                          design = list(type = "oat", levels = 20, 
                                        grid.jump = 2))
  expect_equal(x_oat_onedim3$ee[, , 1, 1], x_oat_orig$ee, 
               check.attributes = FALSE)
  
  # A model function returning a 3-dimensional array with dimension 
  # lengths c(., 1, 1):
  ishigami.fun_onecol_onedim3 <- function(X){
    res_vector <- ishigami.fun(X)
    res_matrix <- matrix(res_vector)
    array(data = res_matrix, dim = c(length(res_vector), 1, 1))
  }
  set.seed(6923)
  x_oat_onecol_onedim3 <- morris(model = ishigami.fun_onecol_onedim3, 
                                 factors = 3, r = 20,
                                 design = list(type = "oat", levels = 20, 
                                               grid.jump = 2))
  expect_equal(x_oat_onecol_onedim3$ee[, , 1, 1], x_oat_orig$ee, 
               check.attributes = FALSE)
  
  # A model function returning an array with more than 3 dimensions should
  # throw an error:
  ishigami.fun_dim4 <- function(X){
    res_vector <- ishigami.fun(X)
    res_matrix <- cbind(res_vector, 2 * res_vector)
    array(data = rep(c(res_matrix, 5 * res_matrix), 3), 
          dim = c(length(res_vector), 2, 2, 3))
  }
  expect_error(morris(model = ishigami.fun_dim4, factors = 3, r = 20,
                      design = list(type = "oat", levels = 20, 
                                    grid.jump = 2)),
               "If the model returns an array, it must not have")
  
  # A model function returning a list should throw an error:
  expect_error(morris(model = function(X) list(ishigami.fun(X)), factors = 3, 
                      r = 20, design = list(type = "oat", levels = 20, 
                                            grid.jump = 2)))
})

test_that("Models with array output work correctly (for simplex design)", {
  # A model function returning a 3-dimensional array with dimension 
  # lengths c(., 1, 1):
  ishigami.fun_onecol_onedim3 <- function(X){
    res_vector <- ishigami.fun(X)
    res_matrix <- matrix(res_vector)
    array(data = res_matrix, dim = c(length(res_vector), 1, 1))
  }
  set.seed(6923)
  oat_dummy <- morris(model = ishigami.fun, factors = 3, r = 20,
                      design = list(type = "oat", levels = 20, grid.jump = 2))
  x_simplex_onecol_onedim3 <- morris(model = ishigami.fun_onecol_onedim3, 
                                     factors = 3, r = 20,
                                     design = list(type = "simplex", 
                                                   scale.factor = 1))
  expect_equal(x_simplex_onecol_onedim3$ee[, , 1, 1], x_simplex_orig$ee, 
               check.attributes = FALSE)
  
  # A model function returning an array with more than 3 dimensions should
  # throw an error:
  ishigami.fun_dim4 <- function(X){
    res_vector <- ishigami.fun(X)
    res_matrix <- cbind(res_vector, 2 * res_vector)
    array(data = rep(c(res_matrix, 5 * res_matrix), 3), 
          dim = c(length(res_vector), 2, 2, 3))
  }
  expect_error(morris(model = ishigami.fun_dim4, factors = 3, r = 20,
                      design = list(type = "simplex", scale.factor = 1)),
               "If the model returns an array, it must not have")
  
  # A model function returning a list should throw an error:
  expect_error(morris(model = function(X) list(ishigami.fun(X)), factors = 3, 
                      r = 20, design = list(type = "simplex", 
                                            scale.factor = 1)))
})
