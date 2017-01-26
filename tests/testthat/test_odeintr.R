library(odeintr)
context("Solvers")

test_that("integrate_sys works", {
  res = integrate_sys(function(x, t) x * (1 - x), 0.01, 40)
  expect_is(res, "data.frame")
  expect_equal(dim(res), c(41, 2))
  expect_equal(res[1, 2], 0.01)
})

# test_that("compile_sys works", {
#   compile_sys("logistic", "dxdt = x * (1 - x)")
#   res = logistic(0.01, 40)
#   expect_is(res, "data.frame")
#   expect_equal(dim(res), c(41, 2))
#   expect_equal(res[1, 2], 0.01)
# })
# 
# test_that("compile_implicit works", {
#   compile_implicit("logi", "dxdt = x * (1 - x)")
#   res = logi(0.01, 40)
#   expect_is(res, "data.frame")
#   expect_equal(dim(res), c(41, 2))
#   expect_equal(res[1, 2], 0.01)
# })
