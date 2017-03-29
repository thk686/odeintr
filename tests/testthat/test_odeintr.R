library(odeintr)
context("Solvers")

Sys.setenv(R_TESTS = "")

# Try to find package header file *somewhere*
wd = getwd()
wd = sub(".tests.testthat$", "", wd)
ipath1 = file.path(wd, "inst", "include")
ipath2 = file.path(wd, "include")
ipath3 = system.file("include", package = "odeintr")
Sys.setenv(PKG_CXXFLAGS = paste(paste0("-I\"", ipath1, "\""),
                                paste0("-I\"", ipath2, "\""),
                                if (nzchar(ipath3)) paste0("-I\"", ipath3, "\"")))

if (file.exists(file.path(ipath1, "odeintr.h")) ||
    file.exists(file.path(ipath2, "odeintr.h")) ||
    file.exists(file.path(ipath3, "odeintr.h")))
{

test_that("integrate_sys works", {
  res = integrate_sys(function(x, t) x * (1 - x), 0.01, 40)
  expect_is(res, "data.frame")
  expect_equal(dim(res), c(41, 2))
  expect_equal(res[1, 2], 0.01)
})

test_that("compile_sys works", {
  compile_sys("logistic", "dxdt = x * (1 - x)")
  res = logistic(0.01, 40)
  expect_is(res, "data.frame")
  expect_equal(dim(res), c(41, 2))
  expect_equal(res[1, 2], 0.01)
})

test_that("large step size does not crash compile_sys", {
  compile_sys("logi", "dxdt[0] = x[0] * (1 - x[0])")
  expect_warning({res = logi(0.01, 40, 100)})
  expect_is(res, "data.frame")
  expect_equal(dim(res), c(4, 2))
  expect_equal(res[1, 2], 0.01)
})

test_that("large step size does not crash compile_implicit", {
  compile_implicit("logi", "dxdt[0] = x[0] * (1 - x[0])")
  expect_warning({res = logi(0.01, 40, 100)})
  expect_is(res, "data.frame")
  expect_equal(dim(res), c(4, 2))
  expect_equal(res[1, 2], 0.01)
})

}