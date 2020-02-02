library(testthat)
source("CorrAggBias.R")


# Fisher Z ----

test_that("Fisher Z", {
  # Use of an alternative algorithm
  inp <- seq(-1, 1, 0.01)
  expect_equal(atanh(inp), FisherZ(inp))
  
  # Taken from the SAS manual
  inp <- 0.22081
  expect_equal(0.22451, round(FisherZ(0.22081)), 5)
  # Taken from https://www.statisticshowto.datasciencecentral.com/fisher-z/  
  expect_equal(0.0000, round(FisherZ(0.0000)), 4)
  expect_equal(0.0100, round(FisherZ(0.0100)), 4)
  expect_equal(0.0200, round(FisherZ(0.0200)), 4)
  expect_equal(0.0300, round(FisherZ(0.0300)), 4)
  expect_equal(0.0400, round(FisherZ(0.0400)), 4)
  expect_equal(0.0500, round(FisherZ(0.0500)), 4)
  expect_equal(0.0601, round(FisherZ(0.0600)), 4)
  #...
  expect_equal(0.1614, round(FisherZ(0.1600)), 4)
  #...
  expect_equal(2.2976, round(FisherZ(0.9800)), 4)
  expect_equal(2.6467, round(FisherZ(0.9900)), 4)
  
  # Taken from https://rdrr.io/cran/GeneNet/man/z.transform.html
  inp <- c(-0.26074194, 0.47251437, 0.23957283,-0.02187209,-0.07699437,
         -0.03809433,-0.06010493, 0.01334491,-0.42383367,-0.25513041)
  e <- c(-0.26690430, 0.51330253, 0.24432088, -0.02187558, -0.07714706, 
         -0.03811277, -0.06017747, 0.01334570, -0.45235595, -0.26089280)
  expect_equal(e, FisherZ(inp))
  
  # not defined for r = 1
  expect_identical(Inf, FisherZ(1))
  expect_identical(-Inf, FisherZ(-1))
  
})


test_that("Fisher Z Inverse", {
  # Use of an alternative algorithm
  inp <- seq(-1, 1, 0.01)
  expect_equal(tanh(inp), FisherZInv(inp))
  
  # Taken from https://rdrr.io/cran/GeneNet/man/z.transform.html
  e <- c(-0.26074194, 0.47251437, 0.23957283,-0.02187209,-0.07699437,
           -0.03809433,-0.06010493, 0.01334491,-0.42383367,-0.25513041)
  inp <- c(-0.26690430, 0.51330253, 0.24432088, -0.02187558, -0.07714706, 
         -0.03811277, -0.06017747, 0.01334570, -0.45235595, -0.26089280)
  expect_equal(e, FisherZInv(inp))
})


test_that("Fisher Z: Reversion Tests", {
  inp <- seq(-1, 1, 0.01)
  expect_equal(inp, FisherZInv(FisherZ(inp)))
  expect_equal(inp, FisherZ(FisherZInv(inp)))
})


# Hotelling Z ----

test_that("Hotelling Z", {
  # Use of an alternative algorithm
  inp <- seq(-1, 1, 0.01)
  
  # Taken from https://rdrr.io/cran/GeneNet/man/z.transform.html
  inp <- c(-0.26074194, 0.47251437, 0.23957283,-0.02187209,-0.07699437,
         -0.03809433,-0.06010493, 0.01334491,-0.42383367,-0.25513041)
  df <- 7
  e <- c(-0.22899520, 0.44143031, 0.20958747, -0.01875062, -0.06613150, 
         -0.03266875, -0.05158328, 0.01143920, -0.38875232, -0.22382820)
  expect_equal(e, HotellingZ(inp, df))
  
  # 
  expect_error(HotellingZ(0), "Degrees of freedom 'df' are missing")
})


test_that("Hotelling Z Inverse", {
  # Use of an alternative algorithm
  inp <- seq(-1, 1, 0.01)
  #expect_equal(tanh(inp), FisherZInv(inp))
  
  # Taken from https://rdrr.io/cran/GeneNet/man/z.transform.html
  e <- c(-0.26074194, 0.47251437, 0.23957283,-0.02187209,-0.07699437,
           -0.03809433,-0.06010493, 0.01334491,-0.42383367,-0.25513041)
  df <- 7
  inp <- c(-0.22899520, 0.44143031, 0.20958747, -0.01875062, -0.06613150, 
         -0.03266875, -0.05158328, 0.01143920, -0.38875232, -0.22382820)
  expect_equal(e, HotellingZInv(inp, e, df))

  # Error: missing df
  expect_error(HotellingZInv(0), "Untransformed correlations 'r' are missing")
  expect_error(HotellingZInv(0, 2), "Degrees of freedom 'df' are missing")
  
})



test_that("Hotelling Z: Reversion Test", {
  # Use of an alternative algorithm
  inp <- seq(-0.99, 0.99, 0.01)
  df <- 2
  expect_equal(inp, HotellingZInv(HotellingZ(inp, df), inp, df))
  # 
})





# MinVar G by Olkin & Pratt G ----

test_that("MinVar G", {
  # Olin & Pratt (1958) - 
  # tolerance had to be lowered because the paper only gives 3 digits
  # Small sample sizes do not work because they do not report values computed
  # by the approximation but by recursive iteration.
  inp.r <- rep(0.5, 3)
  n <- c(11, 21, 31)          # c(5, 11, 21, 31)
  e <- c(0.525, 0.511, 0.507) # c(0.605, 0.525, 0.511, 0.507)
  expect_equal(MinVarZ.hg(inp.r, n), e, tolerance = 1e-4)
  
  inp.r <- rep(0.5, 4)
  n <- c(5, 11, 21, 31)
  e <- c(0.605, 0.525, 0.511, 0.507) 
  expect_equal(MinVarZ.hg(inp.r, n), e, tolerance = 1e-4)
  

  inp.r <- rep(0.1, 3)
  n <- c(11, 21, 31) # c(5, 7, 11, 21, 31)
  e <- c(0.107, 0.103, 0.102) # c(0.148, 0.117, 0.107, 0.103, 0.102)
  expect_equal(MinVarZ(inp.r, n), e, tolerance = 1e-3)
  
  # Taken from https://www.psychometrica.de/correlation.html#fisher
  inp.r <- 0.29
  n <- 45
  e <- 0.2932
  expect_equal(MinVarZ(inp.r, n), e, tolerance = 1e-4)

})

# CorrAggBias.test.R:124: failure: MinVar G
# MinVarZ(inp.r, n) not equal to `e`.
# 1/3 mismatches
# [1] 0.688 - 0.605 == 0.0825
# 
# CorrAggBias.test.R:129: failure: MinVar G
# MinVarZ(inp.r, n) not equal to `e`.
# 1/5 mismatches
# [1] 0.15 - 0.148 == 0.0015




# MeanR ----
test_that("MeanR", {
  # Default averaging
  inp <- seq(-0.5, 0.5, 0.5)
  expect_equal(0, MeanR(inp, 2, "Default"))
  expect_equal(0, MeanR(inp, 10, "Default"))
  expect_equal(0, MeanR(inp, c(2, 5, 2), "Default"))
  
  # Fisher averaging
  inp <- seq(-0.5, 0.5, 0.5)
  expect_equal(0, MeanR(inp, 2, "Fisher"))
  expect_equal(0, MeanR(inp, 10, "Fisher"))
  expect_equal(0, MeanR(inp, c(2, 5, 2), "Fisher"))
  
  # Unequal n
  inp.r <- c(0.2900, 0.4500, 0.3600)
  inp.n <- c(45, 28, 32)
  e <- 0.3558 # ??? 0.35580172832960727
  #expect_equal(MeanR(inp.r, inp.n, "Fisher"), e)
  e <- 1.065834587108014 / 3 #???
  #expect_equal(MeanR(inp.r, inp.n, "MinVar"), e)
  

  # Hotelling averaging
  inp <- seq(-0.5, 0.5, 0.5)
  # expect_equal(0, MeanR(inp, 2, "H"))
  # expect_equal(0, MeanR(inp, 10, "Hotelling"))
  # expect_equal(0, MeanR(inp, c(2, 5, 2), "Hotelling"))
  
  # Check errors
  inp <- seq(-0.5, 0.5, 0.5)
  expect_error(MeanR(inp, 1, "Default"),
               "Sample size 'N' must be larger than 1")
  expect_error(MeanR(inp, 1, "Fisher"),
               "Sample size 'N' must be larger than 1")
  expect_error(MeanR(inp, 1, "Hotelling"),
               "Sample size 'N' must be larger than 1")
})


