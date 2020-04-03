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

  # Taken from https://rdrr.io/cran/GeneNet/man/z.transform.html
  e <- c(-0.26074194, 0.47251437, 0.23957283,-0.02187209,-0.07699437,
           -0.03809433,-0.06010493, 0.01334491,-0.42383367,-0.25513041)
  inp <- c(-0.22899520, 0.44143031, 0.20958747, -0.01875062, -0.06613150, 
         -0.03266875, -0.05158328, 0.01143920, -0.38875232, -0.22382820)
  df <- 7
  
  for(i in 1:10) {
    # Single number input
    expect_equal(e[i], HotellingZInv(inp[i], df),
                 info = paste("i =", i))
  }
  
  # Vector input
  expect_equal(e, HotellingZInv(inp, df),
               info = paste("df =", df))

  # Error: missing df
  expect_error(HotellingZInv(0), "Degrees of freedom 'df' are missing")
})


test_that("Hotelling Z: Reversion Test", {
  # Use of an alternative algorithm
  inp <- seq(-0.99, 0.99, 0.01)
  for (df in c(2,5, 10, 50, 200, 2000))
    expect_equal(inp, HotellingZInv(HotellingZ(inp, df), df))
  # 
})


# Hotellings 2. Equation ----
test_that("Hotelling2", {
  r <- seq(-0.9, 0.9, 0.1)
  z.ori <- seq(-0.9, 0.9, 0.1)
  z <- z.ori
  df <- 5
  zh <- z - (3*z+r)/(4*df) - ((23*z + 33*r - 5*r^3) / (96 * df^2))
  
  #z <- ( 5*df^2*r^3 - 33*df^2*r - 24*df*r ) / ( 23*df^2 + 72*df - 96 )
  #z = (df*r* (df*(5*r^2 - 33) - 24) - 96*zh) / (23*df^2 + 72*df - 96) #and n*(23*n + 72)!=96
  #z <- (96*zh + 24*r*df - 33*r + 5*r^3) / (24*df * (4*df-3) - 23)
  
  # solve o = z - (3*z+r)/(4*n) - (23*z + 33*r - 5*r^3) / (96 * n^2)  for z
  #z <- (24*df*r + 96*df^2*zh - 5*r^3 + 33*r) / ( - 72*df^3 + 96*df^2 - 23 )
  

  
  # Taken from https://rdrr.io/cran/GeneNet/man/z.transform.html
  # These values work for Hotelling - for Hotelling2 they will be imprecise
  inp <- c(-0.26074194, 0.47251437, 0.23957283,-0.02187209,-0.07699437,
           -0.03809433,-0.06010493, 0.01334491,-0.42383367,-0.25513041)
  df <- 7
  e <- c(-0.22899520, 0.44143031, 0.20958747, -0.01875062, -0.06613150, 
         -0.03266875, -0.05158328, 0.01143920, -0.38875232, -0.22382820)
  expect_equal(e, HotellingZ2(inp, df), tolerance = .006)

  
  inp <- seq(-0.99, 0.99, 0.01)
  df <- 2
  #expect_equal(inp, HotellingZ2Inv(HotellingZ2(inp, df), inp, df))

  inp <- seq(-0.99, 0.99, 0.01)
  df <- 10
  #expect_equal(inp, HotellingZ2Inv(HotellingZ2(inp, df), inp, df))

  inp <- seq(-0.99, 0.99, 0.01)
  df <- 20
  #expect_equal(inp, HotellingZ2Inv(HotellingZ2(inp, df), inp, df))
  
})


# MinVar (approx. G by Olkin & Pratt, 1958)  ----

test_that("MinVar", {
  # Olin & Pratt (1958) - 
  # tolerance had to be lowered because the paper only gives 3 digits
  # Small sample sizes do not work because they do not report values computed
  # by the approximation but by recursive iteration.
  inp.r <- rep(0.5, 3)
  n <- c(11, 21, 31)
  e <- c(0.525, 0.511, 0.507)
  expect_equal(MinVarZ(inp.r, n), e, tolerance = 1e-3)
  
  inp.r <- rep(0.5, 3)
  n <- c(13, 21, 31)
  e <- c(0.520, 0.511, 0.507)
  expect_equal(MinVarZ(inp.r, n), e, tolerance = 1e-3)
  

  inp.r <- rep(0.1, 3)
  n <- c(11, 21, 31)
  e <- c(0.107, 0.103, 0.102) 
  expect_equal(MinVarZ(inp.r, n), e, tolerance = 1e-3)
  
  # Taken from https://www.psychometrica.de/correlation.html#fisher
  inp.r <- 0.29
  n <- 45
  e <- 0.2932
  expect_equal(MinVarZ(inp.r, n), e, tolerance = 1e-4)

})



# MinVar true k (approx. G by Olkin & Pratt, 1958) ----

test_that("MinVar TrueK", {
  # No data available to test anything
  expect_equal(0, 0)
})





# Precise MinVar G by Olkin & Pratt G ----

test_that("MinVar G Precise", {
  MinVarZ.pr.alt <- function(r, n) {
    fc2 <- function(t, r, df) {
      ( t^(-0.5) * (1+t)^(1-0.5*df) ) / (1+t*r^2)^0.5 
    }
    Fc2 <- function(r, df) integrate(fc2, 0, Inf, df = df, r = r)$value
    
    if(any(n < 5)) stop("Sample size must be greater than 3")
    df <- n-1 # if all µ and σ are unkonwn (if µ was know it'd be n)
    dfh <- df/2
    # Use value for gamma(0.5) to save computation time - https://oeis.org/A002161
    GammaOfHalf <- 1.772453850905516027298167483341145182797549456122387128213807789852911284591032181374950656738544665
    Component1 <- gamma(dfh-0.5) / GammaOfHalf / gamma(dfh-1)
    Component2 <- mapply(Fc2, df = df, r = r)
    G <- r * Component1 * Component2
    
    return(G)
  }
  
  # This function throws a lot of warnings
  library(hypergeo)
  MinVarZ.hg.udf <- function(r, n) {
    if(any(n < 5)) stop("Sample size must be greater than 3")
    df <- n-1
    G <- r * hypergeo(0.5, 0.5, (df-1)/2, 1-r^2)
    return(Re(G))
  }
  
  # Compare function to alternative form provided by Olkin & Pratt, 1958
  for(r in c(0.01, 0.1, 0.3, 0.5, 0.75, 0.9, 0.99))
    for(n in c(5, 10, 25, 50, 100))
      expect_equal(MinVarZ.pr(r, n), MinVarZ.pr.alt(r, n))    
  for(r in c(0.01, 0.1, 0.3, 0.5, 0.75, 0.9, 0.99) * -1)
    for(n in c(5, 10, 25, 50, 100))
      expect_equal(MinVarZ.pr(r, n), MinVarZ.pr.alt(r, n))    
  
  # Compare function to alternative form provided by Olkin & Pratt, 1958
  # using hypergeometric function
  # r below 0.11 will not converge and throw an error
  for(r in c(0.11, 0.3, 0.5, 0.75, 0.9, 0.99))
    for(n in c(5, 10, 25, 50, 100))
      expect_equal(MinVarZ.pr(r, n), MinVarZ.hg.udf(r, n), tolerance=1e-7)#,
  for(r in c(0.11, 0.3, 0.5, 0.75, 0.9, 0.99) * -1)
    for(n in c(5, 10, 25, 50, 100))
      expect_equal(MinVarZ.pr(r, n), MinVarZ.hg.udf(r, n), tolerance=1e-7)#,
  
  # Test vectorized function call
  r <- seq(0.01, 0.99, 0.01)
  n <- rep(seq(5, 50, 5), length.out = length(r))
  expect_equal(MinVarZ.pr(r, n), MinVarZ.pr.alt(r, n))

  # Test vectorized function call
  r <- seq(0.11, 0.99, 0.01)
  n <- rep(seq(5, 50, 5), length.out = length(r))
  expect_equal(MinVarZ.pr(r, n), MinVarZ.hg.udf(r, n))
  
  # Use particular values from Olkin & Pratt (1958)
  # Tolerance is reduced because the paper provides only 3 digits
  r <- seq(0.1, 0.9, 0.1)
  n <- rep(5, length.out = length(r))
  e <- c(.148, .280, .398, .506, .605, .695, .780, .858, .931)
  expect_equal(MinVarZ.pr(r, n), e, tolerance=5e-4)
  r <- seq(0.1, 0.9, 0.1)
  n <- rep(15, length.out = length(r))
  e <- c(.105, .209, .312, .415, .516, .616, .715, .812, .907)
    expect_equal(MinVarZ.pr(r, n), e, tolerance=5e-4)
  r <- seq(0.1, 0.9, 0.1)
  n <- rep(31, length.out = length(r))
  e <- c(.102, .204, .305, .406, .507, .607, .706, .805, .903)
    expect_equal(MinVarZ.pr(r, n), e, tolerance=5e-4)
})
  





# MeanR ----
test_that("MeanR", {
  # Default averaging
  inp <- seq(-0.5, 0.5, 0.5)
  expect_equal(0, MeanR(inp, 2, "None"))
  expect_equal(0, MeanR(inp, 10, "None"))
  expect_equal(0, MeanR(inp, c(2, 5, 2), "None"))
  #TODO: MORE TESTS NEEDED FOR 'NONE'
  
  
  # Fisher averaging
  inp <- seq(-0.5, 0.5, 0.5)
  expect_equal(0, MeanR(inp, 2, "Fish"))
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
  expect_equal(MeanR(inp, 4, "Hotelling"), 0)
  expect_equal(MeanR(inp, 10, "Hotelling"), 0)
  expect_equal(MeanR(inp, c(4, 10, 4), "Hotelling"), 0)
  
  
  # averaging sgn(r)*r²
  inp <- seq(-0.5, 0.5, 0.1)
  expect_equal(MeanR(inp, 2, "Squared"), 0) #0.31622776601683794)
  expect_equal(MeanR(inp, 20, "Squared"), 0) #0.31622776601683794)
  expect_equal(MeanR(inp, 200, "Squared"), 0) #0.31622776601683794)
  
  inp <- seq(0.0, 1, 0.1)
  expect_equal(MeanR(inp, 2, "Squared"), sqrt(0.35))
  
  inp <- 1
  expect_equal(MeanR(inp, 2, "Squared"), 1)
  
  
  # Check errors
  inp <- seq(-0.5, 0.5, 0.5)
  expect_error(MeanR(inp, 1, "None"),
               "Sample size 'N' must be larger than 1")
  expect_error(MeanR(inp, 1, "Fisher"),
               "Sample size 'N' must be larger than 1")
  expect_error(MeanR(inp, 1, "Hotelling"),
               "Sample size 'N' must be larger than 1")

  expect_error(MeanR(inp, 0, "Squared"),
               "Sample size 'N' must be larger than 1")
  
  expect_error(MeanR(inp, 1, "Hotel"),
               "'arg' should be one of")
  expect_error(MeanR(inp, 1, "TotalerQuatsch"),
               "'arg' should be one of")
  
})


