library(hypergeo)

#' FisherZ
#' Fisher z transformation or it's inverse on a vector 
#' of correlation coefficients.
#' @param r A vectors of correlation coefficients
#' @param z A vectors of fisher/hotelling z transformed correlation 
#' coefficients
#' @param df degrees of freedom of the distribution of the 
#' correlation coefficient
#' @param n Sample size
#' @value Hotelling z is NaN for df <= 1
#' @note MinVarZ.pr is defined for n < 600 (roughly).
#' @author Jan Seifert
FisherZ <- function( r ) {
  # PRECONDITIONS
  if(missing(r)) return(NA)
  if(any(r > 1) || any(r < -1)) 
    stop("Fisher Z transformation is only defined for correlations (i.e. -1 <= r <= 1)")
  
  # Identify +/-1 because transformation cannot handle that
  rp1 <- which(r == 1)
  rm1 <- which(r == -1)
  # Transform
  z <- log( (1+r) / (1-r) ) / 2
  # Set r= 1/-1 to Inf/-Inf
  z[rp1] <- Inf
  z[rm1] <- -Inf
  #
  return( z )
}

#' FisherZInv
#' @describeIn FisherZ Inverse of the fisher z transformation
FisherZInv <- function( z ) {
  # PRECONDITIONS
  if(missing(z)) return(NA)
  if(any(is.nan(z))) z[which(is.nan(z))] <- 1 # 'FisherZ()' not defined for 1
  
  zp1 <- which(z == Inf) # z == plus 1
  zm1 <- which(z ==-Inf) # z == minus 1
  # Tranform inversely
  r <- (exp(2*z)-1) / (exp(2*z)+1)
  # Set r = 1/-1 to Inf/-Inf
  r[zp1] <- 1
  r[zm1] <- -1
  return(r)
}

#' HotellingZ
#' @describeIn FisherZ Improved fisher z transformation by Hotelling (1953)
#' @references Hotelling, H. (1953). New light on the correlation coefficient and
#' its transformations. Journal of the Royal Statistical Society, 15, p. 193-225.
HotellingZ <- function( r, df ) {
  if(missing(df)) stop("Degrees of freedom 'df' are missing")
  if(any(df <= 1)) dfleq1 <- which(df <= 1)
  
  z = FisherZ( r )
  if(any(df <= 1)) z[dfleq1] <- NaN
  return( z - (3*z+r)/(4*df) )
}

#' HotellingZ
#' @describeIn FisherZ Inverse for \code{HotellingZ}
HotellingZInv <- function( z, r, df ) {
  if(missing(r)) stop("Untransformed correlations 'r' are missing")
  if(missing(df)) stop("Degrees of freedom 'df' are missing")
  if(any(df <= 1)) dfleq1 <- which(df <= 1)
  
  Fisher <- (4*df * z + r) / (4*df - 3)
  if(any(df <= 1)) Fisher[dfleq1] <- NaN
  return(FisherZInv(Fisher))
}


#' MinVarZ
#' @describeIn FisherZ Correct sampple r with the equation by Olkin & Pratt (1958).
MinVarZ <- function(r, n, k = (-7 + 9*sqrt(2))/2) {
  if(any(n < 5)) stop("Sample size must be greater than 3")
  if(!is.double(k)) stop("Need a real number for 'k'")

  df <- n-1 # if all µ and σ are unkonwn (if µ was know it'd be n)
  #k <- (-7 + 9*sqrt(2))/2 # k = 2.87 hinted by Olkin & Pratt (1958)
  #k <- 3 # official approximation by Olkin & Pratt (1958), form. 2.7
  G <- r * ( 1 + ((1-r^2) / (2 * (df - k))) ) # formula 2.7
  return(G)
}


#' MinVarZ.pr
#' @describeIn FisherZ Correct sample r with the equation by Olkin & Pratt (1958).
#' @details 
MinVarZ.pr <- function(r, n) {
  fc2 <- function(t, r, df) {
    (1-t)^((df-4)/2)  / (1 - t + t*r^2)^0.5 / sqrt(t) #seems fastest
  }
  Fc2 <- function(r, df) integrate(fc2, 0, 1, r = r, df = df)$value

  if(any(n < 5)) stop("Sample size must be greater than 4")
  df <- n-1 # if all µ and σ are unkonwn (if µ was know it'd be n)
  dfh <- df/2
  # Use value for gamma(0.5) to save computation time - https://oeis.org/A002161
  GammaOfHalf <- log(1.772453850905516027298167483341145182797549456122387128213807789852911284591032181374950656738544665)
  # use lgamma tp prevent overflow
  Component1 <- exp(lgamma(dfh-0.5) - lgamma(dfh-1) - GammaOfHalf)
  Component2 <- mapply(Fc2, df = df, r = r)
  G <- r * Component1 * Component2
  
  return(G)
}


# This function throws a lot of warnings
MinVarZ.hg.udf <- function(r, n) {
  if(any(n < 5)) stop("Sample size must be greater than 3")
  
  G <- r * hypergeo(0.5, 0.5, (n-1)/2, 1-r^2)
  return(Re(G))
}





#' MeanR
#' @param R A vector of correlations
#' @param N A vector of the sample sizes each correlation is based on
#' @details The mean of correlations is weighted by the size of 
#' the sample behind each correlation.
#' - Default uses no correction method
#' - Fisher is the Fisher Z correction
#' - Hotelling is Hotelling's correction methog
#' - MinVar is the correction proposed by Olkin & Pratt (1958)
#' - Precise is a correction proposed by Olkin & Pratt (1958)
MeanR <- function( R, N, Method = c("None", "Fisher", "Hotelling", 
                                    "MinVar", "TrueK", "Precise"), ... ) {
  M <- match.arg(Method)
  Res <- switch(M,
           None = MeanR_None(R, N, ...),
           Fisher = MeanR_Fisher(R, N, ...),
           Hotelling = MeanR_Hotelling(R, N, ...),
           MinVar = MeanR_MinVar(R, N, k = 3, ...),
           TrueK = MeanR_MinVar(R, N, ...),
           Precise = MeanR_Precise(R, N, ...),
           MeanR_None(R, N, ...)
           )
  return(Res)
}


#' MeanR_Default
#' @describeIn MeanR Average of correlations without correction
MeanR_None <- function( R, N, na.rm = FALSE ) {
  if(any(N <= 1)) stop("Sample size 'N' must be larger than 1")
  NaCount <- ifelse(na.rm == TRUE, sum(is.na(R)), 0)
  if(length(N) < length(R)) N <- rep(N, length.out = length(R))
    
  Num <- sum((N-1) * R, na.rm)
  Den <- sum(N, na.rm) - (length(N)-NaCount)
  if(any(Den == 0)) Den[Den == 0] <- NA
  return(Num / Den)
}

#' MeanR_Fisher
#' @describeIn MeanR 
MeanR_Fisher <- function( R, N ) {
  R_ <- FisherZ(R)
  Mean <- MeanR_None(R_, N)
  return(FisherZInv(Mean))
}


#' MeanR_Fisher
#' @describeIn MeanR 
#' @references Hotelling H (1953) New light on the correlation 
#' coefficient and its transforms. J R Stat Soc B 15:193–232.
MeanR_Hotelling <- function( R, N ) {
  df <- N-2
  Z <- HotellingZ(R, df) # df missing
  Mean <- MeanR_None(Z, N)
  return(HotellingZInv(Mean, R, df))
}


#' MeanR_MinVar
#' @describeIn MeanR Minimum variance estimator by Olkin & Pratt (1958)
#' @references TODO
MeanR_MinVar <- function( R, N ) {
  Mean <- MeanR_None(MinVarZ(R, N), N)
  return(Mean)
}

#' MeanR_MinVar
#' @describeIn MeanR Minimum variance estimator by Olkin & Pratt (1958)
#' @references TODO
MeanR_Precise <- function( R, N ) {
  Mean <- MeanR_None(MinVarZ.pr(R, N), N)
  return(Mean)
}
