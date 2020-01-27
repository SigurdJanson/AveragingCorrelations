
#' FisherZ
#' Fisher z transformation or it's inverse on a vector 
#' of correlation coefficients.
#' @param r A vectors of correlation coefficients
#' @param z A vectors of fisher/hotelling z transformed correlation 
#' coefficients
#' @param df degrees of freedom of the distribution of the 
#' correlation coefficient
#' @value Hotelling z is NaN for df <= 1
#' @author Jan Seifert
#' @references Hotelling, H. (1953). New light on the correlation coefficient and
#' its transformations. Journal of the Royal Statistical Society, 15, p. 193-225.
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


#' MeanR
#' @param R A vector of correlations
#' @param N A vector of the sample sizes each correlation is based on
#' @details The mean of correlations is weighted by the size of 
#' the sample behind each correlation.
MeanR <- function( R, N, Method = c("Default", "Fisher", "Hotelling", "MinVar"), ... ) {
  M <- match.arg(Method)
  Res <- switch(M,
           Default = MeanR_Default(R, N, ...),
           Fisher = MeanR_Fisher(R, N, ...),
           Hotelling = MeanR_Hotelling(R, N, ...),
           MinVar = MeanR_MinVar(R, N, ...),
           MeanR_Default(R, N, ...)
           )
  return(Res)
}

#' MeanR_Default
#' @describeIn MeanR Average of correlations without correction
MeanR_Default <- function( R, N, na.rm = FALSE ) {
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
  Mean <- MeanR_Default(R_, N)
  return(FisherZInv(Mean))
}


#' MeanR_Fisher
#' @describeIn MeanR 
MeanR_Hotelling <- function( R, N ) {
  df <- N-2
  Z <- HotellingZ(R, df) # df missing
  Mean <- MeanR_Default(Z, N)
  return(HotellingZInv(Mean, R, df))
}


#' MeanR_MinVar
#' @describeIn MeanR Minimum variance estimator by Olkin & Pratt (1958)
MeanR_MinVar <- function( R, N ) {
  stop("Not implemented")
}

