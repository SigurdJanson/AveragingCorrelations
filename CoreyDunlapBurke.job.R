# Code that replicates the study by Corey, Dunlap & Burke
#
#

library(simstudy)
source("./CorrAggBias.R")


# Variables ----
N <- seq(10, 50, 10)   # sample size
R <- seq(0.00, 0.95, 0.05)  # Rho
D <- 3:10              # Data sets to correlate
#M <- c("None", "Fisher", "Hotelling", "Hotelling2", "MinVar", "TrueK", "Precise", "Squared")
M <- c("Fisher") # Correction method
DataStruc <- "Indie" # structure of the data, one of c("Matrix", "Indie", )
NConditions <- length(D) * length(N) * length(R) * length(M)
NIterations <- 50000

# Create an empty data frame to allocate the memory (saves time (a lot))
Names <- c("Rho", "SampleSize", "Samples", "Method", "Iteration",
           "RObs", "RDelta")
Result <- vector(mode = "list", length = length(Names))
names(Result) <- Names
for(Column in 1:length(Result)) {
  if(Names[Column] %in% c("Method"))
    Result[[Column]] <- character(NConditions*NIterations) 
  else
    if(Names[Column] %in% c("Rho", "RObs", "RDelta"))
      Result[[Column]] <- double(NConditions*NIterations) 
    else
      Result[[Column]] <- integer(NConditions*NIterations)
}
rm(Column)


# SIMULATION ----
#' Computes averaged correlations for a matrix, i.e. these averaged 
#' coefficients ae dependent from one another.
#' @description To enhance performance this function creates a large correlation
#' matrix and then picks a random subset among those to create an average
#' of d in D correlation coefficients.
#' @note This function uses global variables. It is not intended to be used
#' outside this simulation. The global variables are:
#' - D (the number of data sets to be correlated and averaged) {Read}
#' - Result (where the averaged correlations are stored) {Read/Write}
RandomSample_Matrix <- function(Iteration, SampleSize, Rho) {
  # Data
  MaxD <- max(D) # global variable D
  # Generate correlated data - one sample for all sample sizes
  SimDat <- genCorData(SampleSize, mu = rep(0, MaxD),
                       sigma = rep(1, MaxD), rho = Rho, corstr = "cs")
  # Estimate correlation matrix
  CorMat <- cor( SimDat[, 2:(MaxD+1)] )

  for(Samples in D) {
    # Choose subset of the random sample
    SubSample <- sample(1:MaxD, Samples)

    # Use only lower triangle of corr. matrix
    PartialCorMat <- CorMat[SubSample, SubSample]
    Correls <- PartialCorMat[lower.tri(PartialCorMat)]

    #for(m in M) {
      m <- M
      Result$Rho[Index]        <<- Rho #
      Result$SampleSize[Index] <<- SampleSize #
      Result$Samples[Index]    <<- Samples  #
      Result$Iteration[Index]  <<- Iteration #
      
      RObs <- MeanR(Correls, SampleSize, Method = m)
      Result$Method[Index]  <<- m
      Result$RObs[Index]    <<- RObs
      Result$RDelta[Index]  <<- Rho - RObs
      Index <<- Index +1
    #}
  }#for
}



#' Computes averaged correlations for a matrix, i.e. these averaged 
#' coefficients ae dependent from one another.
#' @description To enhance performance this function creates a large correlation
#' matrix and then picks a random subset among those to create an average
#' of d in D correlation coefficients.
#' @note This function uses global variables. It is not intended to be used
#' outside this simulation. The global variables are:
#' - D (the number of data sets to be correlated and averaged) {Read}
#' - Result (where the averaged correlations are stored) {Read/Write}
RandomSample_Indie <- function(Iteration, SampleSize, Rho) {
  genIndieData <- function(n, mu, sigma = 1, rho) {
    k <- length(mu)
    NValues <- n * k
    Comm <- rnorm(NValues, mu[1], sigma)
    ErrX <- rnorm(NValues, mu[1], sigma)
    ErrY <- rnorm(NValues, mu[1], sigma)
    X <- Comm * sqrt(rho) + ErrX * sqrt(1-rho)
    Y <- Comm * sqrt(rho) + ErrY * sqrt(1-rho)
    X <- split(X, rep(1:k, each = n))
    Y <- split(Y, rep(1:k, each = n))
    R <- mapply(cor, X, Y)
    
  }
  # Data
  MaxD <- max(D) # global variable D
  MaxDR <- MaxD * (MaxD-1) / 2
  # Generate correlated data - one sample for all sample sizes
  AllCorrels <- genIndieData(SampleSize, mu = rep(0, MaxDR), rho = Rho)

  for(Samples in D) {
    # Choose subset of the random sample
    SubSample <- sample(1:MaxDR, Samples)
    Correls   <- AllCorrels[SubSample]
    
    #for(m in M) {
    m <- M
    Result$Rho[Index]        <<- Rho #
    Result$SampleSize[Index] <<- SampleSize #
    Result$Samples[Index]    <<- Samples  #
    Result$Iteration[Index]  <<- Iteration #
    
    RObs <- MeanR(Correls, SampleSize, Method = m)
    Result$Method[Index]  <<- m
    Result$RObs[Index]    <<- RObs
    Result$RDelta[Index]  <<- Rho - RObs
    Index <<- Index +1
    #}
  }#for
}



# Computations ----
if(DataStruc == "Indie")
  RandomSample <- RandomSample_Indie
if(DataStruc == "Matrix")
  RandomSample <- RandomSample_Matrix

Index <- 1 # Index to iterate through data
cat(M, " - ", DataStruc, "\n")
for(n in N) {
  cat(n, "\n")
  for(r in R) {
    lapply(1:NIterations, RandomSample, SampleSize = n, Rho = r)
  }
}
class(Result) <- "data.frame"
attr(Result, "row.names") <- .set_row_names(length(Result[[1]]))

Result$Rho <- as.factor(Result$Rho)
Result$SampleSize <- as.factor(Result$SampleSize)
Result$Samples <- as.factor(Result$Samples)
Result$Method <- as.factor(Result$Method)
Result$Iteration <- as.factor(Result$Iteration)

rm(Names, Index, n, r)


# STORAGE ----
save(Result, file=paste("./data/CoreyDunlapBurke", DataStruc, M, "Raw.Rda", sep = "_"))
