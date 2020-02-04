# Code that replicates the study by Corey, Dunlap & Burke
#
#

library(simstudy)
source("./CorrAggBias.R")


# Create data frame
# EmptyDf <- data.frame(Rho   = factor(numeric(), levels = R),
#                  SampleSize = factor(numeric(), levels = N),
#                  Samples    = factor(numeric(), levels = D),
#                  Method     = factor(character(), levels = M),
#                  Iteration  = factor(numeric(), levels = 1:NIterations),
#                  RObs       = numeric(),
#                  RDelta     = numeric(),
#                  stringsAsFactors = FALSE)


# Variables ----
N <- seq(10, 50, 10)   # sample size
R <- seq(0, 0.9, 0.1)  # Rho
D <- 3:10              # Data sets to correlate
#M <- c("None", "Fisher", "Hotelling", "MinVar", "TrueK", "Precise")
M <- c("MinVar") # Correction method
NConditions <- length(D) * length(N) * length(R) * length(M)
NIterations <- 10000

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

    for(m in M) {
      Result$Rho[Index]        <<- Rho #
      Result$SampleSize[Index] <<- SampleSize #
      Result$Samples[Index]    <<- Samples  #
      Result$Iteration[Index]  <<- Iteration #
      
      RObs <- MeanR(Correls, SampleSize, Method = m)
      Result$Method[Index]  <<- m
      Result$RObs[Index]    <<- RObs
      Result$RDelta[Index]  <<- Rho - RObs
      Index <<- Index +1
    }
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
    
    # Add data to lists - append some twice for there are two method values
    Result$Rho[Index]        <<- Rho #
    Result$SampleSize[Index] <<- SampleSize #
    Result$Samples[Index]    <<- Samples  #
    Result$Iteration[Index]  <<- Iteration #
    
    RObs <- MeanR(Correls, SampleSize, Method = "None")
    Result$Method[Index]  <<- "None"
    Result$RObs[Index]    <<- RObs
    Result$RDelta[Index]  <<- Rho - RObs
    Index <<- Index +1
    
    Result$Rho[Index]        <<- Rho #
    Result$SampleSize[Index] <<- SampleSize #
    Result$Samples[Index]    <<- Samples  #
    Result$Iteration[Index]  <<- Iteration #
    
    RObs <- MeanR(Correls, SampleSize, Method = "Fisher")
    Result$Method[Index]  <<- "Fisher"
    Result$RObs[Index]    <<- RObs
    Result$RDelta[Index]  <<- Rho - RObs
    Index <<- Index +1
    
    #RObs <- MeanR( CorMat, n, Method = "Hotelling")
    #...
  }#for
}



# Computations ----
Index <- 1 # Index to iterate through data

for(n in N) {
  cat(n, "\n")
  for(r in R) {
    lapply(1:NIterations, RandomSample_Matrix, SampleSize = n, Rho = r)
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
save(Result, file=paste("./data/CoreyDunlapBurke", M, "Raw.Rda", sep = "_"))
