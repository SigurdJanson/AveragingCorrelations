# Code that replicates the study by Corey, Dunlap & Burke
#
#


library(simstudy)
source("./CorrAggBias.R")

N <- seq(10, 50, 10)
R <- seq(0, 0.9, 0.1)
D <- 3:10
M <- c("Default", "Fisher")
NConditions <- length(D) * length(N) * length(R)
NIterations <- 1000


# Create data frame
# EmptyDf <- data.frame(Rho       = factor(numeric(), levels = R),
#                  SampleSize = factor(numeric(), levels = N),
#                  Samples    = factor(numeric(), levels = D),
#                  Method     = factor(character(), levels = M),
#                  Iteration  = factor(numeric(), levels = 1:NIterations),
#                  RObs       = numeric(),
#                  RDelta     = numeric(),
#                  stringsAsFactors = FALSE)


#' @note This function creates its own list for each piece of data.
#' Afterwards this list of lists (each with equal length) will be
#' coerced into a data frame.
RandomSample <- function( Iteration, SampleSize, Rho ) {
  Names <- c("Rho", "SampleSize", "Samples", "Method", "Iteration", 
             "RObs", "RDelta")
  # Data
  d <- max(D) # global variable D 
  # Generate correlated data - one sample for all sample sizes
  SimDat <- genCorData(SampleSize, mu = rep(0, d), sigma = rep(1, d), rho = Rho, corstr = "cs")
  # Estimate correlation matrix
  CorMat <- cor(SimDat[, 2:(d+1)])
  
  Result <- vector(mode = "list", length = length(Names))
  names(Result) <- Names
  for(Samples in D) {
    # Choose subset of the random sample
    SubSample <- sample(1:d, Samples)

    # Use only lower triangle of corr. matrix
    PartialCorMat <- CorMat[SubSample, SubSample]
    Correls <- PartialCorMat[lower.tri(PartialCorMat)]

    # Add data to lists - append some twice for there are two method values
    Result$Rho  <- c(Result$Rho, rep(Rho, 2)) # 
    Result$SampleSize  <- c(Result$SampleSize, rep(SampleSize, 2)) # 
    Result$Samples  <- c(Result$Samples, rep(Samples, 2)) # 
    Result$Iteration  <- c(Result$Iteration, rep(Iteration, 2)) # 

    RObs <- MeanR( Correls, SampleSize, Method = "Fisher")
    Result$Method  <- c(Result$Method, "Default") 
    Result$RObs    <- c(Result$RObs, RObs)
    Result$RDelta  <- c(Result$RDelta, Rho - RObs)
    
    RObs <- MeanR( Correls, SampleSize, Method = "Fisher")
    Result$Method  <- c(Result$Method, "Fisher") 
    Result$RObs    <- c(Result$RObs, RObs)
    Result$RDelta  <- c(Result$RDelta, Rho - RObs)
    
    #RObs <- MeanR( CorMat, n, Method = "Hotelling")
    #...
  }
  
  class(Result) <- "data.frame"
  attr(Result, "row.names") <- .set_row_names(length(Result[[1]]))
  return(Result)
}

TotalList <- vector("list", length(N)*length(R)) #list()
list_of_dfs <- list()
for(n in N) {
  cat(n, "\n")
  for(r in R) {
      list_of_dfs <- lapply(1:NIterations, RandomSample, SampleSize = n, Rho = r)
      TotalList <- append(TotalList, list_of_dfs)
  }
}
cat("Almost done\n")

rm(list_of_dfs)
df <- Reduce(rbind, TotalList)
rm(TotalList)

df$Rho <- as.factor(df$Rho)
df$SampleSize <- as.factor(df$SampleSize)
df$Samples <- as.factor(df$Samples)
df$Method <- as.factor(df$Method)
df$Iteration <- as.factor(df$Iteration)

# Check, if the data is plausible
# see https://www.littlemissdata.com/blog/simple-eda
# see https://cran.r-project.org/web/packages/DataExplorer/vignettes/dataexplorer-intro.html
#summary(df)
#library(DataExplorer)
#plot_bar(df)



# AddLine <- function(x, y, NameD, NameN...) {
#   #hasName(y)
#   NameD <- "" # ignore it - avoid warning
#   Color <- c("#023FA5","#A1A6C8","#E2E2E2","#CA9CA4","#8E063B")
#   Color <- Color[which(paste0("N", N) == NameN)]
#   lines(x, y, col = Color, ...)
# }


#plot(R, R - R, ylim = c(-0.1, 0.1), pch=".")
