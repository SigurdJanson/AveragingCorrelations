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

Result <- array(0, dim = c(length(D), length(N), length(R), length(M)),
                dimnames = list(paste0("D", D), paste0("N", N), 
                                paste0("Rho=", R), M))
#Result <- provideDimnames(Result, sep="", base = list(c("D", "SampleSize", "Rho", "Method")))
#, sep = "_", base = list('row','col','lev')
for(d in D) {
  cat("D = ", d, "\n")
  for(n in N) {
    for(r in R) {
      for(i in 1:NIterations) {
        # Generate correlated data
        SimDat <- genCorData(n, mu = rep(0, d), sigma = rep(1, d), rho = r, corstr = "cs")
        # Estimate correlation matrix
        CorMat <- cor(SimDat[, 2:(d+1)])
        # Average it - use only lower triangle of corr. matrix
        Correls <- CorMat[lower.tri(CorMat)]
        MeanD <- MeanR( Correls, n, Method = "Default")
        MeanF <- MeanR( Correls, n, Method = "Fisher")
        #MeanH <- MeanR( CorMat, n, Method = "Hotelling")
        
        # Sum up
        # IdxD <- which(D == d)
        # IdxN <- which(N == n)
        # IdxR <- which(R == r)
        # Result[IdxD, IdxN, IdxR, 1] <- Result[IdxD, IdxN, IdxR, 1] + MeanD
        # Result[IdxD, IdxN, IdxR, 2] <- Result[IdxD, IdxN, IdxR, 2] + MeanF
        #Result[IdxD, IdxN, IdxR, 3] <- Result[IdxD, IdxN, IdxR, 3] + MeanH
      }
    }
  }
}
Result <- Result / NIterations
#summary(Result)

AddLine <- function(x, y, NameD, NameN...) {
  #hasName(y)
  NameD <- "" # ignore it - avoid warning
  Color <- c("#023FA5","#A1A6C8","#E2E2E2","#CA9CA4","#8E063B")
  Color <- Color[which(paste0("N", N) == NameN)]
  lines(x, y, col = Color, ...)
}

plot(R, R - R, ylim = c(-0.1, 0.1), pch=".")
#RDiff <- apply(Result[,,, "Default"], c(1, 2), function(x){R-x})
RDiff <- apply(Result, c(1, 2, 4), function(x){R-x})
apply(RDiff[,,, "Default"], c(2, 3), AddLine, x=R, NameN = )
apply(RDiff[,,, "Fisher"], c(2, 3), AddLine, x=R)
