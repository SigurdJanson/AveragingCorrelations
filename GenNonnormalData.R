#' GenData
#' @param Data A data set to bootstrap from. Each column represents a 
#' random variable with any desired random distribution. Correlations are 
#' determined between columns.
#' @param NObs Number of observations used to find the correlation matrix.
#' @param CorMatrix Correlation matrix.
#' @param NFactors 
#' @param Max.Trials 
#' @param Initial.Multiplier 
#' @param Seed Seed of the random number generator to replicate results
#' @details 
#' @note   
# All of the commands shown above draw random samples from specified population 
# distributions. As an alternative, one can reproduce distributions without 
# sampling error. For example, working with a supplied data set, one can replace 
# the 2nd line in this block with:
#    \code{Distributions[,i] <- Supplied.Data[,i]}
#' @examples 
#' N <- 5000
#' # Use chi², normal and log-normal distribution
#' sampleData <- cbind(rchisq(N, df = 2), rnorm(N, 0, 1), rlnorm(N, 0, 1))
#' GenData(SampleData, 10)
#' sampleData <- cbind(runif(N, 0, 1), rexp(N, rate = 1), rpois(N, lambda = 4))
#' GenData(SampleData, 10)
#' @author John Ruscio and Walter Kaczetow; refactored by Jan Seifert
#' @references Ruscio, J., & Kaczetow, W. (2008). Simulating Multivariate Nonnormal Data Using an Iterative Algorithm. 
#' Multivariate Behavioral Research, 43 (3), 355–381. 
GenData <- function(ProbDistr, NObs, CorMatrix, NFactors = 0, 
                    Max.Trials = 5, Initial.Multiplier = 1, Seed = NA)
{
  # Initialize variables and (if applicable) set random number seed (step 1)
  k <- ncol(ProbDistr) #-dim(Data)[2]# Number of variables to generate
  Data <- matrix(0, nrow = NObs, ncol = k) # Matrix to store the simulated data
  Iteration <- 0                     # Iteration counter
  Best.RMSR <- 1                     # Lowest RMSR correlation
  Trials.Without.Improvement <- 0    # Trial counter
  if (!is.na(Seed)) set.seed(Seed)   # If user specified a seed
  
  # Generate distribution for each variable (step 2) 
  Distributions <- matrix(0, nrow = NObs, ncol = k) # Matrix to store each variable’s score distribution
  for (i in 1:k)
    Distributions[,i] <- sort(sample(ProbDistr[,i], NObs, replace = TRUE)) 

  # Calculate and store a copy of the target correlation matrix (step 3)
  #-CorMatrix <- CorMatrix #-cor(Data)
  Intermediate.Corr <- CorMatrix
  
  # If number of latent factors was not specified, determine it through 
  # parallel analysis (step 4)
  if (NFactors == 0) {
    Eigenvalues.Observed <- eigen(Intermediate.Corr)$values
    Eigenvalues.Random <- matrix(0, nrow = 100, ncol = k)
    Random.Data <- matrix(0, nrow = NObs, ncol = k)
    for (i in 1:100) {
      for (j in 1:k) {
        Random.Data[,j] <- sample(Distributions[,j], size = NObs, replace = TRUE)
      }
      Eigenvalues.Random[i,] <- eigen(cor(Random.Data))$values
    }
    # calculate mean eigenvalue for each factor
    Eigenvalues.Random <- apply(Eigenvalues.Random, 2, mean) 
    NFactors <- max(1, sum(Eigenvalues.Observed > Eigenvalues.Random))
  }
  
  # Generate random normal data for shared and unique components, 
  # initialize factor loadings (steps 5, 6) 
  Shared.Comp <- matrix(rnorm(NObs * NFactors, 0, 1), nrow = NObs, ncol = NFactors)
  Unique.Comp <- matrix(rnorm(NObs * k, 0, 1), nrow = NObs, ncol = k)
  Shared.Load <- matrix(0, nrow = k, ncol = NFactors)
  Unique.Load <- matrix(0, nrow = k, ncol = 1)
  
  # Begin loop that ends when specified number of iterations pass without 
  # improvement in RMSR correlation 
  while (Trials.Without.Improvement < Max.Trials) {
    Iteration <- Iteration + 1
    # Calculate factor loadings and apply to reproduce desired correlations (steps 7, 8) 
    Fact.Anal <- Factor.Analysis(Intermediate.Corr, CorrMatrix = TRUE, NFactors = NFactors)
    if (NFactors == 1) Shared.Load[,1] <- Fact.Anal$loadings
    else Shared.Load <- Fact.Anal$loadings
    Shared.Load[Shared.Load > 1] <- 1
    Shared.Load[Shared.Load < -1] <- -1
    if (Shared.Load[1,1] < 0) Shared.Load <- Shared.Load * -1
    Shared.Load.sq <- Shared.Load * Shared.Load
    for (i in 1:k)
      if (sum(Shared.Load.sq[i,]) < 1) Unique.Load[i,1] <- (1 - sum(Shared.Load.sq[i,]))
    else Unique.Load[i,1] <- 0
    Unique.Load <- sqrt(Unique.Load)
    for (i in 1:k) {
      # %*% operator = matrix multiplication, t() = transpose
      Data[,i] <- (Shared.Comp %*% t(Shared.Load))[,i] + Unique.Comp[,i] * Unique.Load[i,1]
    }

    # Replace normal with nonnormal distributions (step 9) 
    for (i in 1:k) {
      Data <- Data[sort.list(Data[,i]),]
      Data[,i] <- Distributions[,i]
    }
    
    # Calculate RMSR correlation, compare to lowest value, take appropriate action (steps 10, 11, 12)
    Reproduced.Corr <- cor(Data)
    Residual.Corr <- CorMatrix - Reproduced.Corr
    RMSR <- sqrt(sum(Residual.Corr[lower.tri(Residual.Corr)] * 
                       Residual.Corr[lower.tri(Residual.Corr)]) /
                   (.5 * (k * k - k)))
    if (RMSR < Best.RMSR) {
      Best.RMSR <- RMSR
      Best.Corr <- Intermediate.Corr
      Best.Res <- Residual.Corr
      Intermediate.Corr <- Intermediate.Corr + Initial.Multiplier * Residual.Corr
      Trials.Without.Improvement <- 0
    } else {
      Trials.Without.Improvement <- Trials.Without.Improvement + 1
      Current.Multiplier <- Initial.Multiplier * .5 ^ Trials.Without.Improvement
      Intermediate.Corr <- Best.Corr + Current.Multiplier * Best.Res
    }
  } # end of the while loop
  
  # Construct the data set with the lowest RMSR correlation (step 13) 
  Fact.Anal <- Factor.Analysis(Best.Corr, Corr.Matrix = TRUE, N.Factors = NFactors)
  if (NFactors == 1) {
    Shared.Load[,1] <- Fact.Anal$loadings
  } else {
    Shared.Load <- Fact.Anal$loadings
  } 
  Shared.Load[Shared.Load > 1] <- 1
  Shared.Load[Shared.Load < -1] <- -1
  if (Shared.Load[1,1] < 0) Shared.Load <- Shared.Load * -1
  Shared.Load.sq <- Shared.Load * Shared.Load
  
  for (i in 1:k)
    if (sum(Shared.Load.sq[i,]) < 1) Unique.Load[i,1] <- (1 - sum(Shared.Load.sq[i,]))
  else Unique.Load[i,1] <- 0
  Unique.Load <- sqrt(Unique.Load)
  for (i in 1:k)
    Data[,i] <- (Shared.Comp %*% t(Shared.Load))[,i] + Unique.Comp[,i] * Unique.Load[i,1]
  Data <- apply(Data, 2, scale) # standardizes each variable in the matrix
  for (i in 1:k) {
    Data <- Data[sort.list(Data[,i]),]
    Data[,i] <- Distributions[,i]
  }
  
  # Report the results and return the simulated data set (step 14) 
  Iteration <- Iteration - Max.Trials
  cat("\nN =", NObs, ", k =", k, ",", Iteration, "Iterations,",
      NFactors, "Factors, RMSR r =", round(Best.RMSR,3), "\n")
  return(Data)
}




#' Factor.Analysis
#' Home made factor analysis as support function for 'GenData'.
#' @param Data Raw data or correlation matrix. Will be coerced to a matrix. Will be 
#' interpreted as correlation matrix only if 'Corr.Matrix' is TRUE.
#' @param CorrMatrix Specifies if Data is raw data or a correlation matrix.
#' @param MaxIter Maximum number of iterations for convergence.
#' @param NFactors Number of factors to extract, default is \code{ncol(Data)}
#' @details If 'NFactors' equals the number of variables (i.e. columns of 'Data') 
#' then this function returns the factors with eigenvalues larger than 1.
#' Otherwise it returns the 'NFactors' factors with the highest eigenvalues.
#' @author John Ruscio and Walter Kaczetow; refactored by Jan Seifert
#' @references Ruscio, J., & Kaczetow, W. (2008). Simulating Multivariate Nonnormal 
#' Data Using an Iterative Algorithm. Multivariate Behavioral Research, 43 (3), 355–381. 
Factor.Analysis <- function(Data, CorrMatrix = FALSE, MaxIter = 50, NFactors = NA) {
  Data <- as.matrix(Data)
  k <- ncol(Data) # dim(Data)[2]
  if (is.na(NFactors)) NFactors <- k
  if (!CorrMatrix)
    Cor.Matrix <- cor(Data)
  else 
    Cor.Matrix <- Data

  Criterion <- 0.001
  Old.H2 <- rep(99, k)
  H2 <- rep(0, k)
  Change <- 1
  Iter <- 0
  Factor.Loadings <- matrix(nrow = k, ncol = NFactors)

  while ((Change >= Criterion) & (Iter < MaxIter)) {
    Iter <- Iter + 1
    Eig <- eigen(Cor.Matrix)
    L <- sqrt(Eig$values[1:NFactors])
    for (i in 1:NFactors)
      Factor.Loadings[,i] <- Eig$vectors[,i] * L[i]
    for (i in 1:k)
      H2[i] <- sum(Factor.Loadings[i,] * Factor.Loadings[i,])
    Change <- max(abs(Old.H2 - H2))
    Old.H2 <- H2
    diag(Cor.Matrix) <- H2
  }
  
  if (NFactors == k) NFactors <- sum(Eig$values > 1)
  return(list(loadings = Factor.Loadings[,1:NFactors], factors = NFactors))
}

